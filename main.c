/*
 * Obbe Tuinenburg
 * Recycling model
*/
#include <omp.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <time.h>
#include <netcdf.h>
#include "recycling.h"
#include "defines.h"

int main(int argc, char *argv[]) {

	// Get the starting date and simulation length from command line arguments
	int year=atoi(argv[1]);
	int month=atoi(argv[2]);
	int day=atoi(argv[3]);
	int release_days=atoi(argv[4]); // Number of days to emit the tracked moisture
	int nday=release_days+24; // Number of days to run the simulation

	struct simulation_settings* allsims= NULL;
    	struct simulation_settings* initial = NULL;

	FILE *ifp;

	int simulationtype,interpolationtype,nparcelsreleased;
	char inputfilepath[100];
	char outputfilepath[100];

    
	// date utils
	int month_days[] = {0, 31, 28, 31, 30, 31, 30, 31 ,31 ,30, 31, 30, 31};
	if(year%4 == 0 && year%100 != 0 || year%400 == 0)
    	{
        	month_days[2] = 29;
	} else {
        	month_days[2] = 28;
	}

	// Looping integers
	int i,j,k,l;

	// OpenMP variables
	int iam = 0 ,np = 1;

	// Allocate memory for ERA5 lsm
	short* lsm = malloc(sizeof(short)*721*1440);
	double* lsmadd = malloc(sizeof(double));
	double* lsmscale = malloc(sizeof(double));
	get_lsm(lsm,lsmadd,lsmscale);
	int num_locations=0;
	float curlatrad,degreelength;
	l=0;

	ifp = fopen("list.txt", "r");

	    if (ifp == NULL)
	    {
	        fprintf(stderr, "Can't open input list.txt\n");
		exit(1);
	    }

    while (fscanf(ifp,"%d %d %d %s %s",
            &simulationtype, &interpolationtype, &nparcelsreleased, inputfilepath,outputfilepath)
            != EOF)
    {
	        if (initial == NULL)
        {
            allsims = initial =  malloc(sizeof(struct simulation_settings));
	}
        else
        {
            allsims->next = malloc(sizeof(struct simulation_settings));
            allsims = allsims->next;
        }
            allsims->simulation_number=simulationtype;
	    if (simulationtype==4){allsims->mixing=6.0;}
	    if (simulationtype==8){allsims->mixing=1.0;}
	    if (simulationtype==9){allsims->mixing=24.0;}
	    if (simulationtype==10){allsims->mixing=120.0;}
	    if (simulationtype==11){allsims->mixing=1.0;}
	    if (simulationtype==12){allsims->mixing=6.0;}
	    if (simulationtype==13){allsims->mixing=24.0;}
	    if (simulationtype==14){allsims->mixing=120.0;}
            allsims->interpolation=interpolationtype;
	    allsims->release_parcel=nparcelsreleased;
	    allsims->runtime=0;
            strcpy(allsims->inputfile,inputfilepath);
            strcpy(allsims->outputfile,outputfilepath);

	    // Allocate memory for input
  	    allsims->release_input = malloc(sizeof(int)*721*1440);
	    get_input(allsims->release_input,inputfilepath);
	
	    // Calculate the number of parcels to release during the simulation
	    num_locations=0;
		for (i=0;i<721*1440;i++){
			if (allsims->release_input[i]==1){
				num_locations++;
			}
		}
		allsims->parcels_in_system=malloc(sizeof(int));
		*allsims->parcels_in_system=0;

		allsims->release_days=release_days;
		allsims->nday=nday;
		allsims->num_release=num_locations*allsims->release_days*24*allsims->release_parcel;
		allsims->grid_size_NW=malloc(sizeof(float)*721);
		allsims->grid_size_EW=malloc(sizeof(float)*721);
		for (i=0;i<721;i++){
			curlatrad=(90-i*0.25)*2*3.14159/360;
			degreelength=((111412.84*cos(curlatrad))-93.5*cos(3*curlatrad)+0.118*cos(5*curlatrad))/1000;
			allsims->grid_size_EW[i]=degreelength;
			degreelength=(111132.92+(-559.82*cos(2*curlatrad))+1.175*cos(4*curlatrad)-0.0023*cos(6*curlatrad))/1000;
			allsims->grid_size_NW[i]=degreelength;}

		// Allocate memory for output
		if (allsims->simulation_number==7){allsims->insystem = malloc(sizeof(double)*25*721*1440);
		memset (allsims->insystem, 0, sizeof (double) * 25*721*1440);} // 3D state for Eulerian
		if (allsims->simulation_number==7||allsims->simulation_number==6||allsims->simulation_number==15){
			allsims->budget= malloc(sizeof(double)*25*721*1440);}
		if (allsims->simulation_number==15){allsims->insystem = malloc(sizeof(double)*2*721*1440);
		memset (allsims->insystem, 0, sizeof (double) * 2*721*1440);} // 3D state for Eulerian

		allsims->insystem_out = malloc(sizeof(double)*721*1440);
		allsims->allocated_out = malloc(sizeof(double)*721*1440);
		allsims->emitted_out = malloc(sizeof(double)*721*1440);
		for (i=0;i<721*1440;i++){
			allsims->insystem_out[i]=0;
			allsims->allocated_out[i]=0;
			allsims->emitted_out[i]=0;
		}

		// Allocated memory for parcels during tracking
		allsims->particle_lat=malloc(sizeof(float)*allsims->num_release);
		allsims->particle_lon=malloc(sizeof(float)*allsims->num_release);
		allsims->particle_lev=malloc(sizeof(float)*allsims->num_release);
		allsims->particle_theta=malloc(sizeof(float)*allsims->num_release);
		allsims->particle_original=malloc(sizeof(float)*allsims->num_release);
		allsims->particle_present=malloc(sizeof(float)*allsims->num_release);
		allsims->particle_curtime=malloc(sizeof(float)*allsims->num_release);

    }

	// Seed the randomizer.//   
	unsigned int iseed = (unsigned int)time(NULL);
	srand (iseed);
   
	// Load data for first two hours from netcdf files
	struct meteoday meteo[2];

	int curmonth=month;
	int curday=day;
	int curyear=year;
	int curhour=0;
	clock_t begin = clock();
	clock_t end = clock();

	// Allocated memory for two hours of data
	meteo[0].q=malloc(sizeof(float) * 25* 721* 1440);
	meteo[0].u=malloc(sizeof(float) * 25* 721* 1440);
	meteo[0].v=malloc(sizeof(float) * 25* 721* 1440);
	meteo[0].w=malloc(sizeof(float) * 25* 721* 1440);
	meteo[0].t=malloc(sizeof(float) * 25* 721* 1440);
	meteo[0].tadd=malloc(sizeof(double));
	meteo[0].tscale=malloc(sizeof(double));
	meteo[0].wadd=malloc(sizeof(double));
	meteo[0].wscale=malloc(sizeof(double));
	meteo[0].vadd=malloc(sizeof(double));
	meteo[0].vscale=malloc(sizeof(double));
	meteo[0].uadd=malloc(sizeof(double));
	meteo[0].uscale=malloc(sizeof(double));
	meteo[0].qadd=malloc(sizeof(double));
	meteo[0].qscale=malloc(sizeof(double));
	meteo[0].pw=malloc(sizeof(float) * 721* 1440);
	meteo[0].nf=malloc(sizeof(float) * 721* 1440);
	meteo[0].ef=malloc(sizeof(float) * 721* 1440);
	meteo[0].E=malloc(sizeof(float) * 721* 1440);
	meteo[0].P=malloc(sizeof(float) * 721* 1440);
	meteo[0].pwadd=malloc(sizeof(double));
	meteo[0].pwscale=malloc(sizeof(double));
	meteo[0].nfadd=malloc(sizeof(double));
	meteo[0].nfscale=malloc(sizeof(double));
	meteo[0].efadd=malloc(sizeof(double));
	meteo[0].efscale=malloc(sizeof(double));
	meteo[0].Eadd=malloc(sizeof(double));
	meteo[0].Escale=malloc(sizeof(double));
	meteo[0].Padd=malloc(sizeof(double));
	meteo[0].Pscale=malloc(sizeof(double));
	meteo[0].lat=malloc(sizeof(float) * 721);
	meteo[0].lon=malloc(sizeof(float) * 1440);
	meteo[1].q=malloc(sizeof(float) * 25* 721* 1440);
	meteo[1].u=malloc(sizeof(float) * 25* 721* 1440);
	meteo[1].v=malloc(sizeof(float) * 25* 721* 1440);
	meteo[1].w=malloc(sizeof(float) * 25* 721* 1440);
	meteo[1].t=malloc(sizeof(float) * 25* 721* 1440);
	meteo[1].tadd=malloc(sizeof(double));
	meteo[1].tscale=malloc(sizeof(double));
	meteo[1].wadd=malloc(sizeof(double));
	meteo[1].wscale=malloc(sizeof(double));
	meteo[1].vadd=malloc(sizeof(double));
	meteo[1].vscale=malloc(sizeof(double));
	meteo[1].uadd=malloc(sizeof(double));
	meteo[1].uscale=malloc(sizeof(double));
	meteo[1].qadd=malloc(sizeof(double));
	meteo[1].qscale=malloc(sizeof(double));
	meteo[1].pw=malloc(sizeof(float) * 721* 1440);
	meteo[1].nf=malloc(sizeof(float) * 721* 1440);
	meteo[1].ef=malloc(sizeof(float) * 721* 1440);
	meteo[1].E=malloc(sizeof(float) * 721* 1440);
	meteo[1].P=malloc(sizeof(float) * 721* 1440);
	meteo[1].pwadd=malloc(sizeof(double));
	meteo[1].pwscale=malloc(sizeof(double));
	meteo[1].nfadd=malloc(sizeof(double));
	meteo[1].nfscale=malloc(sizeof(double));
	meteo[1].efadd=malloc(sizeof(double));
	meteo[1].efscale=malloc(sizeof(double));
	meteo[1].Eadd=malloc(sizeof(double));
	meteo[1].Escale=malloc(sizeof(double));
	meteo[1].Padd=malloc(sizeof(double));
	meteo[1].Pscale=malloc(sizeof(double));
	meteo[1].lat=malloc(sizeof(float) * 721);
	meteo[1].lon=malloc(sizeof(float) * 1440);
	
	
	// Load initial two hours of forcing data into memory
	struct meteoday *p;
	i=load_meteo_hour(curyear,curmonth,curday,curhour,&(meteo[0]));
	i=load_meteo_hour(curyear,curmonth,curday,curhour,&(meteo[1]));

	// ERA5 Latitude and Longitude values
	float* lats=malloc(sizeof(float)*721);
	float* lons=malloc(sizeof(float)*1440);
	for (i=0;i<721;i++){lats[i]=90-0.25*i;}
	for (i=0;i<1440;i++){lons[i]=0.25*i;}


	int daynum=0;
	while (daynum<nday-1){
		// Free array for low time point 
		meteo[0]=meteo[1];

		// Load new day
		curhour++;
		if (curhour>23){
			curhour=0;
			curday++;
			daynum++;
		if (curday>month_days[curmonth]){
			curday=1;
			curmonth++;
			if (curmonth>12){
				curmonth=1;
				curyear++;}}}
		i=load_meteo_hour(curyear,curmonth,curday,curhour,&meteo[1]);
		printf("%d %d %d %d\n",curyear,curmonth,curday,curhour);
        allsims = initial;
        while (allsims)
        {
		begin = clock();

		// For all simulations, do simulation for one hour
		allsims->daynum=daynum;
		allsims->hour=curhour;
		if (allsims->simulation_number<6||(allsims->simulation_number>7&&allsims->simulation_number<15)){
		lagrangian_simulation(meteo,*allsims);}
		if (allsims->simulation_number>5&&allsims->simulation_number<8){
		eulerian_simulation(meteo,*allsims);}
		end = clock();
		allsims->runtime+=(float)(end - begin) / CLOCKS_PER_SEC;
		allsims=allsims->next;

	}

	}


// Output
// Write output arrays to netcdf file
allsims = initial;
while (allsims)
        {
	i=write_output(allsims->emitted_out,allsims->insystem_out,allsims->allocated_out,lats,lons,allsims->nday,allsims->outputfile,allsims->runtime);
	allsims=allsims->next;
	}


   return 0;
}
