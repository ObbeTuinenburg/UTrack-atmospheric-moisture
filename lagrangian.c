#include <omp.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <time.h>
#include "netcdf.h"
#include "recycling.h"
#include "defines.h"

int lagrangian_simulation(struct meteoday* meteo,struct simulation_settings simulation){
	// Simulates the moisture tracking particles for one day using the forcing in meteo and simulation parameters in simulation
	float mean_present,mean_released;
	int daynum;
	int i,j,k,l,iam,np;
	int parcels_released;
	for (daynum=simulation.hour;daynum<simulation.hour+1;daynum++){
		// Reset daily counter for all particles
		for (i=0;i<*simulation.parcels_in_system;i++){
			if (simulation.particle_curtime[i]>23.999){simulation.particle_curtime[i]=0;}} 

		mean_present=0;
		mean_released=0;
		for (i=0;i<*simulation.parcels_in_system;i++){
			mean_present+=simulation.particle_present[i];
			mean_released+=simulation.particle_original[i];}
		if (simulation.daynum>0&&simulation.hour==0){
		printf("Day %d, hour %d: number of parcels: %i, mean fraction present: %f\n",simulation.daynum,simulation.hour,*simulation.parcels_in_system,mean_present/mean_released);}

		// Create new parcels if still releasing
		if (simulation.daynum<simulation.release_days){
			for (i=0;i<721;i++){
				for (j=0;j<1440;j++){
					if (simulation.release_input[i*1440+j]==1){
						for (k=simulation.hour;k<simulation.hour+1;k++){
							parcels_released=(int)ceil(fabs((meteo[0].E[i*1440+j]*(*meteo[0].Escale)+(*meteo[0].Eadd))*1000.0)*simulation.release_parcel);
							for (l=0;l<parcels_released;l++){
		simulation.particle_lat[*simulation.parcels_in_system]=90-i*0.25-0.25*(double)rand()/(RAND_MAX);
		simulation.particle_lon[*simulation.parcels_in_system]=j*0.25+0.125*(double)rand()/(RAND_MAX);
		simulation.particle_original[*simulation.parcels_in_system]=fabs((meteo[0].E[i*1440+j]*(*meteo[0].Escale)+(*meteo[0].Eadd))*1000.0)/parcels_released;
		simulation.particle_lev[*simulation.parcels_in_system]=get_starting_height(meteo,simulation,i,j);
		if (simulation.simulation_number==3){
			simulation.particle_theta[*simulation.parcels_in_system]=get_potential_temp(meteo,i,j,simulation.particle_lev[*simulation.parcels_in_system]);
		}
		simulation.particle_present[*simulation.parcels_in_system]=simulation.particle_original[*simulation.parcels_in_system];
		simulation.particle_curtime[*simulation.parcels_in_system]=(float)k;
		simulation.emitted_out[0*721*1440+i*1440+j]+=simulation.particle_original[*simulation.parcels_in_system];
		(*simulation.parcels_in_system)++;
							}
						}
					}
				}
			}
		}
		
		// Track the particles to the end of the day in an OpenMP loop
		#pragma omp parallel default(shared) private(iam, np,i)
		{
		#if defined (_OPENMP)
		np=omp_get_num_threads();
		iam=omp_get_thread_num();
		#endif
		
		#pragma omp for schedule(dynamic,1) nowait
		for (i=0;i<*simulation.parcels_in_system;i+=1){
			trajectory(meteo,simulation,i);
		}
		}
	}
		return 0;

	}

float get_potential_temp(struct meteoday* meteo,int i,int j,float lev){
	int curlevidx;
	if (lev<750){
	curlevidx=((int)floor((lev)/50.0));}
	if (lev>=750){
	curlevidx=((int)floor((lev-750)/25.0))+15;}
	return (meteo[0].t[curlevidx*1440*721+i*721+j]*(*meteo[0].tscale)+(*meteo[0].tadd))*pow((1000.0/lev),0.286);}


float get_starting_height(struct meteoday* meteo,struct simulation_settings simulation,int lat_idx,int lon_idx){
	float qsum=0;
	float fraction,qlev;
	float store=0;
	int i;
	for (i=0;i<25;i++){
		qsum+=(meteo[0].q[i*721*1440+lat_idx*1440+lon_idx]*(*meteo[0].qscale)+(*meteo[0].qadd));}
	fraction=qsum*(float)rand()/RAND_MAX;
	for (i=0;i<25;i++){
		qlev=(meteo[0].q[(24-i)*721*1440+lat_idx*1440+lon_idx]*(*meteo[0].qscale)+(*meteo[0].qadd));
		if (store>fraction){
			if (i<11){return 1000-i*25.0;}
			if (i>10){return 750-(i-10)*50.0;}
		}

		store+=qlev;}
	return 500.0;
	}




