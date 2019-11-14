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

int eulerian_simulation(struct meteoday* meteo,struct simulation_settings simulation){
	// Simulates the moisture tracking particles for one day using the forcing in meteo and simulation parameters in simulation
	float mean_present,mean_released;
	int daynum=simulation.daynum;
	float time;
	float timefrac;
	float dt=0.1;
	int i,j,k,l,iam,np;
	float NW,EW,PW,P,UW;
	float NWflux,EWflux,UWflux;
	float budgetsum=0;
	float statemax=0;
	int imax;
	float qloc;
	float qsum,qsum2;


	// Clear the budget matrix
	for (i=0;i<25*721*1440;i++){simulation.budget[i]=0;}

	for (time=simulation.hour;time<simulation.hour+1;time+=dt){
		k=(int)floor((time-floor(time))*24);
		timefrac=(time-floor(time))*24-k; // Fraction of the current hour

		// Add moisture to the atmosphere if still releasing
		if (simulation.daynum<simulation.release_days){
			#pragma omp parallel default(shared) private(i,j,l,qsum,qsum2)
			{
			#pragma omp for schedule(dynamic,20) nowait
			for (i=0;i<721;i++){
				for (j=0;j<1440;j++){
					if (simulation.simulation_number==7){
					qsum=0;
					qsum2=0;
					for (l=0;l<25;l++){
						qsum2+=simulation.insystem[l*1440*721+i*1440+j];
						qsum+=(meteo[0].q[l*1440*721+i*1440+j]*(*meteo[0].qscale)+(*meteo[0].qadd));}
					for (l=0;l<25;l++){
						simulation.insystem[l*1440*721+i*1440+j]=qsum2*(meteo[0].q[l*1440*721+i*1440+j]*(*meteo[0].qscale)+(*meteo[0].qadd))/qsum;}
					}

					if (simulation.release_input[i*1440+j]==1){
		if (simulation.simulation_number==7){
		qsum=0;
		for (l=0;l<25;l++){
			if ((meteo[0].q[l*1440*721+i*1440+j]*(*meteo[0].qscale)+(*meteo[0].qadd))<1){
			qsum+=(meteo[0].q[l*1440*721+i*1440+j]*(*meteo[0].qscale)+(*meteo[0].qadd));}}

		for (l=0;l<25;l++){
			if ((meteo[0].q[l*1440*721+i*1440+j]*(*meteo[0].qscale)+(*meteo[0].qadd))<1){
		simulation.insystem[l*1440*721+i*1440+j]+=((meteo[0].q[l*1440*721+i*1440+j]*(*meteo[0].qscale)+(*meteo[0].qadd))/qsum)*fabs((meteo[0].E[i*1440+j]*(*meteo[0].Escale)+(*meteo[0].Eadd))*1000.0)*dt;}
		}
		
		}


		if (simulation.simulation_number==6){
		simulation.insystem_out[0*721*1440+i*1440+j]+=fabs((meteo[0].E[i*1440+j]*(*meteo[0].Escale)+(*meteo[0].Eadd))*1000.0)*dt;}
		simulation.emitted_out[0*721*1440+i*1440+j]+=fabs((meteo[0].E[i*1440+j]*(*meteo[0].Escale)+(*meteo[0].Eadd))*1000.0)*dt;
					}
				}
			}
			}
		}

		// Move moisture around in atmosphere
		// Clear budget matrix
		#pragma omp parallel default(shared) private(i)
		{
		#pragma omp for schedule(dynamic,20) nowait
		for (i=0;i<25*721*1440;i++){simulation.budget[i]=0;}
		}
		
		if (simulation.simulation_number==6){
		// Fill budget matrix
		#pragma omp parallel default(shared) private(i,j,NW,EW,PW,NWflux,EWflux)
		{
		#pragma omp for schedule(dynamic,20) nowait
		for (i=1;i<720;i++){
			for (j=0;j<1440;j++){
				if (simulation.insystem_out[0*721*1440+i*1440+j]>0){
				if (simulation.interpolation==0){
					NW=meteo[0].nf[i*1440+j]*(*meteo[0].nfscale)+(*meteo[0].nfadd);
					PW=meteo[0].pw[i*1440+j]*(*meteo[0].pwscale)+(*meteo[0].pwadd);}
				if (simulation.interpolation==1){
					NW=(1-timefrac)*(meteo[0].nf[i*1440+j]*(*meteo[0].nfscale)+(*meteo[0].nfadd))+timefrac*(meteo[1].nf[i*1440+j]*(*meteo[1].nfscale)+(*meteo[1].nfadd));
					PW=(1-timefrac)*(meteo[0].pw[i*1440+j]*(*meteo[0].pwscale)+(*meteo[0].pwadd))+timefrac*(meteo[1].pw[i*1440+j]*(*meteo[1].pwscale)+(*meteo[1].pwadd));}
				if (PW==0){PW=100;}
				NWflux=(3.6*(NW/PW)/simulation.grid_size_NW[i]*dt);
				if (NWflux>0.5){NWflux=0.5;}
				if (NWflux<-0.5){NWflux=-0.5;}
				if (NWflux>0){
					simulation.budget[(i-1)*1440+j]+=NWflux*simulation.insystem_out[0*721*1440+i*1440+j];
					simulation.budget[i*1440+j]-=NWflux*simulation.insystem_out[0*721*1440+i*1440+j];}
				if (NWflux<0){
					simulation.budget[(i+1)*1440+j]+=fabs(NWflux)*simulation.insystem_out[0*721*1440+i*1440+j];
					simulation.budget[i*1440+j]-=fabs(NWflux)*simulation.insystem_out[0*721*1440+i*1440+j];}
				}
			}
		}
		}
		#pragma omp parallel default(shared) private(i,j,NW,EW,PW,NWflux,EWflux)
		{
		#pragma omp for schedule(dynamic,20) nowait
		for (i=1;i<720;i++){
			for (j=0;j<1440;j++){
				if (simulation.insystem_out[0*721*1440+i*1440+j]>0){
				if (simulation.interpolation==0){
					EW=meteo[0].ef[i*1440+j]*(*meteo[0].efscale)+(*meteo[0].efadd);
					PW=meteo[0].pw[i*1440+j]*(*meteo[0].pwscale)+(*meteo[0].pwadd);}
				if (simulation.interpolation==1){
					EW=(1-timefrac)*(meteo[0].ef[i*1440+j]*(*meteo[0].efscale)+(*meteo[0].efadd))+timefrac*(meteo[1].ef[i*1440+j]*(*meteo[1].efscale)+(*meteo[1].efadd));
					PW=(1-timefrac)*(meteo[0].pw[i*1440+j]*(*meteo[0].pwscale)+(*meteo[0].pwadd))+timefrac*(meteo[1].pw[i*1440+j]*(*meteo[1].pwscale)+(*meteo[1].pwadd));}
				EWflux=(3.6*(EW/PW)/simulation.grid_size_EW[i]*dt);
				if (EWflux>0.5){EWflux=0.5;}
				if (EWflux<-0.5){EWflux=-0.5;}
				if (EWflux>0){
					simulation.budget[i*1440+j+1]+=EWflux*simulation.insystem_out[0*721*1440+i*1440+j];
					simulation.budget[i*1440+j]-=EWflux*simulation.insystem_out[0*721*1440+i*1440+j];
					if (j==1439){
						simulation.budget[i*1440+j+1]-=EWflux*simulation.insystem_out[0*721*1440+i*1440+j];
						simulation.budget[i*1440+j]+=EWflux*simulation.insystem_out[0*721*1440+i*1440+j];
						simulation.budget[i*1440]+=EWflux*simulation.insystem_out[0*721*1440+i*1440+j];
						simulation.budget[i*1440+j]-=EWflux*simulation.insystem_out[0*721*1440+i*1440+j];
					}
				}
				if (EWflux<0){
					simulation.budget[i*1440+j-1]+=fabs(EWflux)*simulation.insystem_out[0*721*1440+i*1440+j];
					simulation.budget[i*1440+j]-=fabs(EWflux)*simulation.insystem_out[0*721*1440+i*1440+j];
					if (j==0){
						simulation.budget[i*1440+j-1]-=fabs(EWflux)*simulation.insystem_out[0*721*1440+i*1440+j];
						simulation.budget[i*1440+j]+=fabs(EWflux)*simulation.insystem_out[0*721*1440+i*1440+j];
						simulation.budget[i*1440+1439]+=fabs(EWflux)*simulation.insystem_out[0*721*1440+i*1440+j];
						simulation.budget[i*1440+j]-=fabs(EWflux)*simulation.insystem_out[0*721*1440+i*1440+j];
					}
				}
				}}}
		}
		}
		if (simulation.simulation_number==7){
		// Fill budget matrix
		#pragma omp parallel default(shared) private(l,i,j,NW,EW,PW,UW,NWflux,EWflux)
		{
		#pragma omp for schedule(dynamic,20) nowait
		for (i=1;i<720;i++){
			for (j=0;j<1440;j++){
				for (l=0;l<25;l++){
					if (simulation.insystem[l*721*1440+i*1440+j]>0){
				if (simulation.interpolation==0){
				NW=meteo[0].v[l*721*1440+i*1440+j]*(*meteo[0].vscale)+(*meteo[0].vadd);}
				if (simulation.interpolation==1){
					NW=(1-timefrac)*(meteo[0].v[l*721*1440+i*1440+j]*(*meteo[0].vscale)+(*meteo[0].vadd))+timefrac*(meteo[1].v[(l)*721*1440+i*1440+j]*(*meteo[1].vscale)+(*meteo[1].vadd));}
				NWflux=(3.6*NW/simulation.grid_size_NW[i]*dt);
				if (NWflux>0.5){NWflux=0.5;}
				if (NWflux<-0.5){NWflux=-0.5;}
				if (NWflux>0){
					simulation.budget[l*721*1440+(i-1)*1440+j]+=NWflux*simulation.insystem[l*721*1440+i*1440+j];
					simulation.budget[l*721*1440+i*1440+j]-=NWflux*simulation.insystem[l*721*1440+i*1440+j];}
				if (NWflux<0){
					simulation.budget[l*721*1440+(i+1)*1440+j]+=fabs(NWflux)*simulation.insystem[l*721*1440+i*1440+j];
					simulation.budget[l*721*1440+i*1440+j]-=fabs(NWflux)*simulation.insystem[l*721*1440+i*1440+j];}
				}
				}
			}
		}
		}
		#pragma omp parallel default(shared) private(l,i,j,NW,EW,PW,UW,NWflux,EWflux)
		{
		#pragma omp for schedule(dynamic,20) nowait
		for (i=0;i<721;i++){
			for (j=0;j<1440;j++){
				for (l=0;l<25;l++){
					if (simulation.insystem[l*721*1440+i*1440+j]>0){
				if (simulation.interpolation==0){
					EW=meteo[0].u[(l)*721*1440+i*1440+j]*(*meteo[0].uscale)+(*meteo[0].uadd);}
				if (simulation.interpolation==1){
					EW=(1-timefrac)*(meteo[0].u[(l)*721*1440+i*1440+j]*(*meteo[0].uscale)+(*meteo[0].uadd))+timefrac*(meteo[0].u[(l)*721*1440+i*1440+j]*(*meteo[0].uscale)+(*meteo[0].uadd));
				}
				EWflux=(3.6*EW/simulation.grid_size_EW[i]*dt);
				if (EWflux>0.5){EWflux=0.5;}
				if (EWflux<-0.5){EWflux=-0.5;}
				if (EWflux>0){
					simulation.budget[l*721*1440+i*1440+j+1]+=EWflux*simulation.insystem[l*721*1440+i*1440+j];
					simulation.budget[l*721*1440+i*1440+j]-=EWflux*simulation.insystem[l*721*1440+i*1440+j];
					if (j==1439){
						simulation.budget[l*721*1440+i*1440+j+1]-=EWflux*simulation.insystem[l*721*1440+i*1440+j];
						simulation.budget[l*721*1440+i*1440+j]+=EWflux*simulation.insystem[l*721*1440+i*1440+j];
						simulation.budget[l*721*1440+i*1440]+=EWflux*simulation.insystem[l*721*1440+i*1440+j];
						simulation.budget[l*721*1440+i*1440+j]-=EWflux*simulation.insystem[l*721*1440+i*1440+j];
					}
				}
				if (EWflux<0){
					simulation.budget[l*721*1440+i*1440+j-1]+=fabs(EWflux)*simulation.insystem[l*721*1440+i*1440+j];
					simulation.budget[l*721*1440+i*1440+j]-=fabs(EWflux)*simulation.insystem[l*721*1440+i*1440+j];
					if (j==0){
						simulation.budget[l*721*1440+i*1440+j-1]-=fabs(EWflux)*simulation.insystem[l*721*1440+i*1440+j];
						simulation.budget[l*721*1440+i*1440+j]+=fabs(EWflux)*simulation.insystem[l*721*1440+i*1440+j];
						simulation.budget[l*721*1440+i*1440+1439]+=fabs(EWflux)*simulation.insystem[l*721*1440+i*1440+j];
						simulation.budget[l*721*1440+i*1440+j]-=fabs(EWflux)*simulation.insystem[l*721*1440+i*1440+j];
					}
				}
				}}
			}
		}
		}
		#pragma omp parallel default(shared) private(l,i,j,NW,EW,PW,UW,NWflux,EWflux)
		{
		#pragma omp for schedule(dynamic,20) nowait
		for (i=0;i<721;i++){
			for (j=0;j<1440;j++){
				for (l=0;l<25;l++){
					if (simulation.insystem[l*721*1440+i*1440+j]>0){
				if (simulation.interpolation==0){
					UW=meteo[0].w[(l)*721*1440+i*1440+j]*(*meteo[0].wscale)+(*meteo[0].wadd);}
				if (simulation.interpolation==1){
					UW=(1-timefrac)*meteo[0].w[(l)*721*1440+i*1440+j]*(*meteo[0].wscale)+(*meteo[0].wadd)+timefrac*meteo[1].w[(l)*721*1440+i*1440+j]*(*meteo[1].wscale)+(*meteo[1].wadd);
				}
				if (l>14){
				UWflux=((UW*3600)/2500*dt);
				}
				if (l<15){
				UWflux=((UW*3600)/5000*dt);
				}
				if (UWflux>0.5){UWflux=0.5;}
				if (UWflux<-0.5){UWflux=-0.5;}
				if (UWflux>0&&l<24){
					simulation.budget[(l+1)*721*1440+i*1440+j]+=UWflux*simulation.insystem[l*721*1440+i*1440+j];
					simulation.budget[l*721*1440+i*1440+j]-=UWflux*simulation.insystem[l*721*1440+i*1440+j];
				}
				if (UWflux<0&&l>0){
					simulation.budget[(l-1)*721*1440+i*1440+j]-=UWflux*simulation.insystem[l*721*1440+i*1440+j];
					simulation.budget[l*721*1440+i*1440+j]+=UWflux*simulation.insystem[l*721*1440+i*1440+j];
				}
				}}
			}
		}
		}
		}


		// Update insystem matrix
		if (simulation.simulation_number==6){
		#pragma omp parallel default(shared) private(i,j)
		{
		#pragma omp for schedule(dynamic,20) nowait
		for (i=0;i<721;i++){
			for (j=0;j<1440;j++){
				simulation.insystem_out[0*721*1440+i*1440+j]+=simulation.budget[i*1440+j];
			}}
		}
		}
		
		if (simulation.simulation_number==7){
		#pragma omp parallel default(shared) private(i,j,l)
		{
		#pragma omp for schedule(dynamic,20) nowait
		for (l=0;l<25;l++){
			for (i=0;i<721;i++){
				for (j=0;j<1440;j++){
					simulation.insystem[l*721*1440+i*1440+j]+=simulation.budget[l*721*1440+i*1440+j];
			}}
		}
		}
		}
		
		// Allocate rainout
		//
		if (simulation.simulation_number==6){
		#pragma omp parallel default(shared) private(i,j,P,PW)
		{
		#pragma omp for schedule(dynamic,20) nowait
		for (i=0;i<721;i++){
			for (j=0;j<1440;j++){
				if (simulation.interpolation==0){
					PW=meteo[0].pw[i*1440+j]*(*meteo[0].pwscale)+(*meteo[0].pwadd);
					P=1000*(meteo[0].P[i*1440+j]*(*meteo[0].Pscale)+(*meteo[0].Padd))*dt;}
				if (simulation.interpolation==1){
					PW=(1-timefrac)*(meteo[0].pw[i*1440+j]*(*meteo[0].pwscale)+(*meteo[0].pwadd))+timefrac*(meteo[1].pw[i*1440+j]*(*meteo[1].pwscale)+(*meteo[1].pwadd));
					P=1000*((1-timefrac)*((meteo[0].P[i*1440+j]*(*meteo[0].Pscale)+(*meteo[0].Padd))*dt)+timefrac*((meteo[1].P[i*1440+j]*(*meteo[1].Pscale)+(*meteo[1].Padd))*dt));
				}
				if (P/PW<1&&P/PW>0){
				simulation.allocated_out[0*1440*721+i*1440+j]+=((P/PW)*simulation.insystem_out[0*721*1440+i*1440+j]);
				simulation.insystem_out[0*721*1440+i*1440+j]-=((P/PW)*simulation.insystem_out[0*721*1440+i*1440+j]);}
			}}
		}
		}
		if (simulation.simulation_number==7){
		#pragma omp parallel default(shared) private(i,j,l,P,PW)
		{
		#pragma omp for schedule(dynamic,20) nowait
		for (i=0;i<721;i++){
			for (j=0;j<1440;j++){
				if (simulation.interpolation==0){
					PW=meteo[0].pw[i*1440+j]*(*meteo[0].pwscale)+(*meteo[0].pwadd);
					P=1000*(meteo[0].P[i*1440+j]*(*meteo[0].Pscale)+(*meteo[0].Padd))*dt;}
				if (simulation.interpolation==1){
					PW=(1-timefrac)*(meteo[0].pw[i*1440+j]*(*meteo[0].pwscale)+(*meteo[0].pwadd))+timefrac*(meteo[1].pw[i*1440+j]*(*meteo[1].pwscale)+(*meteo[1].pwadd));
					P=1000*((1-timefrac)*((meteo[0].P[i*1440+j]*(*meteo[0].Pscale)+(*meteo[0].Padd))*dt)+timefrac*((meteo[1].P[i*1440+j]*(*meteo[1].Pscale)+(*meteo[1].Padd))*dt));
				}
				if (P/PW<1&&P/PW>0){
					for (l=0;l<25;l++){
						simulation.allocated_out[0*1440*721+i*1440+j]+=((P/PW)*simulation.insystem[l*721*1440+i*1440+j]);
						simulation.insystem[l*721*1440+i*1440+j]-=((P/PW)*simulation.insystem[l*721*1440+i*1440+j]);
					}
				}
			}
		}
		}
		}
					
		}
					

		return 0;

	}

int eul2layers(struct meteoday* meteo,struct simulation_settings simulation){
	// Simulates the moisture tracking particles for one day using the forcing in meteo and simulation parameters in simulation
	float mean_present,mean_released;
	//int daynum;
	int daynum=simulation.daynum;
	float time;
	float timefrac;
	float dt=0.1;
	int i,j,k,l,iam,np;
	float NW,EW,PW,P,UW;
	float NWflux,EWflux,UWflux;
	float budgetsum=0;
	float statemax=0;
	int imax;

	// Until level 8
	// Clear the budget matrix
	for (i=0;i<2*721*1440;i++){simulation.budget[i]=0;}

	for (time=simulation.hour;time<simulation.hour+1;time+=dt){
		k=(int)floor((time-floor(time))*24);
		timefrac=(time-floor(time))*24-k; // Fraction of the current hour

		// Add moisture to the atmosphere if still releasing
		if (simulation.daynum<simulation.release_days){
			#pragma omp parallel default(shared) private(i,j)
			{
			#pragma omp for schedule(dynamic,20) nowait
			for (i=0;i<721;i++){
				for (j=0;j<1440;j++){
					if (simulation.release_input[i*1440+j]==1){
		simulation.insystem[i*1440+j]+=fabs((meteo[0].E[i*1440+j]*(*meteo[0].Escale)+(*meteo[0].Eadd))*1000.0)*dt;}
		simulation.emitted_out[0*721*1440+i*1440+j]+=fabs((meteo[0].E[i*1440+j]*(*meteo[0].Escale)+(*meteo[0].Eadd))*1000.0)*dt;
					}
				}
			}
			}

		// Move moisture around in atmosphere
		// Clear budget matrix
		#pragma omp parallel default(shared) private(i)
		{
		#pragma omp for schedule(dynamic,20) nowait
		for (i=0;i<2*721*1440;i++){simulation.budget[i]=0;}
		}
		

		float qsum;
		int k;

		// Northwards transport
		#pragma omp parallel default(shared) private(l,i,j,k,qsum,NW,EW,PW,UW,NWflux,EWflux)
		{
		#pragma omp for schedule(dynamic,20) nowait
		for (i=1;i<720;i++){
			for (j=0;j<1440;j++){
				l=0;
				if (simulation.insystem[l*721*1440+i*1440+j]>0){
					NW=0;
					qsum=0;
					for (k=0;k<17;k++){
						NW+=(meteo[0].v[k*721*1440+i*1440+j]*(*meteo[0].vscale)+(*meteo[0].vadd))*(meteo[0].q[k*721*1440+i*1440+j]*(*meteo[0].qscale)+(*meteo[0].qadd));
						qsum+=(meteo[0].q[k*721*1440+i*1440+j]*(*meteo[0].qscale)+(*meteo[0].qadd));
					}
					if (qsum>0){NW/=qsum;}

				NWflux=(3.6*NW/simulation.grid_size_NW[i]*dt);
				if (NWflux>0.5){NWflux=0.5;}
				if (NWflux<-0.5){NWflux=-0.5;}
				if (NWflux>0){
					simulation.budget[l*721*1440+(i-1)*1440+j]+=NWflux*simulation.insystem[l*721*1440+i*1440+j];
					simulation.budget[l*721*1440+i*1440+j]-=NWflux*simulation.insystem[l*721*1440+i*1440+j];}
				if (NWflux<0){
					simulation.budget[l*721*1440+(i+1)*1440+j]+=fabs(NWflux)*simulation.insystem[l*721*1440+i*1440+j];
					simulation.budget[l*721*1440+i*1440+j]-=fabs(NWflux)*simulation.insystem[l*721*1440+i*1440+j];}

		}
				l=1;
				if (simulation.insystem[l*721*1440+i*1440+j]>0){
					NW=0;
					qsum=0;
					for (k=17;k<25;k++){
						NW+=(meteo[0].v[k*721*1440+i*1440+j]*(*meteo[0].vscale)+(*meteo[0].vadd))*(meteo[0].q[k*721*1440+i*1440+j]*(*meteo[0].qscale)+(*meteo[0].qadd));
						qsum+=(meteo[0].q[k*721*1440+i*1440+j]*(*meteo[0].qscale)+(*meteo[0].qadd));
					}
					if (qsum>0){NW/=qsum;}

				NWflux=(3.6*NW/simulation.grid_size_NW[i]*dt);
				if (NWflux>0.5){NWflux=0.5;}
				if (NWflux<-0.5){NWflux=-0.5;}
				if (NWflux>0){
					simulation.budget[l*721*1440+(i-1)*1440+j]+=NWflux*simulation.insystem[l*721*1440+i*1440+j];
					simulation.budget[l*721*1440+i*1440+j]-=NWflux*simulation.insystem[l*721*1440+i*1440+j];}
				if (NWflux<0){
					simulation.budget[l*721*1440+(i+1)*1440+j]+=fabs(NWflux)*simulation.insystem[l*721*1440+i*1440+j];
					simulation.budget[l*721*1440+i*1440+j]-=fabs(NWflux)*simulation.insystem[l*721*1440+i*1440+j];}

		}
		}
		}
		}

		#pragma omp parallel default(shared) private(l,i,j,k,qsum,NW,EW,PW,UW,NWflux,EWflux)
		{
		#pragma omp for schedule(dynamic,20) nowait
		for (i=0;i<721;i++){
			for (j=0;j<1440;j++){
				l=0;
				if (simulation.insystem[l*721*1440+i*1440+j]>0){
					EW=0;
					qsum=0;
					for (k=0;k<17;k++){
					EW+=(meteo[0].u[k*721*1440+i*1440+j]*(*meteo[0].uscale)+(*meteo[0].uadd))*(meteo[0].q[k*721*1440+i*1440+j]*(*meteo[0].qscale)+(*meteo[0].qadd));
					qsum+=(meteo[0].q[k*721*1440+i*1440+j]*(*meteo[0].qscale)+(*meteo[0].qadd));}
					if (qsum>0){EW/=qsum;}
					EWflux=(3.6*EW/simulation.grid_size_EW[i]*dt);
				if (EWflux>0.5){EWflux=0.5;}
				if (EWflux<-0.5){EWflux=-0.5;}
				if (EWflux>0){
					simulation.budget[l*721*1440+i*1440+j+1]+=EWflux*simulation.insystem[l*721*1440+i*1440+j];
					simulation.budget[l*721*1440+i*1440+j]-=EWflux*simulation.insystem[l*721*1440+i*1440+j];
					if (j==1439){
						simulation.budget[l*721*1440+i*1440+j+1]-=EWflux*simulation.insystem[l*721*1440+i*1440+j];
						simulation.budget[l*721*1440+i*1440+j]+=EWflux*simulation.insystem[l*721*1440+i*1440+j];
						simulation.budget[l*721*1440+i*1440]+=EWflux*simulation.insystem[l*721*1440+i*1440+j];
						simulation.budget[l*721*1440+i*1440+j]-=EWflux*simulation.insystem[l*721*1440+i*1440+j];
					}
				}
				if (EWflux<0){
					simulation.budget[l*721*1440+i*1440+j-1]+=fabs(EWflux)*simulation.insystem[l*721*1440+i*1440+j];
					simulation.budget[l*721*1440+i*1440+j]-=fabs(EWflux)*simulation.insystem[l*721*1440+i*1440+j];
					if (j==0){
						simulation.budget[l*721*1440+i*1440+j-1]-=fabs(EWflux)*simulation.insystem[l*721*1440+i*1440+j];
						simulation.budget[l*721*1440+i*1440+j]+=fabs(EWflux)*simulation.insystem[l*721*1440+i*1440+j];
						simulation.budget[l*721*1440+i*1440+1439]+=fabs(EWflux)*simulation.insystem[l*721*1440+i*1440+j];
						simulation.budget[l*721*1440+i*1440+j]-=fabs(EWflux)*simulation.insystem[l*721*1440+i*1440+j];
					}
				}
				}
				l=1;
				if (simulation.insystem[l*721*1440+i*1440+j]>0){
					EW=0;
					qsum=0;
					for (k=17;k<25;k++){
					EW+=(meteo[0].u[k*721*1440+i*1440+j]*(*meteo[0].uscale)+(*meteo[0].uadd))*(meteo[0].q[k*721*1440+i*1440+j]*(*meteo[0].qscale)+(*meteo[0].qadd));
					qsum+=(meteo[0].q[k*721*1440+i*1440+j]*(*meteo[0].qscale)+(*meteo[0].qadd));}
					if (qsum>0){EW/=qsum;}
					EWflux=(3.6*EW/simulation.grid_size_EW[i]*dt);
				if (EWflux>0.5){EWflux=0.5;}
				if (EWflux<-0.5){EWflux=-0.5;}
				if (EWflux>0){
					simulation.budget[l*721*1440+i*1440+j+1]+=EWflux*simulation.insystem[l*721*1440+i*1440+j];
					simulation.budget[l*721*1440+i*1440+j]-=EWflux*simulation.insystem[l*721*1440+i*1440+j];
					if (j==1439){
						simulation.budget[l*721*1440+i*1440+j+1]-=EWflux*simulation.insystem[l*721*1440+i*1440+j];
						simulation.budget[l*721*1440+i*1440+j]+=EWflux*simulation.insystem[l*721*1440+i*1440+j];
						simulation.budget[l*721*1440+i*1440]+=EWflux*simulation.insystem[l*721*1440+i*1440+j];
						simulation.budget[l*721*1440+i*1440+j]-=EWflux*simulation.insystem[l*721*1440+i*1440+j];
					}
				}
				if (EWflux<0){
					simulation.budget[l*721*1440+i*1440+j-1]+=fabs(EWflux)*simulation.insystem[l*721*1440+i*1440+j];
					simulation.budget[l*721*1440+i*1440+j]-=fabs(EWflux)*simulation.insystem[l*721*1440+i*1440+j];
					if (j==0){
						simulation.budget[l*721*1440+i*1440+j-1]-=fabs(EWflux)*simulation.insystem[l*721*1440+i*1440+j];
						simulation.budget[l*721*1440+i*1440+j]+=fabs(EWflux)*simulation.insystem[l*721*1440+i*1440+j];
						simulation.budget[l*721*1440+i*1440+1439]+=fabs(EWflux)*simulation.insystem[l*721*1440+i*1440+j];
						simulation.budget[l*721*1440+i*1440+j]-=fabs(EWflux)*simulation.insystem[l*721*1440+i*1440+j];
					}
				}
				}
			}
		}
		}

		#pragma omp parallel default(shared) private(l,i,j,NW,EW,PW,UW,NWflux,EWflux)
		{
		#pragma omp for schedule(dynamic,20) nowait
		for (i=0;i<721;i++){
			for (j=0;j<1440;j++){
				UW=meteo[0].w[17*721*1440+i*1440+j]*(*meteo[0].wscale)+(*meteo[0].wadd);
				if (UW>0){UWflux=(UW*3600)/80000.0*dt;}
				if (UW<0){UWflux=(UW*3600)/20000.0*dt;}

				if (UWflux>0.5){UWflux=0.5;}
				if (UWflux<-0.5){UWflux=-0.5;}
				if (UWflux>0){
					simulation.budget[0*721*1440+i*1440+j]+=UWflux*simulation.insystem[1*721*1440+i*1440+j];
					simulation.budget[1*721*1440+i*1440+j]-=UWflux*simulation.insystem[1*721*1440+i*1440+j];
				}
				if (UWflux<0){
					simulation.budget[1*721*1440+i*1440+j]+=UWflux*simulation.insystem[0*721*1440+i*1440+j];
					simulation.budget[0*721*1440+i*1440+j]-=UWflux*simulation.insystem[0*721*1440+i*1440+j];
				}
			}
		}
		}
		


		// Update insystem matrix
		#pragma omp parallel default(shared) private(i,j,l)
		{
		#pragma omp for schedule(dynamic,20) nowait
		for (l=0;l<2;l++){
			for (i=0;i<721;i++){
				for (j=0;j<1440;j++){
					simulation.insystem[l*721*1440+i*1440+j]+=simulation.budget[l*721*1440+i*1440+j];
			}}
		}
		}
		
		// Allocate rainout
		//
		#pragma omp parallel default(shared) private(i,j,l,P,PW)
		{
		#pragma omp for schedule(dynamic,20) nowait
		for (i=0;i<721;i++){
			for (j=0;j<1440;j++){
				PW=meteo[0].pw[i*1440+j]*(*meteo[0].pwscale)+(*meteo[0].pwadd);
				P=1000*(meteo[0].P[i*1440+j]*(*meteo[0].Pscale)+(*meteo[0].Padd))*dt;
				if (P/PW<1&&P/PW>0){
					for (l=0;l<2;l++){
						simulation.allocated_out[0*1440*721+i*1440+j]+=((P/PW)*simulation.insystem[l*721*1440+i*1440+j]);
						simulation.insystem[l*721*1440+i*1440+j]-=((P/PW)*simulation.insystem[l*721*1440+i*1440+j]);
					}
				}
			}
		}
		}
					
		}
					

		return 0;

	}

