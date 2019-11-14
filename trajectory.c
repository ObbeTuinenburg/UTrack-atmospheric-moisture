#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include "defines.h"
#include "recycling.h"

int trajectory(struct meteoday* meteo, struct simulation_settings s,int particle_number){
	if (s.particle_present[particle_number]/s.particle_original[particle_number]<0.01){return 0;}
	struct loc curloc;
	float allocatedwater,curlatrad,degreelength;
	int lonlow,latlow;
	struct state curstate;
	struct state curstateint;
	float dt=0.1;
	curloc.lat=s.particle_lat[particle_number];
	curloc.lon=s.particle_lon[particle_number];
	curloc.day=s.daynum;
	curloc.time=s.particle_curtime[particle_number];
	curloc.lev=s.particle_lev[particle_number];
	curloc.theta=s.particle_theta[particle_number];
	while(curloc.time<s.hour+1){
		// Update state
		if (s.interpolation){
		curstate=interpolate(curloc,meteo,s);}
		else {
		curstate=gridcell_state(curloc,meteo,s);}
	
		// Allocate water
		allocatedwater=s.particle_present[particle_number]*(fabs(curstate.P*dt)/curstate.PW);
		latlow=floor((90-curloc.lat)/0.25);
		lonlow=floor((curloc.lon/0.25));
		s.particle_present[particle_number]-=allocatedwater;
		s.allocated_out[0*1440*721+latlow*1440+lonlow]+=allocatedwater;
		s.insystem_out[0*1440*721+latlow*1440+lonlow]+=s.particle_present[particle_number]*dt;
			
		// Update position
		curloc.lon=curloc.lon+curstate.u*dt*3.6/s.grid_size_EW[latlow];

		// Test if particle went over 0 meridian
		while (curloc.lon<0){curloc.lon+=360;}
		while (curloc.lon>360){curloc.lon-=360;}
		
		curloc.lat=curloc.lat+curstate.v*dt*3.6/s.grid_size_NW[latlow];
		if (curloc.lat>89.9){curloc.lat=89.9;}
		if (curloc.lat<-89.9){curloc.lat=-89.9;}
	
		// Update vertical position - using omega
		if (s.simulation_number==2||s.simulation_number>10){
		curloc.lev+=3600*(curstate.w/100.0);
		if (curloc.lev>1000){curloc.lev=1000;}
		if (curloc.lev<50){curloc.lev=50;}
		}
		
		// Update vertical position - using vertical mixing
		if (s.simulation_number==4||s.simulation_number>7){
			if ((float)rand()/RAND_MAX<(1/s.mixing)*dt){
		curloc.lev=get_starting_height(meteo,s,(int)latlow,(int)lonlow);}
		}
		
		curloc.time+=dt;
		}

	// Write back particle location and time
	s.particle_lat[particle_number]=curloc.lat;
	s.particle_lon[particle_number]=curloc.lon;
	s.particle_curtime[particle_number]=curloc.time;
	s.particle_lev[particle_number]=curloc.lev;
	return 0;}
