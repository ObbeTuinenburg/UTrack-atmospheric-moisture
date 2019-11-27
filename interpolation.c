#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include "netcdf.h"
#include "recycling.h"
#include "defines.h"

struct state gridcell_state(struct loc curloc,struct meteoday* meteo, struct simulation_settings s){
	struct state outstate;
	int curlatidx,curlonidx,curlevidx,i;
	curlatidx=floor((90-curloc.lat)/0.25);
	curlonidx=floor((curloc.lon)/0.25);
	i=curlatidx*1440+curlonidx;	
	outstate.E=(meteo[0].E[i]*(*meteo[0].Escale)+(*meteo[0].Eadd))*1000.0;
	outstate.P=(meteo[0].P[i]*(*meteo[0].Pscale)+(*meteo[0].Padd))*1000.0;
	outstate.PW=meteo[0].pw[i]*(*meteo[0].pwscale)+(*meteo[0].pwadd);

	if (s.simulation_number==1){
	outstate.u=(meteo[0].ef[i]*(*meteo[0].efscale)+(*meteo[0].efadd))/outstate.PW;
	outstate.v=(meteo[0].nf[i]*(*meteo[0].nfscale)+(*meteo[0].nfadd))/outstate.PW;}

	if (s.simulation_number>1){ // 3D simulation
		if (s.simulation_number==2||s.simulation_number==4||s.simulation_number==5||s.simulation_number>7){ // 
			if (curloc.lev<750){
			curlevidx=(int)floor((curloc.lev-50)/50.0);}
			if (curloc.lev>=750){
			curlevidx=(int)(15+floor((curloc.lev-750)/25.0));}

		}
		i+=curlevidx*1440*721;
		outstate.u=(meteo[0].u[i]*(*meteo[0].uscale)+(*meteo[0].uadd));
		outstate.v=(meteo[0].v[i]*(*meteo[0].vscale)+(*meteo[0].vadd));
		outstate.w=(meteo[0].w[i]*(*meteo[0].wscale)+(*meteo[0].wadd));
	}
	return outstate;
	}

struct state interpolate(struct loc curloc,struct meteoday* meteo, struct simulation_settings s){
	struct state outstate;
	struct state *poutstate=&outstate;
	int curlatidx,curlonidx,curlevidx,i;
	float latfrac,lonfrac,levfrac,timefrac;
	curlatidx=floor((90-curloc.lat)/0.25);
	latfrac=((90-curloc.lat)/0.25)-curlatidx;
	curlonidx=floor((curloc.lon)/0.25);
	lonfrac=((curloc.lon)/0.25)-curlonidx;
	timefrac=curloc.time-floor(curloc.time);
	i=curlatidx*1440+curlonidx;	
	
	i=get_2d_interpolation(poutstate,meteo,curlatidx,curlonidx,latfrac,lonfrac,timefrac);
	

	if (s.simulation_number==1){
	outstate.u=outstate.ef/outstate.PW;
	outstate.v=outstate.nf/outstate.PW;}

	if (s.simulation_number>1){ // 3D simulation
			if (s.simulation_number==2||s.simulation_number==4||s.simulation_number==5||s.simulation_number>7){ // 
			curlevidx=(int)floor((curloc.lev-500)/50.0);
			if (curloc.lev<750){
			curlevidx=(int)floor((curloc.lev-50)/50.0);
			levfrac=((curloc.lev-50)/50.0)-curlevidx;}
			if (curloc.lev>=750){
			curlevidx=(int)(15+floor((curloc.lev-750)/25.0));
			levfrac=(15+floor((curloc.lev-750)/25.0))-curlevidx;}
		}
		i+=curlevidx*1440*721;
		i=get_3d_interpolation(poutstate,meteo,curlatidx,curlonidx,curlevidx,latfrac,lonfrac,levfrac,timefrac);
	}
	return outstate;
	}

int get_2d_interpolation(struct state* outstate,struct meteoday* meteo,int curlatidx,int curlonidx,float latfrac,float lonfrac,float timefrac){
	float latlowval,lathighval,lonlowval,lonhighval,timelowval,timehighval;
	int highlat,highlon;
	highlat=curlatidx+1;
	highlon=curlonidx+1;
	if (highlat>720){highlat=720;}
	if (highlon>1439){highlat=0;}
	// E interpolation
	latlowval=meteo[0].E[curlatidx*1440+curlonidx]*(*meteo[0].Escale)+(*meteo[0].Eadd);
	lathighval=meteo[0].E[(highlat)*1440+curlonidx]*(*meteo[0].Escale)+(*meteo[0].Eadd);
	lonlowval=(1-latfrac)*latlowval+latfrac*lathighval;
	latlowval=meteo[0].E[curlatidx*1440+highlon]*(*meteo[0].Escale)+(*meteo[0].Eadd);
	lathighval=meteo[0].E[(highlat)*1440+highlon]*(*meteo[0].Escale)+(*meteo[0].Eadd);
	lonhighval=(1-latfrac)*latlowval+latfrac*lathighval;
	timelowval=(1-lonfrac)*lonlowval+lonfrac*lonhighval;
	latlowval=meteo[1].E[curlatidx*1440+curlonidx]*(*meteo[1].Escale)+(*meteo[1].Eadd);
	lathighval=meteo[1].E[(highlat)*1440+curlonidx]*(*meteo[1].Escale)+(*meteo[1].Eadd);
	lonlowval=(1-latfrac)*latlowval+latfrac*lathighval;
	latlowval=meteo[1].E[curlatidx*1440+highlon]*(*meteo[1].Escale)+(*meteo[1].Eadd);
	lathighval=meteo[1].E[(highlat)*1440+highlon]*(*meteo[1].Escale)+(*meteo[1].Eadd);
	lonhighval=(1-latfrac)*latlowval+latfrac*lathighval;
	timehighval=(1-lonfrac)*lonlowval+lonfrac*lonhighval;
	(*outstate).E=1000*((1-timefrac)*timelowval+timefrac*timehighval);
	
	// PW interpolation
	latlowval=meteo[0].pw[curlatidx*1440+curlonidx]*(*meteo[0].pwscale)+(*meteo[0].pwadd);
	lathighval=meteo[0].pw[(highlat)*1440+curlonidx]*(*meteo[0].pwscale)+(*meteo[0].pwadd);
	//printf("PW: %f %f\n",latlowval,lathighval);
	lonlowval=(1-latfrac)*latlowval+latfrac*lathighval;
	latlowval=meteo[0].pw[curlatidx*1440+highlon]*(*meteo[0].pwscale)+(*meteo[0].pwadd);
	lathighval=meteo[0].pw[(highlat)*1440+highlon]*(*meteo[0].pwscale)+(*meteo[0].pwadd);
	//printf("PW: %f %f\n",latlowval,lathighval);
	lonhighval=(1-latfrac)*latlowval+latfrac*lathighval;
	timelowval=(1-lonfrac)*lonlowval+lonfrac*lonhighval;
	latlowval=meteo[1].pw[curlatidx*1440+curlonidx]*(*meteo[1].pwscale)+(*meteo[1].pwadd);
	lathighval=meteo[1].pw[(highlat)*1440+curlonidx]*(*meteo[1].pwscale)+(*meteo[1].pwadd);
	//printf("PW: %f %f\n",latlowval,lathighval);
	lonlowval=(1-latfrac)*latlowval+latfrac*lathighval;
	latlowval=meteo[1].pw[curlatidx*1440+highlon]*(*meteo[1].pwscale)+(*meteo[1].pwadd);
	lathighval=meteo[1].pw[(highlat)*1440+highlon]*(*meteo[1].pwscale)+(*meteo[1].pwadd);
	//printf("PW: %f %f\n",latlowval,lathighval);
	lonhighval=(1-latfrac)*latlowval+latfrac*lathighval;
	timehighval=(1-lonfrac)*lonlowval+lonfrac*lonhighval;
	(*outstate).PW=((1-timefrac)*timelowval+timefrac*timehighval);
	
	// nf interpolation
	latlowval=meteo[0].nf[curlatidx*1440+curlonidx]*(*meteo[0].nfscale)+(*meteo[0].nfadd);
	lathighval=meteo[0].nf[(highlat)*1440+curlonidx]*(*meteo[0].nfscale)+(*meteo[0].nfadd);
	lonlowval=(1-latfrac)*latlowval+latfrac*lathighval;
	latlowval=meteo[0].nf[curlatidx*1440+highlon]*(*meteo[0].nfscale)+(*meteo[0].nfadd);
	lathighval=meteo[0].nf[(highlat)*1440+highlon]*(*meteo[0].nfscale)+(*meteo[0].nfadd);
	lonhighval=(1-latfrac)*latlowval+latfrac*lathighval;
	timelowval=(1-lonfrac)*lonlowval+lonfrac*lonhighval;
	latlowval=meteo[1].nf[curlatidx*1440+curlonidx]*(*meteo[1].nfscale)+(*meteo[1].nfadd);
	lathighval=meteo[1].nf[(highlat)*1440+curlonidx]*(*meteo[1].nfscale)+(*meteo[1].nfadd);
	lonlowval=(1-latfrac)*latlowval+latfrac*lathighval;
	latlowval=meteo[1].nf[curlatidx*1440+highlon]*(*meteo[1].nfscale)+(*meteo[1].nfadd);
	lathighval=meteo[1].nf[(highlat)*1440+highlon]*(*meteo[1].nfscale)+(*meteo[1].nfadd);
	lonhighval=(1-latfrac)*latlowval+latfrac*lathighval;
	timehighval=(1-lonfrac)*lonlowval+lonfrac*lonhighval;
	(*outstate).nf=((1-timefrac)*timelowval+timefrac*timehighval);
	
	
	 // ef interpolation
	latlowval=meteo[0].ef[curlatidx*1440+curlonidx]*(*meteo[0].efscale)+(*meteo[0].efadd);
	lathighval=meteo[0].ef[(highlat)*1440+curlonidx]*(*meteo[0].efscale)+(*meteo[0].efadd);
	lonlowval=(1-latfrac)*latlowval+latfrac*lathighval;
	latlowval=meteo[0].ef[curlatidx*1440+highlon]*(*meteo[0].efscale)+(*meteo[0].efadd);
	lathighval=meteo[0].ef[(highlat)*1440+highlon]*(*meteo[0].efscale)+(*meteo[0].efadd);
	lonhighval=(1-latfrac)*latlowval+latfrac*lathighval;
	timelowval=(1-lonfrac)*lonlowval+lonfrac*lonhighval;
	latlowval=meteo[1].ef[curlatidx*1440+curlonidx]*(*meteo[1].efscale)+(*meteo[1].efadd);
	lathighval=meteo[1].ef[(highlat)*1440+curlonidx]*(*meteo[1].efscale)+(*meteo[1].efadd);
	lonlowval=(1-latfrac)*latlowval+latfrac*lathighval;
	latlowval=meteo[1].ef[curlatidx*1440+highlon]*(*meteo[1].efscale)+(*meteo[1].efadd);
	lathighval=meteo[1].ef[(highlat)*1440+highlon]*(*meteo[1].efscale)+(*meteo[1].efadd);
	lonhighval=(1-latfrac)*latlowval+latfrac*lathighval;
	timehighval=(1-lonfrac)*lonlowval+lonfrac*lonhighval;
	(*outstate).ef=((1-timefrac)*timelowval+timefrac*timehighval);

	// P interpolation
	latlowval=meteo[0].P[curlatidx*1440+curlonidx]*(*meteo[0].Pscale)+(*meteo[0].Padd);
	lathighval=meteo[0].P[(highlat)*1440+curlonidx]*(*meteo[0].Pscale)+(*meteo[0].Padd);
	lonlowval=(1-latfrac)*latlowval+latfrac*lathighval;
	latlowval=meteo[0].P[curlatidx*1440+highlon]*(*meteo[0].Pscale)+(*meteo[0].Padd);
	lathighval=meteo[0].P[(highlat)*1440+highlon]*(*meteo[0].Pscale)+(*meteo[0].Padd);
	lonhighval=(1-latfrac)*latlowval+latfrac*lathighval;
	timelowval=(1-lonfrac)*lonlowval+lonfrac*lonhighval;
	latlowval=meteo[1].P[curlatidx*1440+curlonidx]*(*meteo[1].Pscale)+(*meteo[1].Padd);
	lathighval=meteo[1].P[(highlat)*1440+curlonidx]*(*meteo[1].Pscale)+(*meteo[1].Padd);
	lonlowval=(1-latfrac)*latlowval+latfrac*lathighval;
	latlowval=meteo[1].P[curlatidx*1440+highlon]*(*meteo[1].Pscale)+(*meteo[1].Padd);
	lathighval=meteo[1].P[(highlat)*1440+highlon]*(*meteo[1].Pscale)+(*meteo[1].Padd);
	lonhighval=(1-latfrac)*latlowval+latfrac*lathighval;
	timehighval=(1-lonfrac)*lonlowval+lonfrac*lonhighval;
	(*outstate).P=1000*((1-timefrac)*timelowval+timefrac*timehighval);


	return 0;}

int get_3d_interpolation(struct state* outstate,struct meteoday* meteo,int curlatidx,int curlonidx,int curlevidx,float latfrac,float lonfrac,float levfrac,float timefrac){
	float latlowval,lathighval,lonlowval,lonhighval,levlowval,levhighval,timelowval,timehighval;
	int highlev,highlat,highlon;
	highlev=curlevidx+1;
	highlat=curlatidx+1;
	highlon=curlonidx+1;
	//if (highlev>9){highlev=9;}
	if (highlev>24){highlev=24;}
	if (highlat>720){highlat=720;}
	if (highlon>1439){highlat=0;}

	// U interpolation
	latlowval=meteo[0].u[curlevidx*1440*721+curlatidx*1440+curlonidx]*(*meteo[0].uscale)+(*meteo[0].uadd);
	lathighval=meteo[0].u[curlevidx*1440*721+(highlat)*1440+curlonidx]*(*meteo[0].uscale)+(*meteo[0].uadd);
	lonlowval=(1-latfrac)*latlowval+latfrac*lathighval;
	//printf("u: %f %f\n",latlowval,lathighval);
	latlowval=meteo[0].u[curlevidx*1440*721+curlatidx*1440+highlon]*(*meteo[0].uscale)+(*meteo[0].uadd);
	lathighval=meteo[0].u[curlevidx*1440*721+(highlat)*1440+highlon]*(*meteo[0].uscale)+(*meteo[0].uadd);
	//printf("u: %f %f\n",latlowval,lathighval);
	lonhighval=(1-latfrac)*latlowval+latfrac*lathighval;
	levlowval=(1-lonfrac)*lonlowval+lonfrac*lonhighval;
	latlowval=meteo[0].u[(highlev)*1440*721+curlatidx*1440+curlonidx]*(*meteo[0].uscale)+(*meteo[0].uadd);
	lathighval=meteo[0].u[(highlev)*1440*721+(highlat)*1440+curlonidx]*(*meteo[0].uscale)+(*meteo[0].uadd);
	//printf("u: %f %f\n",latlowval,lathighval);
	lonlowval=(1-latfrac)*latlowval+latfrac*lathighval;
	latlowval=meteo[0].u[(highlev)*1440*721+curlatidx*1440+highlon]*(*meteo[0].uscale)+(*meteo[0].uadd);
	lathighval=meteo[0].u[(highlev)*1440*721+(highlat)*1440+highlon]*(*meteo[0].uscale)+(*meteo[0].uadd);
	//printf("u: %f %f\n",latlowval,lathighval);
	lonhighval=(1-latfrac)*latlowval+latfrac*lathighval;
	levhighval=(1-lonfrac)*lonlowval+lonfrac*lonhighval;
	timelowval=(1-levfrac)*levlowval+levfrac*levhighval;
	latlowval=meteo[1].u[curlevidx*1440*721+curlatidx*1440+curlonidx]*(*meteo[1].uscale)+(*meteo[1].uadd);
	lathighval=meteo[1].u[curlevidx*1440*721+(highlat)*1440+curlonidx]*(*meteo[1].uscale)+(*meteo[1].uadd);
	//printf("u: %f %f\n",latlowval,lathighval);
	lonlowval=(1-latfrac)*latlowval+latfrac*lathighval;
	latlowval=meteo[1].u[curlevidx*1440*721+curlatidx*1440+highlon]*(*meteo[1].uscale)+(*meteo[1].uadd);
	lathighval=meteo[1].u[curlevidx*1440*721+(highlat)*1440+highlon]*(*meteo[1].uscale)+(*meteo[1].uadd);
	//printf("u: %f %f\n",latlowval,lathighval);
	lonhighval=(1-latfrac)*latlowval+latfrac*lathighval;
	levlowval=(1-lonfrac)*lonlowval+lonfrac*lonhighval;
	latlowval=meteo[1].u[(highlev)*1440*721+curlatidx*1440+curlonidx]*(*meteo[1].uscale)+(*meteo[1].uadd);
	lathighval=meteo[1].u[(highlev)*1440*721+(highlat)*1440+curlonidx]*(*meteo[1].uscale)+(*meteo[1].uadd);
	//printf("u: %f %f\n",latlowval,lathighval);
	lonlowval=(1-latfrac)*latlowval+latfrac*lathighval;
	latlowval=meteo[1].u[(highlev)*1440*721+curlatidx*1440+highlon]*(*meteo[1].uscale)+(*meteo[1].uadd);
	lathighval=meteo[1].u[(highlev)*1440*721+(highlat)*1440+highlon]*(*meteo[1].uscale)+(*meteo[1].uadd);
	//printf("u: %f %f %d\n",latlowval,lathighval,curlevidx);
	lonhighval=(1-latfrac)*latlowval+latfrac*lathighval;
	levhighval=(1-lonfrac)*lonlowval+lonfrac*lonhighval;
	timehighval=(1-levfrac)*levlowval+levfrac*levhighval;
	(*outstate).u=((1-timefrac)*timelowval+timefrac*timehighval);
	
	// V interpolation
	latlowval=meteo[0].v[curlevidx*1440*721+curlatidx*1440+curlonidx]*(*meteo[0].vscale)+(*meteo[0].vadd);
	lathighval=meteo[0].v[curlevidx*1440*721+(highlat)*1440+curlonidx]*(*meteo[0].vscale)+(*meteo[0].vadd);
	lonlowval=(1-latfrac)*latlowval+latfrac*lathighval;
	latlowval=meteo[0].v[curlevidx*1440*721+curlatidx*1440+highlon]*(*meteo[0].vscale)+(*meteo[0].vadd);
	lathighval=meteo[0].v[curlevidx*1440*721+(highlat)*1440+highlon]*(*meteo[0].vscale)+(*meteo[0].vadd);
	lonhighval=(1-latfrac)*latlowval+latfrac*lathighval;
	levlowval=(1-lonfrac)*lonlowval+lonfrac*lonhighval;
	latlowval=meteo[0].v[(highlev)*1440*721+curlatidx*1440+curlonidx]*(*meteo[0].vscale)+(*meteo[0].vadd);
	lathighval=meteo[0].v[(highlev)*1440*721+(highlat)*1440+curlonidx]*(*meteo[0].vscale)+(*meteo[0].vadd);
	lonlowval=(1-latfrac)*latlowval+latfrac*lathighval;
	latlowval=meteo[0].v[(highlev)*1440*721+curlatidx*1440+highlon]*(*meteo[0].vscale)+(*meteo[0].vadd);
	lathighval=meteo[0].v[(highlev)*1440*721+(highlat)*1440+highlon]*(*meteo[0].vscale)+(*meteo[0].vadd);
	lonhighval=(1-latfrac)*latlowval+latfrac*lathighval;
	levhighval=(1-lonfrac)*lonlowval+lonfrac*lonhighval;
	timelowval=(1-levfrac)*levlowval+levfrac*levhighval;
	latlowval=meteo[1].v[curlevidx*1440*721+curlatidx*1440+curlonidx]*(*meteo[1].vscale)+(*meteo[1].vadd);
	lathighval=meteo[1].v[curlevidx*1440*721+(highlat)*1440+curlonidx]*(*meteo[1].vscale)+(*meteo[1].vadd);
	lonlowval=(1-latfrac)*latlowval+latfrac*lathighval;
	latlowval=meteo[1].v[curlevidx*1440*721+curlatidx*1440+highlon]*(*meteo[1].vscale)+(*meteo[1].vadd);
	lathighval=meteo[1].v[curlevidx*1440*721+(highlat)*1440+highlon]*(*meteo[1].vscale)+(*meteo[1].vadd);
	lonhighval=(1-latfrac)*latlowval+latfrac*lathighval;
	levlowval=(1-lonfrac)*lonlowval+lonfrac*lonhighval;
	latlowval=meteo[1].v[(highlev)*1440*721+curlatidx*1440+curlonidx]*(*meteo[1].vscale)+(*meteo[1].vadd);
	lathighval=meteo[1].v[(highlev)*1440*721+(highlat)*1440+curlonidx]*(*meteo[1].vscale)+(*meteo[1].vadd);
	lonlowval=(1-latfrac)*latlowval+latfrac*lathighval;
	latlowval=meteo[1].v[(highlev)*1440*721+curlatidx*1440+highlon]*(*meteo[1].vscale)+(*meteo[1].vadd);
	lathighval=meteo[1].v[(highlev)*1440*721+(highlat)*1440+highlon]*(*meteo[1].vscale)+(*meteo[1].vadd);
	lonhighval=(1-latfrac)*latlowval+latfrac*lathighval;
	levhighval=(1-lonfrac)*lonlowval+lonfrac*lonhighval;
	timehighval=(1-levfrac)*levlowval+levfrac*levhighval;
	(*outstate).v=((1-timefrac)*timelowval+timefrac*timehighval);

	// W interpolation
	latlowval=meteo[0].w[curlevidx*1440*721+curlatidx*1440+curlonidx]*(*meteo[0].wscale)+(*meteo[0].wadd);
	lathighval=meteo[0].w[curlevidx*1440*721+(highlat)*1440+curlonidx]*(*meteo[0].wscale)+(*meteo[0].wadd);
	lonlowval=(1-latfrac)*latlowval+latfrac*lathighval;
	latlowval=meteo[0].w[curlevidx*1440*721+curlatidx*1440+highlon]*(*meteo[0].wscale)+(*meteo[0].wadd);
	lathighval=meteo[0].w[curlevidx*1440*721+(highlat)*1440+highlon]*(*meteo[0].wscale)+(*meteo[0].wadd);
	lonhighval=(1-latfrac)*latlowval+latfrac*lathighval;
	levlowval=(1-lonfrac)*lonlowval+lonfrac*lonhighval;
	latlowval=meteo[0].w[(highlev)*1440*721+curlatidx*1440+curlonidx]*(*meteo[0].wscale)+(*meteo[0].wadd);
	lathighval=meteo[0].w[(highlev)*1440*721+(highlat)*1440+curlonidx]*(*meteo[0].wscale)+(*meteo[0].wadd);
	lonlowval=(1-latfrac)*latlowval+latfrac*lathighval;
	latlowval=meteo[0].w[(highlev)*1440*721+curlatidx*1440+highlon]*(*meteo[0].wscale)+(*meteo[0].wadd);
	lathighval=meteo[0].w[(highlev)*1440*721+(highlat)*1440+highlon]*(*meteo[0].wscale)+(*meteo[0].wadd);
	lonhighval=(1-latfrac)*latlowval+latfrac*lathighval;
	levhighval=(1-lonfrac)*lonlowval+lonfrac*lonhighval;
	timelowval=(1-levfrac)*levlowval+levfrac*levhighval;
	latlowval=meteo[1].w[curlevidx*1440*721+curlatidx*1440+curlonidx]*(*meteo[1].wscale)+(*meteo[1].wadd);
	lathighval=meteo[1].w[curlevidx*1440*721+(highlat)*1440+curlonidx]*(*meteo[1].wscale)+(*meteo[1].wadd);
	lonlowval=(1-latfrac)*latlowval+latfrac*lathighval;
	latlowval=meteo[1].w[curlevidx*1440*721+curlatidx*1440+highlon]*(*meteo[1].wscale)+(*meteo[1].wadd);
	lathighval=meteo[1].w[curlevidx*1440*721+(highlat)*1440+highlon]*(*meteo[1].wscale)+(*meteo[1].wadd);
	lonhighval=(1-latfrac)*latlowval+latfrac*lathighval;
	levlowval=(1-lonfrac)*lonlowval+lonfrac*lonhighval;
	latlowval=meteo[1].w[(highlev)*1440*721+curlatidx*1440+curlonidx]*(*meteo[1].wscale)+(*meteo[1].wadd);
	lathighval=meteo[1].w[(highlev)*1440*721+(highlat)*1440+curlonidx]*(*meteo[1].wscale)+(*meteo[1].wadd);
	lonlowval=(1-latfrac)*latlowval+latfrac*lathighval;
	latlowval=meteo[1].w[(highlev)*1440*721+curlatidx*1440+highlon]*(*meteo[1].wscale)+(*meteo[1].wadd);
	lathighval=meteo[1].w[(highlev)*1440*721+(highlat)*1440+highlon]*(*meteo[1].wscale)+(*meteo[1].wadd);
	lonhighval=(1-latfrac)*latlowval+latfrac*lathighval;
	levhighval=(1-lonfrac)*lonlowval+lonfrac*lonhighval;
	timehighval=(1-levfrac)*levlowval+levfrac*levhighval;
	(*outstate).w=((1-timefrac)*timelowval+timefrac*timehighval);
	return 0;}

