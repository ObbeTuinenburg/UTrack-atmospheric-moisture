#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <stdlib.h>
#include "netcdf.h"
#include "recycling.h"
#include "defines.h"

int load_meteo_hour(int curyear,int curmonth,int curday,int curhour,struct meteoday* meteo){
	meteo->year=curyear;
	meteo->month=curmonth;
	meteo->day=curday;
	get_era5_hour(curyear,curmonth,curday,curhour,*meteo);
	get_era5_3dt(curyear,curmonth,curday,curhour,*meteo);
	get_era5_3dq(curyear,curmonth,curday,curhour,*meteo);
	get_era5_3dw(curyear,curmonth,curday,curhour,*meteo);
	get_era5_3duv(curyear,curmonth,curday,curhour,*meteo);
	return 0;}

int free_meteo_hour(struct meteoday* meteo){
	free(meteo->ef);
	free(meteo->u);
	free(meteo->v);
	free(meteo->w);
	free(meteo->t);
	free(meteo->q);
	free(meteo->pw);
	free(meteo->nf);
	free(meteo->E);
	free(meteo->P);
	free(meteo->lat);
	free(meteo->lon);
	return 0;}
