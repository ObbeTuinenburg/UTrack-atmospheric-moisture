#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include "netcdf.h"
#include "recycling.h"
#include "defines.h"

int get_input(int *array,char* path){
        int ncid, varid, status;//

    /* The start and count arrays will tell the netCDF library where to
              read our data. */
        size_t start[2], count[2];
           
	/* Error handling. */
           int retval;

           /* Open the file. */

           if ((retval = nc_open(path, NC_NOWRITE, &ncid)))
              ERR(retval);

            count[0] = 721;
            count[1] = 1440;
            start[0] = 0;
            start[1] = 0;


	   if ((retval = nc_inq_varid(ncid, "release" , &varid)))
	      ERR(retval);

               retval = nc_get_vara_int(ncid, varid, start,count,&array[0]);

	return 0;
	}

int get_lsm(short *array,double*add,double*scale){
        int ncid, varid, status;//

    /* The start and count arrays will tell the netCDF library where to
              read our data. */
        size_t start[3], count[3];
	char path[]="./forcing/landmask.nc";
           
	/* Error handling. */
           int retval;

           /* Open the file. */

           if ((retval = nc_open(path, NC_NOWRITE, &ncid)))
              ERR(retval);

            count[0] = 1;
            count[1] = 721;
            count[2] = 1440;
            start[0] = 0;
            start[1] = 0;
            start[2] = 0;


	   if ((retval = nc_inq_varid(ncid, "lsm" , &varid)))
	      ERR(retval);
           retval = nc_get_att_double(ncid, varid, "add_offset", add);
           retval = nc_get_att_double(ncid, varid, "scale_factor", scale);

               retval = nc_get_vara_short(ncid, varid, start,count,&array[0]);

	return 0;
	}

int get_era5_hour(int year,int month,int day,int hour,struct meteoday meteo){
        int ncid, varid, status;
        int lat_varid, lon_varid, lev_varid, time_varid;

        size_t start[3], count[3];
	char path[]=".............................................................................................................................................";
	FILE *ptr;

	if (month>9){
		if (day>9){
		   	sprintf(path, "%s%d%d%d%s", "./forcing/ERA5_2d",day,month,year,".nc");
		}
		else {
		   	sprintf(path, "%s%d%d%d%s", "./forcing/ERA5_2d0",day,month,year,".nc");
		}
	}
		else {
		if (day>9){
		   	sprintf(path, "%s%d%s%d%d%s", "./forcing/ERA5_2d",day,"0",month,year,".nc");}
		else {
		   	sprintf(path, "%s%d%s%d%d%s", "./forcing/ERA5_2d0",day,"0",month,year,".nc");}
		}


           /* Error handling. */
           int retval;

           /* Open the file. */

           if ((retval = nc_open(path, NC_NOWRITE, &ncid)))
              ERR(retval);


           /* Get the varids of the latitude and longitude coordinate
            * variables. */
           if ((retval = nc_inq_varid(ncid, "latitude", &lat_varid)))
              ERR(retval);
           if ((retval = nc_inq_varid(ncid, "longitude", &lon_varid)))
              ERR(retval);

           /* Read the coordinate variable data. */
           if ((retval = nc_get_var_float(ncid, lat_varid, meteo.lat)))
              ERR(retval);
           if ((retval = nc_get_var_float(ncid, lon_varid, meteo.lon)))
              ERR(retval);


            count[0] = 1;
            count[1] = 721;
            count[2] = 1440;
            start[0] = hour;
            start[1] = 0;
            start[2] = 0;


	   if ((retval = nc_inq_varid(ncid, "tcw" , &varid)))
	      ERR(retval);
           retval = nc_get_att_double(ncid, varid, "add_offset", meteo.pwadd);
           retval = nc_get_att_double(ncid, varid, "scale_factor", meteo.pwscale);

               if ((retval = nc_get_vara_float(ncid, varid, start,count,meteo.pw)))
                   ERR(retval);



           /* Get the varids of the latitude and longitude coordinate
            * variables. */
           if ((retval = nc_inq_varid(ncid, "latitude", &lat_varid)))
              ERR(retval);
           if ((retval = nc_inq_varid(ncid, "longitude", &lon_varid)))
              ERR(retval);

           /* Read the coordinate variable data. */
           if ((retval = nc_get_var_float(ncid, lat_varid, meteo.lat)))
              ERR(retval);
           if ((retval = nc_get_var_float(ncid, lon_varid, meteo.lon)))
              ERR(retval);


	   if ((retval = nc_inq_varid(ncid, "p71.162" , &varid)))
	      ERR(retval);
           retval = nc_get_att_double(ncid, varid, "add_offset", meteo.efadd);
           retval = nc_get_att_double(ncid, varid, "scale_factor", meteo.efscale);

               if ((retval = nc_get_vara_float(ncid, varid, start,count,meteo.ef)))
                   ERR(retval);
	   if ((retval = nc_inq_varid(ncid, "p72.162" , &varid)))
	      ERR(retval);
           retval = nc_get_att_double(ncid, varid, "add_offset", meteo.nfadd);
           retval = nc_get_att_double(ncid, varid, "scale_factor", meteo.nfscale);

               if ((retval = nc_get_vara_float(ncid, varid, start,count,meteo.nf)))
                   ERR(retval);
	   
	       if ((retval = nc_inq_varid(ncid, "e" , &varid)))
	      ERR(retval);
           retval = nc_get_att_double(ncid, varid, "add_offset", meteo.Eadd);
           retval = nc_get_att_double(ncid, varid, "scale_factor", meteo.Escale);

               if ((retval = nc_get_vara_float(ncid, varid, start,count,meteo.E)))
                   ERR(retval);
	       
	       if ((retval = nc_inq_varid(ncid, "tp" , &varid)))
	      ERR(retval);
           retval = nc_get_att_double(ncid, varid, "add_offset", meteo.Padd);
           retval = nc_get_att_double(ncid, varid, "scale_factor", meteo.Pscale);

               if ((retval = nc_get_vara_float(ncid, varid, start,count,meteo.P)))
                   ERR(retval);
	       retval = nc_close(ncid);
    
	return 0;}

int get_era5_3dq(int year,int month,int day,int hour,struct meteoday meteo){
        int ncid, varid, status;//
        int lat_varid, lon_varid, lev_varid, time_varid;

    /* The start and count arrays will tell the netCDF library where to
              read our data. */
        size_t start[4], count[4];
	char path[]=".............................................................................................................................................";
	FILE *ptr;
            count[0] = 1;
            count[1] = 25;
            count[2] = 721;
            count[3] = 1440;
            start[0] = hour;
            start[1] = 0;
            start[2] = 0;
            start[3] = 0;
           /* Error handling. */
           int retval;


	   /* Read in q data */
	if (month>9){
		if (day>9){
		   	sprintf(path, "%s%d%d%d%s", "./forcing/ERA5_q",day,month,year,".nc");
		}
		else {
		   	sprintf(path, "%s%d%d%d%s", "./forcing/ERA5_q0",day,month,year,".nc");
		}
	}
		else {
		if (day>9){
		   	sprintf(path, "%s%d%s%d%d%s", "./forcing/ERA5_q",day,"0",month,year,".nc");}
		else {
		   	sprintf(path, "%s%d%s%d%d%s", "./forcing/ERA5_q0",day,"0",month,year,".nc");}
		}


           /* Open the file. */

           if ((retval = nc_open(path, NC_NOWRITE, &ncid)))
              ERR(retval);

	   if ((retval = nc_inq_varid(ncid, "q" , &varid)))
	      ERR(retval);
           retval = nc_get_att_double(ncid, varid, "add_offset", meteo.qadd);
           retval = nc_get_att_double(ncid, varid, "scale_factor", meteo.qscale);

               if ((retval = nc_get_vara_float(ncid, varid, start,count,meteo.q)))
                   ERR(retval);
	       /* Close the file */
	       retval = nc_close(ncid);
	return 0;}

int get_era5_3dt(int year,int month,int day,int hour,struct meteoday meteo){
        int ncid, varid, status;//
        int lat_varid, lon_varid, lev_varid, time_varid;

    /* The start and count arrays will tell the netCDF library where to
              read our data. */
        size_t start[4], count[4];
	char path[]=".............................................................................................................................................";
	FILE *ptr;
            count[0] = 1;
            count[1] = 25;
            count[2] = 721;
            count[3] = 1440;
            start[0] = hour;
            start[1] = 0;
            start[2] = 0;
            start[3] = 0;
           /* Error handling. */
           int retval;


	   /* Read in w data */
	if (month>9){
		if (day>9){
		   	sprintf(path, "%s%d%d%d%s", "./forcing/ERA5_t",day,month,year,".nc");
		}
		else {
		   	sprintf(path, "%s%d%d%d%s", "./forcing/ERA5_t0",day,month,year,".nc");
		}
	}
		else {
		if (day>9){
		   	sprintf(path, "%s%d%s%d%d%s", "./forcing/ERA5_t",day,"0",month,year,".nc");}
		else {
		   	sprintf(path, "%s%d%s%d%d%s", "./forcing/ERA5_t0",day,"0",month,year,".nc");}
		}


           /* Open the file. */

           if ((retval = nc_open(path, NC_NOWRITE, &ncid)))
              ERR(retval);

	   if ((retval = nc_inq_varid(ncid, "t" , &varid)))
	      ERR(retval);
           retval = nc_get_att_double(ncid, varid, "add_offset", meteo.tadd);
           retval = nc_get_att_double(ncid, varid, "scale_factor", meteo.tscale);

               if ((retval = nc_get_vara_float(ncid, varid, start,count,meteo.t)))
                   ERR(retval);
	       /* Close the file */
	       retval = nc_close(ncid);
	return 0;}

int get_era5_3dw(int year,int month,int day,int hour,struct meteoday meteo){
        int ncid, varid, status;//
        int lat_varid, lon_varid, lev_varid, time_varid;

    /* The start and count arrays will tell the netCDF library where to
              read our data. */
        size_t start[4], count[4];
	char path[]=".............................................................................................................................................";
	FILE *ptr;
            count[0] = 1;
            count[1] = 25;
            count[2] = 721;
            count[3] = 1440;
            start[0] = hour;
            start[1] = 0;
            start[2] = 0;
            start[3] = 0;
           /* Error handling. */
           int retval;


	   /* Read in w data */
	if (month>9){
		if (day>9){
		   	sprintf(path, "%s%d%d%d%s", "./forcing/ERA5_w",day,month,year,".nc");
		}
		else {
		   	sprintf(path, "%s%d%d%d%s", "./forcing/ERA5_w0",day,month,year,".nc");
		}
	}
		else {
		if (day>9){
		   	sprintf(path, "%s%d%s%d%d%s", "./forcing/ERA5_w",day,"0",month,year,".nc");}
		else {
		   	sprintf(path, "%s%d%s%d%d%s", "./forcing/ERA5_w0",day,"0",month,year,".nc");}
		}


           /* Open the file. */

           if ((retval = nc_open(path, NC_NOWRITE, &ncid)))
              ERR(retval);

	   if ((retval = nc_inq_varid(ncid, "w" , &varid)))
	      ERR(retval);
           retval = nc_get_att_double(ncid, varid, "add_offset", meteo.wadd);
           retval = nc_get_att_double(ncid, varid, "scale_factor", meteo.wscale);

               if ((retval = nc_get_vara_float(ncid, varid, start,count,meteo.w)))
                   ERR(retval);
	       /* Close the file */
	       retval = nc_close(ncid);
	return 0;}

int get_era5_3duv(int year,int month,int day,int hour,struct meteoday meteo){
        int ncid, varid, status;//
        int lat_varid, lon_varid, lev_varid, time_varid;

    /* The start and count arrays will tell the netCDF library where to
              read our data. */
        size_t start[4], count[4];
	char path[]=".............................................................................................................................................";
	FILE *ptr;
            count[0] = 1;
            count[1] = 25;
            count[2] = 721;
            count[3] = 1440;
            start[0] = hour;
            start[1] = 0;
            start[2] = 0;
            start[3] = 0;
           /* Error handling. */
           int retval;

	   /* Read in uv data */
	if (month>9){
		if (day>9){
		   	sprintf(path, "%s%d%d%d%s", "./forcing/ERA5_uv",day,month,year,".nc");
		}
		else {
		   	sprintf(path, "%s%d%d%d%s", "./forcing/ERA5_uv0",day,month,year,".nc");
		}
	}
		else {
		if (day>9){
		   	sprintf(path, "%s%d%s%d%d%s", "./forcing/ERA5_uv",day,"0",month,year,".nc");}
		else {
		   	sprintf(path, "%s%d%s%d%d%s", "./forcing/ERA5_uv0",day,"0",month,year,".nc");}
		}


           /* Open the file. */

           if ((retval = nc_open(path, NC_NOWRITE, &ncid)))
              ERR(retval);
	   if ((retval = nc_inq_varid(ncid, "v" , &varid)))
	      ERR(retval);
           retval = nc_get_att_double(ncid, varid, "add_offset", meteo.vadd);
           retval = nc_get_att_double(ncid, varid, "scale_factor", meteo.vscale);

               if ((retval = nc_get_vara_float(ncid, varid, start,count,meteo.v)))
                   ERR(retval);

	   if ((retval = nc_inq_varid(ncid, "u" , &varid)))
	      ERR(retval);
           retval = nc_get_att_double(ncid, varid, "add_offset", meteo.uadd);
           retval = nc_get_att_double(ncid, varid, "scale_factor", meteo.uscale);

               if ((retval = nc_get_vara_float(ncid, varid, start,count,meteo.u)))
                   ERR(retval);
	       /* Close the file */
	//       retval = nc_close(ncid);

	return 0;}
