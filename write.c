#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include "netcdf.h"
#include "recycling.h"
#include "defines.h"

int write_output(double *released,double *insystem,double *allocated,float *lats,float *lons,int ndays,char* outfile, float runtime){
   
/* IDs for the netCDF file, dimensions, and variables. */
   int ncid, lon_dimid, lat_dimid, rec_dimid;
   int lat_varid, lon_varid, pres_varid, temp_varid,rel_varid;
   int dimids[2];

   /* The start and count arrays will tell the netCDF library where to
      write our data. */
   size_t start[2], count[2];


   /* Loop indexes. */
   int lat, lon, rec, i = 0;
   
   /* Error handling. */
   int retval;
   
	/* Create the file. */
   if ((retval = nc_create(outfile, NC_CLOBBER, &ncid)))
      ERR(retval);

   /* Define the dimensions. The record dimension is defined to have
    * unlimited length - it can grow as needed. In this example it is
    * the time dimension.*/
   if ((retval = nc_def_dim(ncid, "latitude", 721, &lat_dimid)))
      ERR(retval);
   if ((retval = nc_def_dim(ncid, "longitude",1440, &lon_dimid)))
      ERR(retval);
   if ((retval = nc_def_dim(ncid, "time",ndays, &rec_dimid)))
      ERR(retval);

   /* Define the coordinate variables. We will only define coordinate
      variables for lat and lon.  Ordinarily we would need to provide
      an array of dimension IDs for each variable's dimensions, but
      since coordinate variables only have one dimension, we can
      simply provide the address of that dimension ID (&lat_dimid) and
      similarly for (&lon_dimid). */
   if ((retval = nc_def_var(ncid, "latitude", NC_FLOAT, 1, &lat_dimid, 
			    &lat_varid)))
      ERR(retval);
   if ((retval = nc_def_var(ncid, "longitude", NC_FLOAT, 1, &lon_dimid, 
			    &lon_varid)))
      ERR(retval);

   /* Assign units attributes to coordinate variables. */
   if ((retval = nc_put_att_text(ncid, lat_varid, "units", 
				 strlen("degrees_north"), "degrees_north")))
      ERR(retval);
   if ((retval = nc_put_att_text(ncid, lon_varid, "units", 
				 strlen("degrees_east"), "degrees_east")))
      ERR(retval);

   /* The dimids array is used to pass the dimids of the dimensions of
      the netCDF variables. Both of the netCDF variables we are
      creating share the same four dimensions. In C, the
      unlimited dimension must come first on the list of dimids. */
   dimids[0] = lat_dimid;
   dimids[1] = lon_dimid;

   /* Define the netCDF variables for the pressure and temperature
    * data. */
   if ((retval = nc_def_var(ncid, "insystem", NC_FLOAT, 2, dimids, &temp_varid)))
      ERR(retval);
   if ((retval = nc_def_var(ncid, "allocated", NC_FLOAT, 2, dimids, &pres_varid)))
      ERR(retval);
       char buf[100];

    int j = snprintf(buf, 20, "%f", runtime);
   if ((retval = nc_put_att_text(ncid, pres_varid, "runtime", 
				 strlen(buf), buf)))
      ERR(retval);
   if ((retval = nc_def_var(ncid, "released", NC_FLOAT, 2, dimids, &rel_varid)))
      ERR(retval);
/*
   if ((retval = nc_put_att_text(ncid, pres_varid, "units", 
				 strlen(PRES_"units"), PRES_"units")))
      ERR(retval);
   if ((retval = nc_put_att_text(ncid, temp_varid, "units", 
				 strlen(TEMP_"units"), TEMP_"units")))
      ERR(retval);
*/

   /* End define mode. */
   if ((retval = nc_enddef(ncid)))
      ERR(retval);

   /* Write the coordinate variable data. This will put the latitudes
      and longitudes of our data grid into the netCDF file. */
   if ((retval = nc_put_var_float(ncid, lat_varid, &lats[0])))
      ERR(retval);
   if ((retval = nc_put_var_float(ncid, lon_varid, &lons[0])))
      ERR(retval);

   /* These settings tell netcdf to write one timestep of data. (The
     setting of start[0] inside the loop below tells netCDF which
     timestep to write.) */
   count[0] = 721;
   count[1] = 1440;
   start[0] = 0;
   start[1] = 0;

   /* Write the pretend data. This will write our surface pressure and
      surface temperature data. The arrays only hold one timestep worth
      of data. We will just rewrite the same data for each timestep. In
      a real application, the data would change between timesteps. */
      if ((retval = nc_put_vara_double(ncid, temp_varid, start, count,&insystem[0])))
	 ERR(retval);
      if ((retval = nc_put_vara_double(ncid, pres_varid, start, count,&allocated[0])))
	 ERR(retval);
      if ((retval = nc_put_vara_double(ncid, rel_varid, start, count,&released[0])))
	 ERR(retval);

   /* Close the file. */
   if ((retval = nc_close(ncid)))
      ERR(retval);
   
   printf("*** Output written ***\n");
   return 0;
}
