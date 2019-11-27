
struct meteoday {
	float* ef;
	double* efadd;
	double* efscale;
	float* nf;
	double* nfadd;
	double* nfscale;
	float* pw;
	double* pwadd;
	double* pwscale;
	float* q;
	double* qadd;
	double* qscale;
	float* t;
	double* tadd;
	double* tscale;
	float* u;
	double* uadd;
	double* uscale;
	float* v;
	double* vadd;
	double* vscale;
	float* w;
	double* wadd;
	double* wscale;
	float* E;
	double* Eadd;
	double* Escale;
	float* P;
	double* Padd;
	double* Pscale;
	float* lat;
	float* lon;
	float* lev;
	int   day;
	int   month;
	int   year;
};

struct simmeteo {
	struct meteoday* meteo[2];
};

struct state {
	float u;
	float v;
	float w;
	float E;
	float P;
	float PW;
	float nf;
	float ef;
};

struct loc {
	float lat;
	float lon;
	float lev;
	float theta;
	float time;
	int day;
	int curlevidx;
	int curlevidx2;
};

struct simulation_settings {
	int simulation_number;
	float mixing;
	int interpolation;
	int daynum;
	int nday;
	int* parcels_in_system;
	int num_release;
	int release_days;
	int release_parcel;
	int hour;
	float runtime;

	// Allocated memory for parcels during tracking
	float* particle_lat;
	float* particle_lon;
	float* particle_lev;
	float* particle_theta;
	float* particle_original;
	float* particle_present;
	float* particle_curtime;
	
	// Grid cell size
	float* grid_size_NW;
	float* grid_size_EW;

	// Input
	int* release_input;
	char inputfile[100];

	// Output
	char outputfile[100];
	double* emitted_out;
	double* budget;
	double* insystem_out;
	double* insystem;
	double* allocated_out;

	struct simulation_settings *next;

};



int readdata(float* array, char* s, const char *path, double* add, double* offset, float* lats, float* lons, float* levs, float* time);
int read2d(float* array, char* s, const char *path, double* add, double* scale, float* lats, float* lons, float* time);
int readgldas(float* array, char* s, const char *path, float* lats, float* lons, double* time);
int readcellnum(int* array, char* s, const char *path, float* lats, float* lons);
int readisogsm(float* array, char* s, const char *path, float* lats, float* lons, float* time);
int readirrigation(double* array, char* s, const char *path, float* lats, float* lons);
int readforest(float* array, char* s, const char *path, float* lats, float* lons);
int readcat(float* array, char* s, const char *path, float* lats, float* lons);
int readbasins(double* array, char* s, const char *path, float* lats, float* lons);
int readisomean(float* arraye, float* arrayp, const char *path);

struct state interpolate(struct loc curloc,struct meteoday* meteo, struct simulation_settings s);
int get_lsm(short*array,double*add,double*scale);
int get_era5_hour(int year,int month,int day,int hour,struct meteoday meteo);
int get_era5_3d(int year,int month,int day,int hour,struct meteoday meteo);
int get_era5_3dt(int year,int month,int day,int hour,struct meteoday meteo);
int get_era5_3dq(int year,int month,int day,int hour,struct meteoday meteo);
int get_era5_3dw(int year,int month,int day,int hour,struct meteoday meteo);
int get_era5_3duv(int year,int month,int day,int hour,struct meteoday meteo);

int get_input(int *array,char* path);

int write_output(double *released,double *insystem,double *allocated,float *lats,float *lons,int ndays,char* outfile, float runtime);
int free_meteo_hour(struct meteoday* meteo);
int load_meteo_hour(int curyear,int curmonth,int curday,int curhour,struct meteoday* meteo);
int lagrangian_simulation(struct meteoday* meteo,struct simulation_settings simulation);
int trajectory(struct meteoday* meteo, struct simulation_settings s,int particle_number);
int eulerian_simulation(struct meteoday* meteo,struct simulation_settings simulation);
float get_starting_height(struct meteoday* meteo,struct simulation_settings simulation,int lat_idx,int lon_idx);
float get_potential_temp(struct meteoday* meteo,int i,int j,float lev);
int get_3d_interpolation(struct state* outstate,struct meteoday* meteo,int curlatidx,int curlonidx,int curlevidx,float latfrac,float lonfrac,float levfrac,float timefrac);
int get_2d_interpolation(struct state* outstate,struct meteoday* meteo,int curlatidx,int curlonidx,float latfrac,float lonfrac,float timefrac);
struct state gridcell_state(struct loc curloc,struct meteoday* meteo, struct simulation_settings s);
int eul2layers(struct meteoday* meteo,struct simulation_settings simulation);
