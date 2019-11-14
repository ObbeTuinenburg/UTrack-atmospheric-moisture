LIBS=-lnetcdf
clean:
	rm *.o

all:
	gcc -fopenmp -O3 -o run_recycling main.c interpolation.c read.c write.c load_meteo.c eulerian.c lagrangian.c trajectory.c -lm -lnetcdf 
