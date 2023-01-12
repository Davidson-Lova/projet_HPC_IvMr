#CC = g++
CC = mpic++

CFLAGS = -Wall -fPIC -O2 -mtune=native -fopenmp
LIBS = -lm -lrt  -lmpi 

#OBJS = Chronometer.o Advection_diffusion_seq.o  
OBJS = Advection_diffusion_seq.o  
OBJS_omp = Advection_diffusion_omp.o  
OBJS_mpi1 = Advection_diffusion_mpiv1.o  
OBJS_mpi2 = Advection_diffusion_mpiv2.o  

%.o: %.cpp
	$(CC) $(CFLAGS) -c $<

all:	advectiondiff.exe advectiondiff_omp.exe advectiondiff_mpiv1.exe advectiondiff_mpiv2.exe

clean:
	@rm -fr *.o *.exe *~ res.txt Sort* res_Sort*

advectiondiff.exe: $(OBJS)
	$(CC) $(CFLAGS) -o $@ $(OBJS) $(LIBS)

advectiondiff_omp.exe: $(OBJS_omp)
	$(CC) $(CFLAGS) -o $@ $(OBJS_omp) $(LIBS)

advectiondiff_mpiv1.exe: $(OBJS_mpi1)
	$(CC) $(CFLAGS) -o $@ $(OBJS_mpi1) $(LIBS)

advectiondiff_mpiv2.exe: $(OBJS_mpi2)
	$(CC) $(CFLAGS) -o $@ $(OBJS_mpi2) $(LIBS)
