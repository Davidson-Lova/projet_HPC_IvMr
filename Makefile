#CC = g++
CC = mpic++

CFLAGS = -Wall -fPIC -O2 -mtune=native -fopenmp
LIBS = -lm -lrt  -lmpi 

#OBJS = Chronometer.o Advection_diffusion_seq.o  
OBJS = Advection_diffusion_seq.o  

%.o: %.cpp
	$(CC) $(CFLAGS) -c $<

all:	advectiondiff.exe

clean:
	@rm -fr *.o *.exe *~

advectiondiff.exe: $(OBJS)
	$(CC) $(CFLAGS) -o $@ $(OBJS) $(LIBS)
