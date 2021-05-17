#export LD_LIBRARY_PATH=/srv/djw1g16/ray_dev/sph_raytrace/lib:$LD_LIBRARY_PATH

#CONFIG = -D BENCHMARK
CONFIG =

TEST_OBJS  =   ArrayParticlePositionCoupler.o  testdata.o main.o

LIB_OBJS  =   coupling.o absorbtree.o agn_kernel.o raytracing.o

CC	 =  mpicc
CXX      =  mpicxx

#OPTIMIZE =  -O2 -Wall -std=c++11 -flto -funroll-loops  -march=native
#OPTIMIZE =  -pg -O2 -Wall -std=c++11 -flto -funroll-loops  -march=native
OPTIMIZE =  -pg -std=c++11 -flto

MPICHLIB =  -L/usr/lib64/openmpi/lib

LIBS   = -g $(MPICHLIB) 
INC = -I../src/

EXEC = RAYTEST
LIBOUT = sph_raytrace


test: $(TEST_OBJS)
	$(CXX) $(OPTIMIZE) $(TEST_OBJS) $(LIBS) -L../lib/ -l $(LIBOUT) $(INC)  -o ../$(EXEC)

$(TEST_OBJS): %.o: %.cpp
	$(CXX) $(OPTIMIZE) $(CONFIG) $(INC) -c $< -o $@


clean:
	rm -f $(TEST_OBJS) ../$(EXEC)
