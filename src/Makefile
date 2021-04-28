#export LD_LIBRARY_PATH=/srv/djw1g16/ray_dev/sph_raytrace/lib:$LD_LIBRARY_PATH

TEST_OBJS  =   ArrayParticlePositionCoupler.o  testdata.o main.o

LIB_OBJS  =   coupling.o absorbtree.o agn_kernel.o raytracing.o

CC	 =  mpicc
CXX      =  mpicxx

#OPTIMIZE =  -O2 -Wall -std=c++11 -flto -funroll-loops  -march=native
#OPTIMIZE =  -pg -O2 -Wall -std=c++11 -flto -funroll-loops  -march=native
OPTIMIZE =  -pg -std=c++11 -flto
LIB_FLAGS = -shared $(OPTIMIZE) -fPIC

MPICHLIB =  -L/usr/lib64/openmpi/lib

LIBS   = -g $(MPICHLIB) 

EXEC = RAYTEST
LIBOUT = sph_raytrace.so

lib: $(LIB_OBJS)
	$(CXX) -shared $(LIB_FLAGS) $(LIB_OBJS) $(LIBS) -o ../lib/$(LIBOUT)

test: $(TEST_OBJS)
	$(CXX) $(OPTIMIZE) $(TEST_OBJS) $(LIBS) $(LIBOUT) -o ../$(EXEC)

$(TEST_OBJS): %.o: %.cpp
	$(CXX) $(OPTIMIZE) -c $< -o $@

$(LIB_OBJS): %.o: %.cpp
	$(CXX) $(OPTIMIZE) -c $< -o $@


clean:
	rm -f $(TEST_OBJS) ../$(EXEC) *.o $(LIBOUT).so
