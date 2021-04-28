TEST_OBJS  =   coupling.o testdata.o absorbtree.o agn_kernel.o raytracing.o main.o

CC	 =  mpicc
CXX      =  mpicxx

#OPTIMIZE =  -O2 -Wall -std=c++11 -flto -funroll-loops  -march=native
#OPTIMIZE =  -pg -O2 -Wall -std=c++11 -flto -funroll-loops  -march=native
OPTIMIZE =  -pg -std=c++11 -flto
MPICHLIB =  -L/usr/lib64/openmpi/lib

LIBS   = -g $(MPICHLIB) 

EXEC = RAYTEST

raytest: $(TEST_OBJS)
	$(CXX) $(OPTIMIZE) $(TEST_OBJS) $(LIBS) -o ../$(EXEC)

$(TEST_OBJS): %.o: %.cpp
	$(CXX) $(OPTIMIZE) -c $< -o $@


clean:
	rm -f $(TEST_OBJS) ../$(EXEC) *.o
