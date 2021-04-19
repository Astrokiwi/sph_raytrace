EXEC   = RAYTEST

OBJS  =   coupling.o testdata.o absorbtree.o agn_kernel.o raytracing.o main.o

CC	 =  mpicc
CXX      =  mpicxx

#OPTIMIZE =  -O2 -Wall -std=c++11 -flto -funroll-loops  -march=native
#OPTIMIZE =  -pg -O2 -Wall -std=c++11 -flto -funroll-loops  -march=native
OPTIMIZE =  -pg -std=c++11 -flto
MPICHLIB =  -L/usr/lib64/openmpi/lib

LIBS   = -g $(MPICHLIB) 

$(EXEC): $(OBJS)
	$(CXX) $(OPTIMIZE) $(OBJS) $(LIBS) -o ../$(EXEC)

$(OBJS): %.o: %.cpp
	$(CXX) $(OPTIMIZE) -c $< -o $@


clean:
	rm -f $(OBJS) $(FOBJS) ../$(EXEC) *.o
