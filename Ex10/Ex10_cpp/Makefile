randomdir := ../GenRand

MPICXX = mpicxx
CFLAGS = -Wall -O3 --std=c++11 -I$(randomdir)

EXECUTABLES = main1.exe main2.exe

all: 
		make $(EXECUTABLES)

%.exe : %.o random.o
	$(MPICXX) $^ -o $@

%.o : %.cpp
	$(MPICXX) -c $< -o $@ $(CFLAGS)

random.o : $(randomdir)/random.cpp $(randomdir)/random.h
	$(MPICXX) -c $(randomdir)/random.cpp -o random.o $(CFLAGS)

clean :
	rm -f *.o *.exe seed.out 

remove :
	rm -f *.txt *.dat

