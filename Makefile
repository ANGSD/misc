FLAGS=-O3


CFLAGS += $(FLAGS) -std=gnu99
CXXFLAGS += $(FLAGS)

PRG=cmpFasta getinsertsize readplink safsubsampler simPileup


all: $(PRG)

.PHONY: clean

CSRC = $(wildcard *.c) 
CXXSRC = $(wildcard *.cpp)
OBJ = $(CSRC:.c=.o) $(CXXSRC:.cpp=.o)

-include $(OBJ:.o=.d)

%.o: %.c
	$(CC) -c  $(CFLAGS) $*.c
	$(CC) -MM $(CFLAGS) $*.c >$*.d
%.o: %.cpp
	$(CXX) -c  $(CXXFLAGS) $*.cpp
	$(CXX) -MM $(CXXFLAGS) $*.cpp >$*.d

cmpFastaOBJ = razf.o faidx.o cmpFasta.o
cmpFasta: $(cmpFastaOBJ)
	$(CXX) $(FLAGS)  -o cmpFasta $(cmpFastaOBJ) -lz 

getinsertsize: getinsertsize.o
	$(CXX) $(FLAGS)  -o getinsertsize getinsertsize.o -lz

simPileup: simPileup.o
	$(CXX) $(FLAGS)  -o simPileup simPileup.o -lz

safsubsampler: safsubsampler.o
	$(CXX) $(FLAGS)  -o safsubsampler safsubsampler.o -lz

readplink: readplink.c
	$(CC) $(FLAGS)  -o readplink readplink.c -lz -D__WITH_MAIN__

clean:
	rm  -f *.o *.d getinsertsize cmpFasta readplink safsubsampler simPileup *~
