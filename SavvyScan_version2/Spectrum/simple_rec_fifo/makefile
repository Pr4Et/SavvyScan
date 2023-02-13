COMPILER =   	g++
EXECUTABLE = 	simple_rec_fifo
LIBS += 	-lspcm_linux -lpthread -lrt

OBJECTS = 	simple_rec_fifo.o\
		../common/ostools/spcm_ostools_linux.o

all: $(EXECUTABLE)

$(EXECUTABLE):	$(OBJECTS)
		$(COMPILER) $(CFLAGS) -o $(EXECUTABLE) $(OBJECTS) $(LIBS) 

%.o:		%.cpp
		$(COMPILER) $(CFLAGS) -o $*.o -c $*.cpp

install:

clean:
	rm -f $(EXECUTABLE)
	rm -f $(OBJECTS)

