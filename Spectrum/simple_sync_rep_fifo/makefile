COMPILER =   	g++
EXECUTABLE = 	simple_sync_rep_fifo
LIBS += 	-lspcm_linux -lpthread -lrt

OBJECTS = 	simple_sync_rep_fifo.o\
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

