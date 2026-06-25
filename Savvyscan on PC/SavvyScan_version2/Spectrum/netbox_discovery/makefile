COMPILER =   	g++
EXECUTABLE = 	netbox_discovery
LIBS += 	-lspcm_linux -lpthread -lrt

OBJECTS = 	netbox_discovery.o\
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

