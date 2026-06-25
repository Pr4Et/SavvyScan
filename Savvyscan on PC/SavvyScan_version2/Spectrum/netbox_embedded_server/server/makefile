COMPILER =   	g++
EXECUTABLE = 	netbox_embedded_server
LIBS += 	-lspcm_linux -lpthread -lrt
CFLAGS +=   -DSPCM_NAMESPACE=server

OBJECTS = 	server_main.o\
		../../common/ostools/spcm_ostools_linux.o \
		../../common/ostools/spcm_network_winLin.o \

all: $(EXECUTABLE)

$(EXECUTABLE):	$(OBJECTS)
		$(COMPILER) $(CFLAGS) -o $(EXECUTABLE) $(OBJECTS) $(LIBS) 

%.o:		%.cpp
		$(COMPILER) $(CFLAGS) -o $*.o -c $*.cpp

install:

clean:
	rm -f $(EXECUTABLE)
	rm -f $(OBJECTS)

