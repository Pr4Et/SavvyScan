COMPILER =   	g++
EXECUTABLE = 	spectrum_simple_template
LIBS += 	-lspcm_linux -lpthread -lrt

OBJECTS = 	spectrum_simple_template.o\
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

