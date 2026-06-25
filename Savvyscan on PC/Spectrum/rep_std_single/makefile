COMPILER =   g++
EXECUTABLE = rep_std_single
LIBS += -lspcm_linux -lpthread -lrt

OBJECTS = \
    rep_std_single.o\
    	../common/ostools/spcm_ostools_linux.o\
	../common/spcm_lib_card.o\
	../common/spcm_lib_data.o

all: $(EXECUTABLE)

$(EXECUTABLE):	$(OBJECTS)
		$(COMPILER) $(CFLAGS) -o $(EXECUTABLE) $(OBJECTS) $(LIBS)

%.o:		%.cpp
		$(COMPILER) $(CFLAGS) -o $*.o -c $*.cpp

install:

clean:
	rm -f $(EXECUTABLE)
	rm -f $(OBJECTS)

