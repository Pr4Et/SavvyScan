COMPILER =   g++
EXECUTABLE = rep_sequence
LIBS += -lspcm_linux -lpthread -lrt
CFLAGS += -DSPCM_NAMESPACE=spcmdrv

OBJECTS = \
    rep_sequence.o\
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

