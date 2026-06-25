COMPILER =   g++
EXECUTABLE = rep_fifo_multi
LIBS += -lspcm_linux -lpthread -lrt

OBJECTS = \
	rep_fifo_multi.o\
	../common/spcm_lib_card.o\
	../common/spcm_lib_data.o\
	../common/ostools/spcm_ostools_linux.o

all: $(EXECUTABLE)

$(EXECUTABLE):	$(OBJECTS)
		$(COMPILER) $(CFLAGS) -o $(EXECUTABLE) $(OBJECTS) $(LIBS)

%.o:		%.cpp
		$(COMPILER) $(CFLAGS) -o $*.o -c $*.cpp

install:

clean:
	rm -f $(OBJECTS)
	rm -f $(EXECUTABLE)

