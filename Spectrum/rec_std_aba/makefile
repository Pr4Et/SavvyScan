COMPILER =   g++
EXECUTABLE = rec_std_aba
LIBS += -lspcm_linux -lpthread -lrt

OBJECTS = \
    rec_std_aba.o\
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
	rm -f $(EXECUTABLE)
	rm -f $(OBJECTS)

