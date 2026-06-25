COMPILER =   	g++
EXECUTABLE = 	closed_loop_ad_da
LIBS += 	-lspcm_linux -lpthread -lrt
CFLAGS += -O2

OBJECTS = 	closed_loop_ad_da.o\
		../../common/ostools/spcm_ostools_linux.o

all: $(EXECUTABLE)

$(EXECUTABLE):	$(OBJECTS)
		$(COMPILER) $(CFLAGS) -o $(EXECUTABLE) $(OBJECTS) $(LIBS) 

%.o:		%.cpp
		$(COMPILER) $(CFLAGS) -o $*.o -c $*.cpp

install:

clean:
	rm -f $(EXECUTABLE)
	rm -f $(OBJECTS)

