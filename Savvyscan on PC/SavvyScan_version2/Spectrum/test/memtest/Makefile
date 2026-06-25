CFLAGS =     -fPIC -Wall -O1 -DISP_SUPPORT_M2IKERNEL  -DSPCM_NAMESPACE=spcm_memtest -DNDEBUG

COMPILER =   g++
EXECUTABLE = spcm_memtest

#IMPORTANT: add files with .o ending, otherwise linux will produce run-time errors when library is used

OBJECTS =       main.o \
				../../common/ostools/spcm_ostools_linux.o	


LIBS +=       -lrt -lstdc++ -lpthread -lspcm_linux

all: $(EXECUTABLE)

$(EXECUTABLE):	$(OBJECTS)
		$(COMPILER) $(CFLAGS) $(GCC3VISIBILITYFLAGS) -o $(EXECUTABLE) $(OBJECTS) $(LIBS)

%.o:		%.cpp
		$(COMPILER) $(CFLAGS) -o $*.o -c $*.cpp

install:

clean:
	rm -f $(EXECUTABLE)
	rm -f $(OBJECTS)

