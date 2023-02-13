COMPILER =   g++
EXECUTABLE = rec_fifo_aba_thread
LIBS += -lspcm_linux -lpthread -lrt

OBJECTS = \
    rec_fifo_aba_thread.o\
	../common/spcm_lib_card.o\
	../common/spcm_lib_data.o\
	../common/spcm_lib_thread.o\
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

