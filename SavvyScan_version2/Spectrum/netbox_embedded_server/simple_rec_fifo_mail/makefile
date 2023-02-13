COMPILER =   	g++
EXECUTABLE = 	simple_rec_fifo_mail
LIBS += 	-lspcm_linux -lpthread -lrt -lcurl

OBJECTS = 	simple_rec_fifo_mail.o\
		../../common/ostools/spcm_ostools_linux.o\
		../../common//spcm_netbox_sendmail.o\

all: $(EXECUTABLE)

$(EXECUTABLE):	$(OBJECTS)
		$(COMPILER) $(CFLAGS) -o $(EXECUTABLE) $(OBJECTS) $(LIBS) 

%.o:		%.cpp
		$(COMPILER) $(CFLAGS) -o $*.o -c $*.cpp

install:

clean:
	rm -f $(EXECUTABLE)
	rm -f $(OBJECTS)

