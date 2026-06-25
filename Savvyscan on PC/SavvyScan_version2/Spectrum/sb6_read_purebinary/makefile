COMPILER =   	g++
EXECUTABLE = 	sb6_read_purebinary
LIBS += 	-lpthread -lrt

OBJECTS = 	main.o\
            sb6_read_purebinary.o

all: $(EXECUTABLE)

$(EXECUTABLE):	$(OBJECTS)
		$(COMPILER) $(CFLAGS) -o $(EXECUTABLE) $(OBJECTS) $(LIBS) 

%.o:		%.cpp
		$(COMPILER) $(CFLAGS) -o $*.o -c $*.cpp

install:

clean:
	rm -f $(EXECUTABLE)
	rm -f $(OBJECTS)

