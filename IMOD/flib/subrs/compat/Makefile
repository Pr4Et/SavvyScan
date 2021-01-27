include ../../configure

CFLAGS = $(LIBCFLAGS)
FFLAGS = $(LIBFFLAGS)

OBJS = degtrig.$(OBJEXT) $(DATEOBJ)

all: $(OBJS)

clean:
	$(FIND) . -type f -name "*.$(OBJEXT)" -exec rm "{}" \;
	$(FIND) . -type f -name "*~" -exec rm "{}" \;
	$(FIND) . -type f -name "#*#" -exec rm "{}" \;
