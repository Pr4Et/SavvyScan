include ../../configure

FFLAGS = $(LIBFFLAGS)

O        = $(OBJEXT)
OBJS = xflincom.$(O) \
 xfrdall.$(O) xfread.$(O) xfwrite.$(O) readdistortions.$(O) findxf.$(O) \
 writedistortions.$(O)

all: $(OBJS)

clean:
	$(FIND) . -type f -name "*.$(OBJEXT)" -exec rm "{}" \;
	$(FIND) . -type f -name "*~" -exec rm "{}" \;
	$(FIND) . -type f -name "#*#" -exec rm "{}" \;
