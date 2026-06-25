#
# Makefile for libimxml
#
# $Id$
#
include ../configure

LIBI   = libimxml.$(DLLEXT)
O     = $(OBJEXT)
CFLAGS = $(LIBCFLAGS) $(MXMLFLAGS)

MXMLOBJS = mxml-attr.$(O) mxml-file.$(O) mxml-index.$(O) mxml-private.$(O) mxml-set.$(O) \
	mxml-entity.$(O) mxml-get.$(O) mxml-node.$(O) mxml-search.$(O) mxml-string.$(O)

all : $(LIBI)

$(LIBI): $(MXMLOBJS)
	$(FIND) . -type f -name $(LIBI) -exec rm "{}" \;
	$(MAKELIB)$(LIBI) $(LIBIMXMLDEF) $(MXMLOBJS)
	$(RANLIB) $(LIBI)
	$(CP) libimxml.$(LIBEXT) $(LIBDIR)

install : $(LIBI)
	$(CP) $(LIBI) $(INSTLIBDIR)

$(MXMLOBJS) : config.h ../include/mxml.h mxml-private.h

clean:
	$(FIND) . -type f -name "*.$(OBJEXT)" -exec rm "{}" \;
	$(FIND) . -type f -name "*~" -exec rm "{}" \;
	$(FIND) . -type f -name "#*#" -exec rm "{}" \;
