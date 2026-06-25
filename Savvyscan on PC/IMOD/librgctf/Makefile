#
# Makefile for libctffind and test program
#
# $Id$
#

include ../configure

O = $(OBJEXT)
E = $(BINEXT)
LIBCTF   = libctffind.$(DLLEXT)
CXXFLAGS = $(LIBCFLAGS) -DEXPORT_FOR_DLL


PROGS     = testctffind$(E)

IILIBS    =  $(LIIMOD) $(LCFSHR) $(LIMXML) $(TIFFLIBS)

CTFOBJS = angles_and_shifts.$(O) brute_force_search.$(O) ctf.$(O) curve.$(O) \
 conjugate_gradient.$(O) empirical_distribution.$(O) functions.$(O) globals.$(O) \
 image_trim.$(O) matrix.$(O) ctffind.$(O) randomnumbergenerator.$(O) va04.$(O)


all : $(LIBCTF) $(PROGS)

functions.$(O) : functions.cpp
	$(CXX) $(LIBCFLAGS) $(OPENMP) -c functions.cpp

brute_force_search.$(O) : brute_force_search.cpp
	$(CXX) $(LIBCFLAGS) $(OPENMP) -c brute_force_search.cpp

$(LIBCTF) : $(CTFOBJS)
	$(FIND) . -type f -name $(LIBCTF) -exec rm "{}" \;
	$(MAKELIBCPP)$(LIBCTF) $(LIBCTFFINDDEF) $(CTFOBJS)
	$(RANLIB) $(LIBCTF)
	$(CP) libctffind.$(LIBEXT) $(LIBDIR)

testctffind$(E) : testctffind.$(O) ../include/ctffind.h $(LIBCTF)
	$(CXX) -o testctffind$(E) testctffind.$(O) $(LDFLAGS) $(LCTFFIND) $(IILIBS) \
	$(LIFFT) $(CLIBS)
testctffind.$(O): testctffind.cpp
	$(CXX) $(LIBCFLAGS) -c testctffind.cpp -o testctffind.$(O)

install : $(LIBCTF)
	$(CP) $(LIBCTF) $(INSTLIBDIR)
	$(STRIP) $(PROGS)
	$(CP) $(PROGS) $(INSTBINDIR)

clean : 
	$(FIND) . -type f -name "*.$(OBJEXT)" -exec rm "{}" \;
	$(FIND) . -name "*~" -exec rm "{}" \;
	$(FIND) . -type f -name "*.exe" -exec rm "{}" \;
	$(FIND) . -type f -name "*.so" -exec rm "{}" \;
	$(FIND) . -type f -name "*.a" -exec rm "{}" \;
	$(FIND) . -type f -name "*.*lib" -exec rm "{}" \;
	$(FIND) . -type f -name "*.dll" -exec rm "{}" \;
	$(FIND) . -type f -name "*.exp" -exec rm "{}" \;
	-rm -f $(PROGS)
