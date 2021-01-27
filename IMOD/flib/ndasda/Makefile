include ../configure
O        = $(OBJEXT)
E        = $(BINEXT)
LDFLG    = $(LDFFLAGS)

# DNM 3/11/01: renamed MODIMLIBS to GRAPHLIBS

PROGS = nda$(E) mtk$(E) sda$(E)

NDAOBJS = nda.$(O) ndmodel.$(O) graphden.$(O) angdist.$(O) ndrandom.$(O) \
 ndsubs.$(O) ndmtkrand.$(O) ndmtksubs.$(O) manipgraphs.$(O)

MTKOBJS = mtk.$(O) mtkmodel.$(O) closedist.$(O) mtksubs.$(O) mtkrandom.$(O)\
 mtkransubs.$(O) ndmtksubs.$(O) ndmtkrand.$(O) bundledist.$(O) \
 manipgraphs.$(O) mtkvars.$(O)

SDAOBJS = sda.$(O) sdmodel.$(O) sddist.$(O) sd3dsubs.$(O) surfden.$(O) \
 kerneldens.$(O) find_points_area.$(O) ndmtksubs.$(O) ndmtkrand.$(O) \
 ndsubs.$(O) sdrandom.$(O) manipgraphs.$(O) sdlength.$(O)

INCOBJS = sda.$(O) sdmodel.$(O) sddist.$(O) surfden.$(O) kerneldens.$(O) \
 find_points_area.$(O)



#CFLAGS = -O -I$(INCDIR)
#FFLAGS = -O -I$(INCDIR) 

LIBS = $(LWMOD) $(LIMOD) $(GRAPHLIBS) $(FORTCLIBS)

all: $(PROGS)

nda$(E): $(NDAOBJS)
	$(FC)  -o nda $(NDAOBJS) $(LDFLG) $(LIBS)

mtk$(E): $(MTKOBJS)
	$(FC)  -o mtk $(MTKOBJS) $(LDFLG) $(LIBS)

sda$(E): $(SDAOBJS)
	$(FC)  -o sda $(SDAOBJS) $(LDFLG) $(LIBS)

ndmtksubs.$(O): scrnvars.mod

mtkvars.mod: mtkvars.f90
	$(FC) -c $(LIBFFLAGS) mtkvars.f90

scrnvars.mod: ../subrs/graphics/screenpak.f90
	$(FC) -c $(LIBFFLAGS) ../subrs/graphics/screenpak.f90

$(INCOBJS): sda.inc
mtkmodel.$(O): $(INCDIR)/model.inc mtkvars.mod
mtkransubs.$(O): $(INCDIR)/model.inc mtkvars.mod
ndmodel.$(O): $(INCDIR)/model.inc
sdmodel.$(O): $(INCDIR)/model.inc
closedist.$(O): mtkvars.mod
mtksubs.$(O): mtkvars.mod
mtkrandom.$(O): mtkvars.mod
mtk.$(O): mtkvars.mod


install: $(PROGS)
	$(STRIP) $(PROGS)
	$(CP) $(PROGS) $(INSTBINDIR)

clean:
	-\rm -f *.$(OBJEXT) *.exe *~ nda mtk sda *.mod


