#
# Makefile for mrc programs.
#
# $Id$
#

include ../configure

O = $(OBJEXT)
E = $(BINEXT)


PROGS     = raw2mrc$(E) tif2mrc$(E) mrcbyte$(E) mrcx$(E) preNID$(E) alignframes$(E) \
 mrctaper$(E) dm3props$(E) nad_eed_3d$(E) ctfphaseflip$(E) preNAD$(E) $(FRAMEGPU)

IILIBS    =  $(LIIMOD) $(LCFSHR) $(LIMXML) $(TIFFLIBS)

FRAMEOBJS = alignframes.$(O) framealign.$(O) CorrectDefects.$(O) frameutil.$(O) \
	$(FRAMEGPUOBJ)

all : $(PROGS)

##############################################################################
# make programs.
#

raw2mrc$(E) : raw2mrc.$(O) 
	$(CC) -o raw2mrc raw2mrc.$(O)  $(LDFLAGS) $(IILIBS) $(CLIBS)

tif2mrc$(E) : tif2mrc.$(O) tiff.$(O) 
	$(CC) -o tif2mrc tif2mrc.$(O) tiff.$(O) $(LDFLAGS) $(IILIBS) $(CLIBS)

mrcx$(E) : mrcx.$(O)
	$(CC) -o mrcx mrcx.$(O) $(LDFLAGS) $(IILIBS) $(CLIBS)

mrcbyte$(E) : mrcbyte.$(O) 
	$(CC) -o mrcbyte mrcbyte.$(O) $(LDFLAGS) $(IILIBS)  $(CLIBS)

mrctaper$(E) : mrctaper.$(O) 
	$(CC) -o mrctaper mrctaper.$(O) $(LDFLAGS) $(IILIBS)  $(CLIBS)

dm3props$(E) : dm3props.$(O) 
	$(CC) -o dm3props dm3props.$(O) $(LDFLAGS) $(IILIBS)  $(CLIBS)

nad_eed_3d$(E) : nad_eed_3d.$(O) nrutil.$(O)
	$(CC) -o nad_eed_3d nad_eed_3d.$(O) nrutil.$(O) $(LDFLAGS) $(IILIBS) \
	$(CLIBS)

ctfphaseflip$(E) : ctfphaseflip.$(O) ctfutils.$(O) gpuctf.$(O) frameutil.$(O)
	$(CXX) -o ctfphaseflip ctfphaseflip.$(O) ctfutils.$(O) gpuctf.$(O) \
	frameutil.$(O) $(LDFLAGS) $(CUDALIBS) $(IILIBS) $(LIFFT) $(CLIBS)

ctfutils.$(O) : ../ctfplotter/ctfutils.cpp
	$(CXX) -c -o ctfutils.$(O) $(CXXFLAGS) ../ctfplotter/ctfutils.cpp

gpuctf.$(O) : gpuctf.cu frameutil.h
	$(NVCC_TOOL_ENV) nvcc -c $(NVCC_FLAGS) -I../include gpuctf.cu

preNAD$(E) : preNAD.$(O) ctfutils.$(O) recline.$(O)
	$(LAPACKLINK) -o preNAD preNAD.$(O) ctfutils.$(O) recline.$(O) $(LDFLAGS) \
	$(IILIBS) $(LAPACKLIBS) $(LAPACK_FLIBS) $(FORT_CPPLIBS) $(CLIBS)

preNID$(E) : preNID.$(O) ctfutils.$(O) recline.$(O)
	$(LAPACKLINK) -o preNID preNID.$(O)  ctfutils.$(O) recline.$(O) $(LDFLAGS) \
	$(IILIBS) $(LAPACKLIBS) $(LAPACK_FLIBS) $(FORT_CPPLIBS) $(CLIBS)


alignframes$(E) : $(FRAMEOBJS)
	$(CXX) -o alignframes $(FRAMEOBJS) $(LDFRAMEFLAGS) $(IILIBS) \
	$(LIFFT) $(CLIBS)

gpuframe.$(O) : gpuframe.cu gpuframe.h framealign.h frameutil.h
	$(NVCC_TOOL_ENV) nvcc -c $(NVCC_FLAGS) -I../include gpuframe.cu

alignframes.$(O) : alignframes.cpp  framealign.h gpuframe.h
	$(CXX) -c -o alignframes.$(O) $(CFRAMEFLAGS) alignframes.cpp

framealign.$(O) : framealign.cpp  framealign.h gpuframe.h
	$(CXX) -c -o framealign.$(O) $(OPENMP) $(CFRAMEFLAGS) framealign.cpp

frameutil.$(O) : frameutil.cpp  frameutil.h
	$(CXX) -c -o frameutil.$(O) $(CFRAMEFLAGS) frameutil.cpp

frameutildll.$(O) : frameutil.cpp frameutil.h
	nvcc -c $(NVCC_FLAGS) -I../include frameutil.cpp -o frameutildll.obj

CorrectDefects.$(O) : ../clip/CorrectDefects.cpp
	$(CXX) -c -o CorrectDefects.$(O) $(CFRAMEFLAGS) ../clip/CorrectDefects.cpp

ShrMemClient.$(O) : ShrMemClient.cpp ShrMemClient.h shrmemframe.h
	$(CXX) -c -o ShrMemClient.$(O) $(CFRAMEFLAGS) ShrMemClient.cpp

FrameGPU.dll : gpuframe.obj  frameutildll.obj
	$(FIND) . -type f -name FrameGPU.dll -exec rm "{}" \;
	$(MAKELIB)FrameGPU.dll gpuframe.obj frameutildll.obj \
	$(FGPULIBS) $(CUDALIBS)

shrmemframe$(E) : framealign.$(O) CorrectDefects.$(O) frameutil.$(O) shrmemframe.$(O)
	$(CXX) -o shrmemframe framealign.$(O) CorrectDefects.$(O) frameutil.$(O) \
	shrmemframe.$(O) $(LDFRAMEFLAGS) $(IILIBS) $(LIFFT) $(CLIBS)	

shrmemframe.$(O) : shrmemframe.cpp  shrmemframe.h framealign.h gpuframe.h
	$(CXX) -c -o shrmemframe.$(O) $(CFRAMEFLAGS) shrmemframe.cpp

##############################################################################
# install programs

install : $(PROGS)
	$(STRIP) $(PROGS)
	$(CP) $(PROGS) $(INSTBINDIR)

clean : 
	$(FIND) . -type f -name "*.$(OBJEXT)" -exec rm "{}" \;
	$(FIND) . -name "*~" -exec rm "{}" \;
	$(FIND) . -type f -name "core*" -exec rm "{}" \;
	$(FIND) . -type f -name "*.exe" -exec rm "{}" \;
	-rm -f $(PROGS)

cleanframe :
	-rm -f $(FRAMEOBJS) alignframes$(E) gpuframe.$(O) frameutildll.$(O) \
	shrmemframe$(E) shrmemframe.$(O) FrameGPU.dll ShrMemClient.$(O)

##############################################################################
#
# Not in bl3dfs distribution:
# older, test, and unsupported programs to make.
#

mrclog$(E) : mrclog.$(O) 
	$(CC) -o mrclog mrclog.$(O) -L$(LIBDIR) $(LIBS)
