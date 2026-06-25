The Visual Studio solution file here contains project files for all of the C
libraries.  They were originally set up to build static libraries for
compiling with Cygwin compilers.  They have been modified for their current
main use, building libraries for SerialEM, in several ways:

1) The 64-bit output files for libcfshr and libiimod are named libcfshr-64.lib
and libiimod-64.lib, respectively.  To use these for a general compilation,
the project configurations would have to be changed to rename these back to
libcfshr.lib and libiimod.lib; or the libraries could simply be renamed after
they are built.

2) A new configuration, DLL-Release, was added so that libifft could be built
as a DLL with the Intel Math Kernel Libraries statically linked.  In addition,
this configuration makes static libs (not DDLs) for libcfshr and libiimod that
are linked against libiomp5md.lib instead of vcomp.lib.

3) The Ctffind library, libctffind, needed to be a DLL in Windows to avoid
some multiply defined symbols.  SerialEM is distributed with the DLLs from the
regular IMOD build with the Intel compilers because the 32-bit version is so
much faster than the one built in Visual Studio.  The Release configurations
in this solution file produce libraries named libctffind-VCOMP.dll that depend
on vcomp.dll; they can be used for building a SerialEM that does not depend on
libiomp5md.dll.  To use this for a general compilation, the project
configuration would have to be changed to produce libctffind.dll and
libctffind.lib; or the -VCOMP files could simply be renamed after they are
built; in either case, libctffind.lib would need to be placed in buildlib.
The DLL-Release configurations produce libraries named libctffind.dll that
depend on libiomp5md.lib.  These can be used interchangeably with the
Intel-built libraries in building SerialEM against libiomp5md.lib, for testing
purposes. 
