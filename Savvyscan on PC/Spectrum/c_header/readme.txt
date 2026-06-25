**************************************************************************
c_header directory                                       (c) Spectrum GmbH
**************************************************************************

The directory contains all header and library files for all Spectrum
drivers for ISA, PCI, PCI-X, PCIe, cPCI and PXI cards.

**************************************************************************



Common header files used by all drivers
---------------------------------------

dlltyp.h: definitions common for all Spectrum drivers and card types. This
		header tries to examine the type of compiler and then defines common
		data types that have the same length under all compilers and all
		operating systems.

regs.h: software register and constants definition for all Spectrum
		drivers.

spcerr.h: error codes of all Spectrum drivers. Until may 2004 this file was
		errors.h. Name has been changed because errors.h has been already in
		use by windows.



Library and Header files of driver for ISA/PCI/MI/MC/MX cards
-------------------------------------------------------------

spcioctl.inc:	linux include file to access driver functions via kernel
		calls. Is needed by all linux based programs that access one of the
		ISA/PCI/MI/MC/MX Spectrum cards

errors.h: former error file. Thsi file is just included because of
		compatibility reasons with old projects. Please use spcerr.h

spectrum.h: header file that contains all the prototypes of the driver
		functions

spectrum.lib:	library file for Microsoft Visual C++ for the spectrum
		driver DLL. Calling type is c-call.

SpcStdNT.lib: library file for other compilers for the spectrum
		driver DLL. Calling type is stdcall.

spclib_bcc.lib: library for Borland C++ Builder for the spectrum
		driver DLL.



Library and Header files of driver for SPCM driver based cards
-------------------------------------------------------------

spcm_drv.h: header file that contains all the prototypes of the
		driver functions of the spcm driver

spcm_win32_msvcpp.lib: library file for the Microsoft Visual C++
		compiler. Calling type is stdcall.

spcm_win32_bcppb.lib: library file for the Borland C++ Builder
		compiler

spcm_win32_cvi.lib: library file for National Instruments
		LabWindows/CVI compiler

spectrum_comp.lib: library file of the compatibility DLL that
		simulates MI cards when findng M2i cards. Please include
		this file instead of spectrum.lib for all projects that
		should use M2i cards with the MI software interface
