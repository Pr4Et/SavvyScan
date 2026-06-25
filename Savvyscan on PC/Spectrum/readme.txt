******************************************************************************

spcm_drv_c - c/c++ examples for the Spectrum SpcM driver interface.

******************************************************************************

(c) Spectrum GmbH
    Ahrensfelder Weg 13-17
	22927 Grosshansdorf
	Germany

    www.spectrum-instrumentation.com

******************************************************************************

This readme file should give you an overview of the different examples and
the content of the different directories.

Documentation for the API as well as a detailed description of the hardware
can be found in the manual for each device which can be found on our website:
www.spectrum-instrumentation.com/en/downloads

Further information can be found online in the Knowledge Base:
www.spectrum-instrumentation.com/en/knowledge-base-overview

The SpcMDrv C/C++ examples from Spectrum are only using general c commands
and have been tested with Microsoft Visual C++, Borland C++ Builder under
Windows and and gnu c++ under Linux. Each sub directory contains the source
code of one example, Visual C++ and Borland C++ Builder project files and
one or more makefiles for gnu c++.

******************************************************************************

directory c_header
-----------------
contains all the standard header and library files that are needed to
access the driver. Please do not modify any of these files as they
provide the link to the driver itself. Keep this directory updated
with every driver update you get to have access to new features.

file dlltyp.h:
	defines standard types for different platforms and compiler. All our
	examples only use these standard types to be sure to have a defined type
	declaration. If using other C/C++ compilers than we you should make a copy
	of this file and setup the type declatation matching the new C/C++
	compiler.

file regs.h:
	defines all software registers. This file contains all software registers
	for a couple of different cards and also for older driver versions. Please
	stay to the manual to see which registers can be used with your card.

file spcerr.h:
	contains the error codes that the driver can return. Again inhere are also
	error codes that are only used for other card types or for older driver
	versions.

file spcm_drv.h:
	defines the driver interface of the SpcMDrv

file spcm_win32_msvcpp.lib
	Microsoft Visual C++ library file for 32 bit windows platforms

file spcm_win32_bcppb.lib
	Borland C++ Builder library file for 32 bit windows platforms

All other files of this directory are not used with the SpcMDrv

******************************************************************************



directory common
----------------
This directory contains some c-files and header files that are used by all
examples. These common source files mainly group driver acces functions to
make the examples more simple to understand. Feel free to use the common
functions for your own programs. All these functions are distributed "as is".

******************************************************************************



directory rec_std_single
------------------------
Example for acquistion cards using the standard singleshot mode.

******************************************************************************



directory rec_std_multi
-----------------------
Example for acquistion cards using the standard Multiple Recording mode. If
timestamp is installed on the card the timestamps of each segment are also
read and displayed.

******************************************************************************



directory rec_std_gate
----------------------
Example for acquisition cards using the standard Gated Sampling mode together
with timestamps. Position and length of each gate segment is calculated and
display and access to the gate segment is shown.

******************************************************************************



directory rec_fifo_single
-------------------------
Examples for acquistion cards using the FIFO mode with continuous acquisition
of data. Control of the data transfer is located in a thread. The examples show
the buffer handling, data separation, calculation and writing of data to file

******************************************************************************



directory rec_fifo_multi
------------------------
Examples for acquisition cards using the FIFO mode with continuous acquisition
of data and Multiple Recording option. If timestamp is also installed
timestamps are recorded and can be also written to file. The examples shows
the buffer handling (including timestamps buffer), data separation, calculation
and writing of data to file

******************************************************************************



directory rec_fifo_aba
------------------------
Examples for acquisition cards using the FIFO mode with continuous acquisition
of data, Multiple Recording option, ABA option and timestamp option. The
examples shows the buffer handling with three totally independent working
threads each responsible for one type of transfer. A fourth main thread
is collecting results from the working threads and displays average and
timestamp distance

******************************************************************************



directory rep_std_single
------------------------
Example for single replay modes. This example covers analog and digital output
cards and allows to test single shots, continous output and single restart.
The signal shape can be selected.

******************************************************************************



directory rep_sequence
------------------------
Example for sequence replay mode. This example covers analog and digital output
cards and shows sequence replay mode as simple sequence and with sequence change
at runtime.


******************************************************************************


directory rec_std_single_sync
---------------------
Contains simple examples showing the synchronization of multiple cards with
the star-hub. The examples only show the synchronization part without any
further data processing.
It also contains examples for usage of system star-hub.

******************************************************************************



directory dll_loading
---------------------
The example shows how one can access the dll functions from a non supported
C/C++ compiler. The example shows how to load a Windows library and how to
access DLL functions without using a lib-file. However it is recommended to
use a lib conversion tool that is normally included into the compiler to
directly access the dll.



directory netbox_embedded_server
--------------------------------
This directory contains examples that show a client-server communication where 
the netbox acts as server and the user's PC as client.
Another example show how to send emails with acquired data using libcurl.



directory sse
--------------------------------
This directory contains some helpful functions that use SSE commands for
higher speed.



directory test
--------------------------------
This directory contains several small tools for testing the cards.
