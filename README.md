# SavvyScan
## A flexibility upgrade for multi-channel scanning transmission electron microscopy

Savvyscan software source files. Jan 2021.
Editted/ written by Shahar Seifer, Elbaum lab, Weizmann Institute of Science, Israel.

Details in: Shahar Seifer, Lothar Houben, Michael Elbaum, "Flexible STEM with Simultaneous Phase and Depth Contrast", Microscopy and Microanalysis (2021), 27, 1476–1487.

Compilation and installation:
The main project file is \ServerCamera\SERVER-SEMCamServer.vcxproj.
Project can be built in VisualStudio community 2019 with full V142 platform toolset (including MFC), in 64 bits environment.
Edit files \ServerCamera\BaseServer.cpp and \Spectrum\shahar\RecRepMulti.cpp to change the directory "d:\SavvyscanData" into which is your output data directory.
Install SeiralEM interface using the settings suggested in directory Install_in_SerialEM.
The system runs by starting the camera server, then the shadow GUI and SerialEM.
Output file types: *.MRC and *.MAT. 
MAT files can be displayed with Matlab programs found in directory Post_Processing. 
