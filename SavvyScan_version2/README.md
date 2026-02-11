# SavvyScan
## A flexibility upgrade for multi-channel scanning transmission electron microscopy

Savvyscan software source files. Jan 2021.
Editted/ written by Shahar Seifer, Elbaum lab, Weizmann Institute of Science, Israel.

Details in: Shahar Seifer, Lothar Houben, Michael Elbaum, "Flexible STEM with Simultaneous Phase and Depth Contrast", Microscopy and Microanalysis (2021), 27, 1476â€“1487.

Compilation and installation:
The main project file is \ServerCamera\SERVER-SEMCamServer.vcxproj.
Project can be built in VisualStudio community 2019 with full V142 platform toolset (including MFC), in 64 bits environment.
The plugin requires also installtion Microsoft Visual C++ Redistributable (v14, x64) for VS 2015–2022.
Edit files \ServerCamera\BaseServer.cpp and \Spectrum\shahar\RecRepMulti.cpp to change the directory "d:\SavvyscanData" into which is your output data directory.
Install SeiralEM interface using the settings suggested in directory Install_in_SerialEM.
The system runs by starting the camera server, then the shadow GUI and SerialEM.
Output file types: *.MRC and *.MAT. 
MAT files can be displayed with Matlab programs found in directory Post_Processing. 

Version 2: Copy and overwrite the appropriate folders with the incremental project found in folder SavvyScan_version2.

## What's new in SavvyScan version 2 (Feb 2023)
The software has been extensively upgraded based on hundreds of hours of practice in cryo-STEM tasks.
* Supports OPAL detector as well as Dectris Arina ultrafast camera for 4DSTEM.
* Spectrum cards run in FIFO mode, with latest firmware as for end of 2022.
* Spectrum cards synchronization is chosen either via STAR-HUB or cables.
* STEM camera extended to 8K x 8k pixels.
* Acquisition at highest sampling rate and digital averaging to improve signal sampling over integration time.
* Added scan patterns and improvement in existing ones. The raster scan is recommended for general purpose.
* Use of SerialEM python socket to temporary take control over serialEM and the microscope for special tasks, such as automated diffraction disc alignment of a quadrature detector.
* Full support of SeiralEM tilt series acquisition function. Scans longer than threshold time set in the GUI are considered record scans and stored in final MRC stack.
* Many bug fixes and robust performance.

