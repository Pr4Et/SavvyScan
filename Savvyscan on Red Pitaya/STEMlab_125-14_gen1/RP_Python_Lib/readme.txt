To allow connection between Red Pitaya board and Dectris server create (mkdir) the following folder in Red Pitaya Linux 
/opt/velan/Lib
Then copy (scp in the Windows computer) the file Lib.zip and extract (unzip in Red Pitaya Linux) the content to this folder.

To enable easy downloading of hdf5 files from Dectris server copy the velan folder also to the Windows computer under the Start_SerialEM folder, install python3.9, and add a target folder (d:\ArinaData).
Download the h5 files by executing run_Download.bat on a Windows computer.