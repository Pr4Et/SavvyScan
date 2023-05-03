How to run the LiberTEM server:
conda activate libertem
./LiberTEM_passive.py

How to set SavvyScan ver2b to transfer images from LiberTEM server to SerialEM:
Click on LiberTEM dropdown menu: CONNECT, click <SEND>, then choose one of the process types and click <SEND>.
Select <Channel sent to SerialEM>: 8

How to install on Linux (assuming anaconda pre-installed):
conda activate
conda install gh --channel conda-forge
gh auth login
conda create -n libertem python=3.10
conda activate libertem
gh repo clone LiberTEM/LiberTEM
gh repo clone LiberTEM/LiberTEM-live
cd LiberTEM-live
gh pr checkout 74
cd ..
pip install -e ./LiberTEM  
pip install -e ./LiberTEM-live  

The file LiberTEM_passive.py should be placed in the same directory, then using text editor: 
modify the ip address of Dectris server (instead of "192.168.100.70")
modify the folder name /home/stem to one on your system.
If there is GPU and CUDA installed then you can increase the size limit for SSB processing in the line:
elif ptype=="SSB" and nx>128: 