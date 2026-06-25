## How to run the LiberTEM server: <br>
Start the DECTRIS server and pixelated detector. Type in a terminal of the LiberTEM server computer: <br>
conda activate libertem  <br>
./LiberTEM_passive.py  <br>

## Setting SavvyScan ver2b to transfer images from LiberTEM server to SerialEM: <br>
Choose in the LiberTEM dropdown menu: CONNECT, click SEND. Then choose one of the process types, edit the parameters, and click SEND. <br>
Select "Channel sent to SerialEM": 8   <br>
Click "Arm 4D STEM" before the scan (optionally: set multislice for a sequence of record scans).<br>
Only scans with duration larger than the chosen "threshold" time will activate the pixelated detector. <br>


## How to install on LiberTEM server on a Linux computer (with anaconda pre-installed): <br>
conda activate <br>
conda install gh --channel conda-forge <br>
gh auth login <br>
conda create -n libertem python=3.10 <br>
conda activate libertem <br>
gh repo clone LiberTEM/LiberTEM <br>
gh repo clone LiberTEM/LiberTEM-live <br>
cd LiberTEM-live <br>
gh pr checkout 74 <br>
cd .. <br>
pip install -e ./LiberTEM   <br>
pip install -e ./LiberTEM-live   <br>

The file LiberTEM_passive.py should be placed in the same directory, then using text editor: <br>
modify the ip address of Dectris server (instead of "192.168.100.70") <br>
modify the folder name /home/stem to one on your system. <br>
If there is GPU and CUDA installed then you can increase the size limit for SSB processing in the line: <br>
elif ptype=="SSB" and nx>128:  <br>