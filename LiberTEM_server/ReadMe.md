How to run the LiberTEM server: <br>
conda activate libertem  <br>
./LiberTEM_passive.py  <br>

How to set SavvyScan ver2b to transfer images from LiberTEM server to SerialEM: <br>
Click on LiberTEM dropdown menu: CONNECT, click SEND, then choose one of the process types and click SEND. <br>
Select "Channel sent to SerialEM": 8   <br>

How to install on Linux (assuming anaconda pre-installed): <br>
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