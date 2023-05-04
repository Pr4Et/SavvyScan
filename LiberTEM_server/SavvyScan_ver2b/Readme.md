SavvyScan ver2b  (3 May 2023)

The files included in the subfolder are replacement for version 2.
What is new:
1. Client for LiberTEM server with DECTRIS pixelated detector (ARINA / Simplon): added channel 8 for images from LiberTEM-live processing of 4D-STEM data. Connection via ZMQ TCP/IP sockets.
2. Optional suppression of channels 1-6 (channel 7 is meant for HAADF detector, which can work in parallel to the 4D-STEM pixelated detector).