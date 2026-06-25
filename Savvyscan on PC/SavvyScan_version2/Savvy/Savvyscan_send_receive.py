'''
Uses Python server for SerialEM and communicates with shadow1 GUI as server (replier)

'''
import numpy as np
import mrcfile
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy import fftpack
import math
import sys
import configparser
import ast
import argparse
import zmq;
import time;

sys.path.insert(0,'C:\\Program Files\\SerialEM\\PythonModules')
import serialem as sem
connect_single=True

context=zmq.Context()
socket=context.socket(zmq.REP) #replier
socket.bind("tcp://*:5580")

reply_txt="00"


while True:
        # wait for request from Jeolcam server via zeromq
        message=socket.recv()
        #message_str=message.decode("utf-8")
        print("Received request: %s" % message)
        if message==b'data':
            socket.send_string(reply_txt)
        elif message==b'continous':
            connect_single=False
            try:
                sem.ConnectToSEM(48888,'192.168.100.20')
            except:
                connect_single=True
            time.sleep(0.5)
            socket.send_string("OK")
        elif message==b'single':
            connect_single=True
            try:
                sem.Exit(1)
            except:
                print("Already in single mode")
            socket.send_string("OK")
        else:
            try:
                if connect_single:
                    sem.ConnectToSEM(48888,'192.168.100.20')
                    time.sleep(0.5)
                time.sleep(0.3)
                #connect_once=False
                x=0
                y=0
                exec(message)
                reply_txt="x={0},y={1}".format(x,y)
                if connect_single:
                    sem.Exit(1)
            except:
                print("Script command failed")
                sem.Exit(1)
                message=b' '
            socket.send_string("OK")
#(xnow,ynow)=sem.ReportDiffractionShift()
#print('xnow, ynow ',xnow,ynow)

        
