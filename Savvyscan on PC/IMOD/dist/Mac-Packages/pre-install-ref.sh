#!/bin/bash

#-------------------------------------------------------------------------------
#Install IMOD

sh imod_4.10.32_osx64_10.8_CUDA4.1.sh -yes > $HOME/IMOD-install.log
if [ $? -ne 0 ] ; then exit 1 ; fi
rm -f $HOME/IMOD-install.log
exit 0
