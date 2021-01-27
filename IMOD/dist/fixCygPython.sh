#!/bin/bash
# This script will make sure that there is a current python.exe in cygwin
# Newer cygwin has only a cygwin link from python to python2.x.exe
# Due to some arcane memory operations related to security vulnerabilities,
# the file now needs ot be a hard link instead of a copy
# This will not work when called from non-cygwin applications (java)
PATH="/bin:$PATH"
makeLink=0
if [[ -e /bin/python.exe && -h /bin/python ]] ; then

    # The IMOD installStub has been copying or linking the link target to python.exe
    # If it already exists, check whether that copy is the same as the current link target
    diff -q /bin/python.exe /bin/python > /dev/null 2>&1
    if [ $? -ne 0 ] ; then 
        echo "The file /bin/python.exe appears to be out of date"
        echo "Making a new hard link of the current python to /bin/python.exe"
        makeLink=1

        # Then check whether it is still a copy instead of a hard link
    elif [[ "$(stat -c %h python.exe)" -eq 1 ]] ; then
        echo "The file /bin/python.exe is a copy, not a hard link"
        echo "Making a hard link instead, of the current python to /bin/python.exe"
        makeLink=1
    fi
    
elif [[ ! -e /bin/python.exe ]] ; then

    # Or, if there is no python.exe and there is a link, make it a hard link too
    if [ -h /bin/python ] ; then 
        echo "The Cygwin link to python will not work for IMOD"
        echo "Making a hard link of the current python to /bin/python.exe"
        makeLink=1
    else
        echo "There is no Python installed in Cygwin"
        exit 1
    fi
fi
if [ $makeLink -eq 1 ] ; then
    ln -Lf /bin/python /bin/python.exe
    if [ $? -ne 0 ] ; then 
        echo "There was an error copying to /bin/python.exe"
        exit 1
    fi
fi
exit 0
