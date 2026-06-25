#!/bin/sh
#
# Stub for self-extracting IMOD installer
# Use -help to see the options
#

installer=installIMOD
this="$0"
rootname=${this##*/}
packname=${rootname%.*}.tar.gz
tempdir=IMODtempDir
installdir=""
namedir=""
nameopt=""
custom=0
skipopt=""
scriptopt=""
sysdir=""
yesopt=""
debopt=""
ubuntu=0
ubscripts1=/etc/csh.cshrc
ubscripts2=/etc/bash.bashrc

# Find the end of this script
#
stubline=`awk '/END OF STUB/ {if (gotit > 0) {print NR + 1 ; exit} else { gotit = 1}}' "$this"`
if [ -z $stubline ]; then
    echo "Could not find package in this file!"
    exit 1
fi

extract=0
test=0
while [ $# -gt 0 ]
do
    key="$1"

    case $key in
        -ex*)
            extract=1
            ;;

        -test)
            test=1
            ;;

        -di*)
            installdir=$2
            shift
            ;;

        -na*)
            namedir=$2
            nameopt="-name $namedir"
            spacetst=`echo $namedir | grep '[ /\\]'`
            if [ -n "$spacetst" ]; then
                echo "The entry for -name must be a single directory name"
                echo "with no spaces or directory separators"
                exit 1
            fi
            shift
            ;;

        -sc*)
            sysdir=$2
            shift
            ;;
        
        -sk*)
            skipopt="-skip"
            ;;

        -y*)
            yesopt="-yes"
            ;;
         
        -de*)
            debopt="-debian"
            ubuntu=1
            ;;

        -h*|--help)
            cat <<EOF
This is a self-installing IMOD package which will install IMOD in a default
or specified location and set up startup files.
Options can be added after the name of the package on the command line.  
Options may be abbreviated.  The options are:
  -dir dir     Install IMOD within the given directory.  The directory must
                 exist (unless -name is entered also) and you must have
                 permission to write to it.  
  -name dir    Rename the IMOD directory to the given name.  The script will
                 automatically remove an existing version by this name and not
                 offer to clean up of other versions
  -extract     Just extract the package and the real IMOD install script into
                 a directory named $tempdir, without installing
  -script dir  Place the IMOD startup scripts in the given directory.  On Linux
                 and Windows, scripts will be copied there instead of to 
                 /etc/profile.d.  On Mac, scripts will be copied there
                 and the system startup files will not be modified.
  -skip        Skip copying startup scripts (Linux and Windows) or modifying
                  system startup files (Mac)
  -debian      Modify $ubscripts1 and $ubscripts2 instead of 
                  copying startup scripts to /etc/profile.d
                   - this option is set automatically for Ubuntu versions < 10
  -yes         Install this version and remove old versions without asking for
                  confirmation
  -test        Run this self-installing script but not the real install script
  -h, --help   Print this help
EOF
            exit 0
            ;;

        *)
            echo "Illegal option $1"
            exit 1
            ;;
    esac
    shift
done    

sysfiles1=
sysfiles2=    
sysfiles3=

system=`uname -s`
defaultdir=/usr/local

case $system in
    *Linux*)
         os=linux
         scripts1=IMOD-linux.csh
         scripts2=IMOD-linux.sh
         if [ $ubuntu -eq 0 ]; then
             relfiles=`\find /etc -maxdepth 1 -type f -name '*release' -print`
             if [ -n "$relfiles" ]; then
                 ubtest=`grep -i ubuntu $relfiles`
                 if [ -n "$ubtest" ]; then 
                     if [ -e /etc/lsb-release ]; then
                         major=`sed -n '/RELEASE/s/.*=\([0-9]*\).*/\1/p' /etc/lsb-release`
                         if [ $major -lt 10 ]; then ubuntu=1 ; fi
                     else
                         ubuntu=1
                     fi
                 fi
             fi
         fi
         if [ $ubuntu -eq 1 ]; then
             sysfiles2=$ubscripts2
         fi
         ;;

    *CYGWIN*)
        os=windows
        scripts1=IMOD-cygwin.csh
        scripts2=IMOD-cygwin.sh
        ;;

    *Darwin*)
        os=osx
        defaultdir=/Applications
        scripts1=mac.cshrc
        scripts2=mac.profile
        sysfiles1=/etc/csh.login
        sysfiles2=/etc/profile
        if [ -e /etc/zprofile ] ; then sysfiles3=/etc/zprofile ; fi
        ;;

    *)
        if [ $extract -eq 0 ]; then
            echo "IMOD will not run on this system ($system)"
            echo "You can add the option -extract to just extract the tar file"
            exit 1
        fi
        ;;    

esac

# For install and script dir, need to make an absolute path before
# going into temp dir and running install script
#
if [ -n "$installdir" ]; then
    if [ $os = windows ]; then installdir=`cygpath "$installdir"` ; fi
    absolute=`echo $installdir | grep '^/'`
    if [ $? -ne 0 ]; then installdir="`pwd`""/$installdir" ; fi
fi

copyto=/etc/profile.d
if [ -n "$sysdir" ]; then
    if [ $os = windows ]; then sysdir=`cygpath "$sysdir"` ; fi
    absolute=`echo $sysdir | grep '^/'`
    if [ $? -ne 0 ]; then sysdir="`pwd`""/$sysdir" ; fi
    scriptopt="-script $sysdir"
    copyto="$sysdir"
fi


if [ -z "$installdir" ]; then installdir=$defaultdir ; fi
if [ "$installdir" != "$defaultdir" ]; then custom=1 ; fi

if [ $extract -eq 0 ]; then
    echo ""
    if [ -z "$namedir" ]; then
        echo "This script will install IMOD in $installdir and rename"
        echo "any previous version, or remove another copy of this version."
    else
        echo "This script will install IMOD in $installdir and name it $namedir."
        echo "It will first remove $namedir if it exists, and remove anything that"
        echo "matches the imod_n.n.n name that this IMOD will unpack as"
    fi
    echo ""
    if [ -z "$skipopt" ]; then
        if [ -z "$sysfiles1" ]; then
            echo "It will copy $scripts1 and $scripts2 to $copyto"
        else
            echo -n "It will edit $sysfiles1" 
            if [ -n "$sysfiles3" ] ; then echo -n ", ${sysfiles3}," ; fi
            echo " and $sysfiles2 to source a startup script"
            echo " in $installdir/IMOD unless they already appear to do so"
        fi
    fi
    
    if [ -z "$yesopt" ]; then
        echo ""
        echo "You can add the option -h to see a full list of options"
        echo ""
        /bin/echo -n "Enter Y if you want to proceed: "
        read yesno
        if [ "$yesno" != "Y" ]; then
            if [ "$yesno" != "y" ]; then exit 0 ; fi
        fi
    fi
fi

mkdir -p $tempdir    

echo "Extracting $packname ..."
tail -n +$stubline "$this" > $tempdir/$packname

version=`echo $rootname | sed '/\(imod_[0-9.]*\).*/s//\1/' | sed '/\.$/s///'`

cd $tempdir

echo "Extracting $installer"

tar -xzf $packname $version/$installer

if [ ! -e $version/$installer ]; then
    echo "Could not extract $installer"
    exit 1
fi

mv -f $version/$installer .
rmdir $version

if [ $extract -eq 1 ]; then
    echo "Left $packname and $installer in $tempdir"
    exit 0
fi

# Fix python link in recent Cygwin distributions
python=python
if [ $test -eq 0 ]; then
    if [ $os = windows ]; then
        if [ -e /bin/python.exe ]; then
            if [ -L /bin/python ]; then
                diff -q /bin/python.exe /bin/python > /dev/null 2>&1
                if [ $? -ne 0 ]; then
                    echo " "
                    echo "The copy of python in /bin/python.exe appears to be out of date"
                    echo "Making a new copy of the current python to /bin/python.exe"
                    cp -Lf /bin/python /bin/python.exe
                fi
            fi
        else
            echo " "
            if [ -L /bin/python ]; then
                echo "The Cygwin link to python will not work for IMOD"
                echo "Making a copy of its target that is an actual /bin/python.exe"
                cp -Lf /bin/python /bin/python.exe
            else
                echo "You must have Python installed in Cygwin to run this installer"
                test=1
            fi
        fi
    else
        which python > /dev/null 2>&1
        if [ $? -ne 0 ]; then
            python=
            if [ -e /usr/bin/python3 ] ; then
                echo " "
                echo "There is no /usr/bin/python, just /usr/bin/python3."
                echo "You will need a /usr/bin/python to run scripts in IMOD."
                echo "You can either make a link in /usr/bin to python3 and use that in IMOD,"
                echo "or use python3 just for this install and install python 2 later."
                echo " "
                echo -n "Enter Y to make a link and use python3 in IMOD: "
                read yesno
                if [ "$yesno" != "Y" ]; then
                    if [ "$yesno" != "y" ]; then python=python3 ; fi
                fi
                if [ -z $python ] ; then
                    ln -s /usr/bin/python3 /usr/bin/python
                    if [ $? -ne 0 ]; then
                        python=python3
                        echo "You need to run this installer as root to make the link"
                    else
                        python=python
                    fi
                fi
            else
                echo " "
                echo "There is no Python on the search path that runs with the command \"python\"."
                echo "IMOD cannot be installed without Python."
            fi
        fi
    fi
fi
if [ $os = windows ]; then export PATH="/bin:$PATH" ; fi

retval=0
if [ $test -eq 0 ]; then
    if [ -n $python ] ; then
        $python -u $installer -dir "$installdir" $nameopt $scriptopt $skipopt $debopt $yesopt $packname
        retval=$?
    fi
fi

echo " "
echo "Cleaning up $packname, $installer, and $tempdir"
cd ..
if [ -n "$yesopt" ]; then sleep 2 ; fi
rm -f $tempdir/$installer $tempdir/$packname
rm -rf $tempdir
if [ $? -ne 0 ]; then echo "Sorry, you will have to do that yourself: rm -rf $tempdir" ; fi

exit $retval
    
#  $Id$
#

END OF STUB
