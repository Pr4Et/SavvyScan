#!/bin/sh

LIST_OF_FAILED=()

# ***** compile + error check *****
compile_and_check() # (path, makefile-name)
    {
    echo
    echo
    echo "***** Compiling `pwd $1`/"$1
    make -f $1 clean
    if [ $? != 0 ]; then
        LIST_OF_FAILED[${#LIST_OF_FAILED[*]}]="$(pwd $1)/$1"
        return 1
    fi
    make -f $1
    if [ $? != 0 ]; then
        LIST_OF_FAILED[${#LIST_OF_FAILED[*]}]="$(pwd $1)/$1"
        return 1
    fi
    make -f $1 clean
    return $?
    }

# check C++ examples
for i in `find . -type d -not -iwholename '*.svn*' -not -iwholename '*cuda*'`; do
    pushd $i
    for j in `ls makefile*`; do
        compile_and_check $j
    done
    popd
done

echo
echo ********************************************************************************
echo Failed Builds:

echo $LIST_OF_FAILED
