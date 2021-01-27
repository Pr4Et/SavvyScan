This directory contains a copy of the Mini-XML library version 2.10, with
mxml-node.c, which was updated with a fixed version on 4/8/17..
None of the mxml-*.c files should be modified from the distributed form except
for their $Id: lines; doing so would trigger various requirements in the
Library GPL.

The config.h was generated under Linux and modified to work in Windows.

testmxml can be compiled to run with the IMOD library with a command such as:
gcc -o testmxml-IMOD testmxml.c -L $IMOD_DIR/lib -limxml
