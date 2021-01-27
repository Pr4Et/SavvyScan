This directory contains a library based on cttfind4 from Rohou and
Grigorieff.  It was originally derived from verson 4.1.9, as contained in the
cistem-1.0.0 package.  For most of the components, the complete file was used
from that package, with a few small modifications for IMOD.  
image_trim.cpp was derived from image.cpp by extracting the sections marked as
needed for CTFFIND.  Some additional functions unneeded for the ctffind
function abstracted here were also commented out, to avoid obvious additional
dependencies.  All modifications of these files are marked with "IMOD"
comments.  The more extensive changes in ctffind.cpp are also marked.

This software is covered by the license in IMOD/dist/Janelia.txt
