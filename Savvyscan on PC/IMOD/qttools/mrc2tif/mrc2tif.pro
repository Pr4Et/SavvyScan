SOURCES += mrc2tif.cpp iiqimage.cpp tiff.c
TEMPLATE = app
CONFIG += qt
INCLUDEPATH += . ../../include

include (qconfigure)

INSTALLS += target

tiffc.target = 
tiffc.depends = ../../mrc/tiff.c
tiffc.commands = cp -f ../../mrc/tiff.c .

iiqim.target = 
iiqim.depends = ../../3dmod/iiqimage.cpp
iiqim.commands = cp -f ../../3dmod/iiqimage.cpp .

QMAKE_EXTRA_TARGETS += tiffc iiqim
