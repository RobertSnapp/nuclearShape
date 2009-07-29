# makefile for project liblsm
# Created by Robert R. Snapp
# Copyright (c) 2008

INSTALL_LIBDIR = ~/lib
INSTALL_INCLUDEDIR = ~/include

HEADERS = lsm.h  lsmVoxelImage.h
SOURCE =  lsm.cc lsmVoxelImage.cc
OBJECTS = lsm.o  lsmVoxelImage.o
INCLUDE_LSM = ../include
CFLAGS =  -I$(INCLUDE_LSM) -I/usr/local/include -g -Dmac_osx -fno-common

GFXLIB =  -framework GLUT -framework OpenGL \
          -framework Foundation  -framework vecLib 

LIB = liblsm.1.0.0.dylib
#LIB = liblsm.a

CPP = g++

lsm.o : lsm.cc lsm.h
	$(CPP) $(CFLAGS) -c lsm.cc

lsmVoxelImage.o : lsmVoxelImage.cc lsm.h
	$(CPP) $(CFLAGS) -c lsmVoxelImage.cc

liblsm.1.0.0.dylib : $(OBJECTS)
	$(CPP) -dynamiclib -install_name ~/lib/liblsm.1.dylib \
	-compatibility_version 1.0 -current_version 1.0.0 \
	-undefined suppress -flat_namespace -o liblsm.1.0.0.dylib $(OBJECTS) 

liblsm.a : $(OBJECTS)
	ar r $(LIB) $(OBJECTS)
	ranlib $(LIB)

clean:
	rm $(OBJECTS) $(LIB)

all: $(LIB)

install:
	cp $(LIB) $(INSTALL_LIBDIR)
	cd $(INSTALL_LIBDIR); rm -f liblsm.dylib; ln -s $(LIB) liblsm.dylib;\
		rm -f liblsm.1.dylib; ln -s $(LIB) liblsm.1.dylib
	echo -n "Installed library $(LIB) in directory $(INSTALL_LIBDIR)."
	cp $(HEADERS) $(INSTALL_INCLUDEDIR)
	echo -n "Installed headers $(HEADERS) in directory $(INSTALL_INCLUDEDIR)."
	