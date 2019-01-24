# $Id: GNUmakefile,v 1.1 2007/08/08 10:04:26 bregeon Exp $
# --------------------------------------------------------------
# GNUmakefile for examples module.  Gabriele Cosmo, 06/04/98.
# --------------------------------------------------------------

name := CUTower
G4TARGET := $(name)
G4EXLIB := true

#CPPVERBOSE := 1

ifndef G4INSTALL
  G4INSTALL = ../../..
endif

.PHONY: all
all: lib bin


# ROOT stuff
CPPFLAGS += `root-config --cflags`
EXTRALIBS += `root-config --glibs`

include $(G4INSTALL)/config/binmake.gmk

visclean:
	rm -f g4*.prim g4*.eps g4*.wrl
	rm -f .DAWN_*

