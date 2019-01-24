#!/bin/tcsh
source $HOME/glast/geant4/v9.3/geant4v9.3.cshrc
setenv G4OUTFILEDIR $G4WORKDIR/results

setenv G4CURRENTMACRO $1
setenv G4OUTFILENAME $G4OUTFILEDIR/$G4CURRENTMACRO
$G4WORKDIR/bin/Linux-g++/CUTower ${G4WORKDIR}/CUTower/mac/$G4CURRENTMACRO.mac | tee /tmp/$G4CURRENTMACRO.txt

cat /tmp/$G4CURRENTMACRO.txt | mail -s 'Job Done' johan.bregeon@pi.infn.it
rm -f /tmp/$G4CURRENTMACRO.txt
