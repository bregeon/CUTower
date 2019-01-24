//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// $Id: CUTowerDetectorMessenger.hh,v 1.1 2007/08/08 10:07:52 bregeon Exp $
// GEANT4 tag $Name:  $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef CUTowerDetectorMessenger_h
#define CUTowerDetectorMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class CUTowerDetectorConstruction;
class G4UIdirectory;
class G4UIcmdWithAString;
class G4UIcmdWithAnInteger;
class G4UIcmdWithADouble;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithoutParameter;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class CUTowerDetectorMessenger: public G4UImessenger
{
  public:
    CUTowerDetectorMessenger(CUTowerDetectorConstruction* );
   ~CUTowerDetectorMessenger();
    
    void SetNewValue(G4UIcommand*, G4String);
    
  private:
    CUTowerDetectorConstruction* CUTowerDetector;
    
    G4UIdirectory*             N03Dir;
    G4UIdirectory*             detDir;
    G4UIcmdWithAnInteger*      NbLayersCmd;    
    G4UIcmdWithAnInteger*      NbColumnsCmd;    
    G4UIcmdWithADouble*        LayerGapCmd;    
    G4UIcmdWithADouble*        ColumnGapCmd;    
    G4UIcmdWithADouble*        AlTkrWidthX0Cmd;    
    G4UIcmdWithAnInteger*      HoneyCombTypeCmd;    
    G4UIcmdWithoutParameter*   UpdateCmd;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

