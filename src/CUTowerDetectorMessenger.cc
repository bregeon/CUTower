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
// $Id: CUTowerDetectorMessenger.cc,v 1.1 2007/08/08 10:07:55 bregeon Exp $
// GEANT4 tag $Name:  $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "CUTowerDetectorMessenger.hh"

#include "CUTowerDetectorConstruction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithoutParameter.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

CUTowerDetectorMessenger::CUTowerDetectorMessenger(
                                           CUTowerDetectorConstruction* CUTowerDet)
:CUTowerDetector(CUTowerDet)
{ 
  N03Dir = new G4UIdirectory("/N03/");
  N03Dir->SetGuidance("UI commands of this example");
  
  detDir = new G4UIdirectory("/N03/det/");
  detDir->SetGuidance("detector control");
       
   
  NbLayersCmd = new G4UIcmdWithAnInteger("/N03/det/setNbOfLayers",this);
  NbLayersCmd->SetGuidance("Set number of layers.");
  NbLayersCmd->SetParameterName("NbLayers",false);
  NbLayersCmd->SetRange("NbLayers>=0 && NbLayers<100");
  NbLayersCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  NbColumnsCmd = new G4UIcmdWithAnInteger("/N03/det/setNbOfColumns",this);
  NbColumnsCmd->SetGuidance("Set number of Columns.");
  NbColumnsCmd->SetParameterName("NbColumns",false);
  NbColumnsCmd->SetRange("NbColumns>=0 && NbColumns<100");
  NbColumnsCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  LayerGapCmd = new G4UIcmdWithADouble("/N03/det/setLayerGap",this);
  LayerGapCmd->SetGuidance("Set the gap between CAL layer.");
  LayerGapCmd->SetParameterName("LayerGap",false);
  LayerGapCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  ColumnGapCmd = new G4UIcmdWithADouble("/N03/det/setColumnGap",this);
  ColumnGapCmd->SetGuidance("Set the gap between CAL column");
  ColumnGapCmd->SetParameterName("ColumnGap",false);
  ColumnGapCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  AlTkrWidthX0Cmd = new G4UIcmdWithADouble("/N03/det/setAlTkrWidthX0",this);
  AlTkrWidthX0Cmd->SetGuidance("Set a fake monoblock Al Tkr of X0 width.");
  AlTkrWidthX0Cmd->SetParameterName("AlTkrWidthX0",false);
  AlTkrWidthX0Cmd->SetRange("AlTkrWidthX0>=0 && AlTkrWidthX0<10");
  AlTkrWidthX0Cmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  HoneyCombTypeCmd = new G4UIcmdWithAnInteger("/N03/det/setHoneyCombType",this);
  HoneyCombTypeCmd->SetGuidance("Set HoneyComb Type:");
  HoneyCombTypeCmd->SetGuidance("Type 0 - nothing");
  HoneyCombTypeCmd->SetGuidance("Type 1 - homogeneous");
  HoneyCombTypeCmd->SetGuidance("Type 2 - realistic structure (default)");
  HoneyCombTypeCmd->SetParameterName("HoneyCombType",false);
  HoneyCombTypeCmd->SetRange("HoneyCombType>=0 && HoneyCombType<3");
  HoneyCombTypeCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  UpdateCmd = new G4UIcmdWithoutParameter("/N03/det/update",this);
  UpdateCmd->SetGuidance("Update calorimeter geometry.");
  UpdateCmd->SetGuidance("This command MUST be applied before \"beamOn\" ");
  UpdateCmd->SetGuidance("if you changed geometrical value(s).");
  UpdateCmd->AvailableForStates(G4State_Idle);
      
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

CUTowerDetectorMessenger::~CUTowerDetectorMessenger()
{
  delete NbLayersCmd;
  delete NbColumnsCmd;
  delete LayerGapCmd;
  delete ColumnGapCmd;
  delete AlTkrWidthX0Cmd;
  delete HoneyCombTypeCmd;
  delete UpdateCmd;
  delete detDir;
  delete N03Dir;  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void CUTowerDetectorMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{ 
  if( command == NbLayersCmd )
   { CUTowerDetector->SetNbOfLayers(NbLayersCmd->GetNewIntValue(newValue));}
  
  if( command == NbColumnsCmd )
   { CUTowerDetector->SetNbOfColumns(NbColumnsCmd->GetNewIntValue(newValue));}
  
  if( command == LayerGapCmd )
   { CUTowerDetector->SetLayerGap(LayerGapCmd->GetNewDoubleValue(newValue));}

  if( command == ColumnGapCmd )
   { CUTowerDetector->SetColumnGap(ColumnGapCmd->GetNewDoubleValue(newValue));}

  if( command == AlTkrWidthX0Cmd )
   { CUTowerDetector->SetAlTkrWidthX0(AlTkrWidthX0Cmd->GetNewDoubleValue(newValue));}

  if( command == HoneyCombTypeCmd )
   { CUTowerDetector->SetHoneyCombType(HoneyCombTypeCmd->GetNewIntValue(newValue));}
  
  if( command == UpdateCmd )
   { CUTowerDetector->UpdateGeometry(); }

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
