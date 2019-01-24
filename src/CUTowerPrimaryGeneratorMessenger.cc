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
// $Id: CUTowerPrimaryGeneratorMessenger.cc,v 1.2 2007/08/20 10:19:49 bregeon Exp $
// GEANT4 tag $Name:  $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "CUTowerPrimaryGeneratorMessenger.hh"

#include "CUTowerPrimaryGeneratorAction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWith3VectorAndUnit.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithAnInteger.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

CUTowerPrimaryGeneratorMessenger::CUTowerPrimaryGeneratorMessenger(
                                          CUTowerPrimaryGeneratorAction* CUTowerGun)
:CUTowerAction(CUTowerGun)
{
  gunDir = new G4UIdirectory("/N03/gun/");
  gunDir->SetGuidance("PrimaryGenerator control");
   
  PosCmd = new G4UIcmdWith3VectorAndUnit("/N03/gun/pos",this);
  PosCmd->SetGuidance("Beam Position");
  PosCmd->SetGuidance("  Choice :  X Y Z");
  PosCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  WidthCmd = new G4UIcmdWithADoubleAndUnit("/N03/gun/width",this);
  WidthCmd->SetGuidance("Beam Width (when random on)");
  WidthCmd->SetParameterName("BeamWidth",false);
  WidthCmd->SetRange("BeamWidth>=0 && BeamWidth<1000");
  WidthCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  AngleCmd = new G4UIcmdWithADoubleAndUnit("/N03/gun/angle",this);
  AngleCmd->SetGuidance("Beam Angle with unit");
  AngleCmd->SetParameterName("BeamAngle",false);
  AngleCmd->SetRange("BeamAngle>=0 && BeamAngle<90");
  AngleCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  DivCmd = new G4UIcmdWithADoubleAndUnit("/N03/gun/div",this);
  DivCmd->SetGuidance("Beam Div in mrad");
  DivCmd->SetParameterName("BeamDiv",false);
  DivCmd->SetRange("BeamDiv>=0 && BeamDiv<50");
  DivCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  EnergyCmd = new G4UIcmdWithADoubleAndUnit("/N03/gun/energy",this);
  EnergyCmd->SetGuidance("Beam Energy");
  EnergyCmd->SetParameterName("BeamEnergy",false);
  EnergyCmd->SetRange("BeamEnergy>0");
  EnergyCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  PartCmd = new G4UIcmdWithAString("/N03/gun/particle",this);
  PartCmd->SetGuidance("Choose gun particle.");
  PartCmd->SetParameterName("particle",true);
  PartCmd->SetDefaultValue("mu+");
  PartCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  RndmCmd = new G4UIcmdWithAString("/N03/gun/rndm",this);
  RndmCmd->SetGuidance("Shoot randomly the incident particle.");
  RndmCmd->SetGuidance("  Choice : on(default), off");
  RndmCmd->SetParameterName("choice",true);
  RndmCmd->SetDefaultValue("on");
  RndmCmd->SetCandidates("on off");
  RndmCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  BeamCmd = new G4UIcmdWithAString("/N03/gun/beam",this);
  BeamCmd->SetGuidance("Choose shoot particle from beam06");
  BeamCmd->SetGuidance("  Choice : off(default), on");
  BeamCmd->SetParameterName("choice",true);
  BeamCmd->SetDefaultValue("off");
  BeamCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  RunNumberCmd = new G4UIcmdWithAnInteger("/Cern/random/run", this);
  RunNumberCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

CUTowerPrimaryGeneratorMessenger::~CUTowerPrimaryGeneratorMessenger()
{
  delete RndmCmd;
  delete gunDir;
  delete PosCmd;
  delete WidthCmd;
  delete EnergyCmd;
  delete PartCmd;
  delete BeamCmd;
  delete AngleCmd;
  delete DivCmd;
  delete RunNumberCmd;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void CUTowerPrimaryGeneratorMessenger::SetNewValue(
                                        G4UIcommand* command, G4String newValue)
{ 
  if( command == RndmCmd )
   { CUTowerAction->SetRndmFlag(newValue);}

  if( command == PosCmd )
   { CUTowerAction->SetGunPosition(PosCmd->GetNew3VectorValue(newValue));}

  if( command == WidthCmd )
   { CUTowerAction->SetGunWidth(WidthCmd->GetNewDoubleValue(newValue));}

  if( command == AngleCmd )
   { CUTowerAction->SetGunAngle(AngleCmd->GetNewDoubleValue(newValue));}

  if( command == DivCmd )
   { CUTowerAction->SetGunDivergence(DivCmd->GetNewDoubleValue(newValue));}

  if( command == EnergyCmd )
   { CUTowerAction->SetGunEnergy(EnergyCmd->GetNewDoubleValue(newValue));}

  if( command == PartCmd )
   { CUTowerAction->SetGunParticle(newValue);}

  if( command == BeamCmd )
    { CUTowerAction->SetBeamFlag(newValue);}

 if( command == RunNumberCmd )
   { CUTowerAction->SetRunNumber(RunNumberCmd->GetNewIntValue(newValue));}


}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

