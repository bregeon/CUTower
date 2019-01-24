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
// $Id: CUTowerSteppingAction.cc,v 1.1 2007/08/08 10:07:55 bregeon Exp $
// GEANT4 tag $Name:  $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "CUTowerSteppingAction.hh"

#include "CUTowerDetectorConstruction.hh"
#include "CUTowerEventAction.hh"

#include "G4Step.hh"
#include "G4VProcess.hh"

////#include "G4RunManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

CUTowerSteppingAction::CUTowerSteppingAction(CUTowerDetectorConstruction* det,
                                         CUTowerEventAction* evt)
:detector(det), eventaction(evt)					 
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

CUTowerSteppingAction::~CUTowerSteppingAction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void CUTowerSteppingAction::UserSteppingAction(const G4Step* aStep)
{
  // get volume of the current step
  G4VPhysicalVolume* volume 
  = aStep->GetPreStepPoint()->GetTouchableHandle()->GetVolume();

//  G4String KnownProcesses[19] = {"Transportation", "LowEnRayleigh", "LowEnPhotoElec", "LowEnergyIoni",
//  				 "LowEnBrem", "msc", "LowEnCompton", "LowEnConversion", "UserMaxStep",
//				 "eBrem", "eIoni", "annihil", "PhotonInelastic", "ElectroNuclear",
//				 "Decay" ,"PositronNuclear", "hElastic", "hInelastic", "nCapture"};
//  G4String CurrentProcess = aStep->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName();
//  bool NewProc = true;
//  for(int i=0; i<19; i++)
//    {
//    if(CurrentProcess == KnownProcesses[i])
//       NewProc = false;
//    }
//  if(NewProc)G4cout <<"Process : "<<CurrentProcess<< G4endl; 

  
  // collect energy and track length step by step
  G4double edep = aStep->GetTotalEnergyDeposit();
    
  //Get energy deposit in piece of cylinder
  G4int ilay = -1, icol = -1, copyNo = -1;
  if(volume->GetLogicalVolume()->GetMaterial() ==  detector->GetCalMaterial())
    {
    copyNo = volume->GetCopyNo();// CopyNo = 100*iLayer + j*iCol
    ilay = (int)floor(copyNo/100.);
    icol = copyNo%100;
    //G4cout<<copyNo<<" "<<ilay<<" "<<icyl<<G4endl; 
    eventaction->AddLayer(edep,ilay,icol);
    }
      

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

