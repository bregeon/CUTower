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
// $Id: CUTowerRunAction.cc,v 1.2 2007/08/20 10:19:49 bregeon Exp $
// GEANT4 tag $Name:  $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "CUTowerRunAction.hh"

#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4UnitsTable.hh"

extern std::ifstream inFile;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

CUTowerRunAction::CUTowerRunAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

CUTowerRunAction::~CUTowerRunAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void CUTowerRunAction::BeginOfRunAction(const G4Run* aRun)
{ 
  G4cout << "### Run " << aRun->GetRunID() << " start." << G4endl;
  
  char name[15];
  sprintf(name,"sps_part.dat");
  inFile.open(name);

  //inform the runManager to save random number seed
  G4RunManager::GetRunManager()->SetRandomNumberStore(true);
    
  //create the Root Tree to save data
  theRootTree = new CUTowerRootTree();
  //theRootTree->SetFileNameByTime();
  theRootTree->SetFileNameByEnv();
  theRootTree->OpenFile();
  theRootTree->CreateTree();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void CUTowerRunAction::EndOfRunAction(const G4Run* aRun)
{
  G4int NbOfEvents = aRun->GetNumberOfEvent();
  if (NbOfEvents == 0) return;
   
  //Write file to disk
  theRootTree->Save();
  
  //print
  //
  G4cout<< "\n--------------------End of Run------------------------------\n"
        << G4endl;
     
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
