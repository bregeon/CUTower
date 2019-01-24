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
// $Id: CUTower.cc,v 1.7 2007/08/20 12:54:14 bregeon Exp $
// GEANT4 tag $Name:  $
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4UIterminal.hh"
#include "G4UItcsh.hh"

#ifdef G4UI_USE_XM
#include "G4UIXm.hh"
#endif

#ifdef G4UI_USE_WIN32
#include "G4UIWin32.hh"
#endif

#include "Randomize.hh"

#ifdef G4VIS_USE
#include "G4VisExecutive.hh"
#endif

#include "PhysicsList.hh"
#include "CUTowerDetectorConstruction.hh"
#include "CUTowerPrimaryGeneratorAction.hh"
#include "CUTowerRunAction.hh"
#include "CUTowerEventAction.hh"
#include "CUTowerSteppingAction.hh"
#include "CUTowerSteppingVerbose.hh"

std::ifstream inFile;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int main(int argc,char** argv)
{
  // Choose the Random engine
  //
  CLHEP::HepRandom::setTheEngine(new CLHEP::HepJamesRandom);

  // User Verbose output class
  //
  G4VSteppingVerbose* verbosity = new CUTowerSteppingVerbose;
  G4VSteppingVerbose::SetInstance(verbosity);
     
  // Construct the default run manager
  //
  G4RunManager * runManager = new G4RunManager;

  // Set mandatory initialization classes
  //
  CUTowerDetectorConstruction* detector = new CUTowerDetectorConstruction;
  runManager->SetUserInitialization(detector);
  //
  G4VUserPhysicsList* physics = new PhysicsList;
  runManager->SetUserInitialization(physics);

  G4UIsession* session=0;
  
  if (argc==1)   // Define UI session for interactive mode.
    {
      // G4UIterminal is a (dumb) terminal
      //
#if defined(G4UI_USE_XM)
      session = new G4UIXm(argc,argv);
#elif defined(G4UI_USE_WIN32)
      session = new G4UIWin32();
#elif defined(G4UI_USE_TCSH)
      session = new G4UIterminal(new G4UItcsh);      
#else
      session = new G4UIterminal();
#endif
    }
  
#ifdef G4VIS_USE
  // Visualization manager
  //
  G4VisManager* visManager = new G4VisExecutive;
  visManager->Initialize();
#endif
    
  // Set user action classes
  //
  G4VUserPrimaryGeneratorAction* gen_action = new CUTowerPrimaryGeneratorAction(detector);
  runManager->SetUserAction(gen_action);
  //
  CUTowerRunAction* run_action = new CUTowerRunAction;  
  runManager->SetUserAction(run_action);
  //
  CUTowerEventAction* event_action = new CUTowerEventAction(run_action);
  runManager->SetUserAction(event_action);
  //
  G4UserSteppingAction* stepping_action =
    new CUTowerSteppingAction(detector, event_action);
  runManager->SetUserAction(stepping_action);
  
  // Initialize G4 kernel
  //
  //runManager->Initialize();
    
  // Get the pointer to the User Interface manager
  //
  G4UImanager* UI = G4UImanager::GetUIpointer();  

  if (session)   // Define UI session for interactive mode
    {
      // G4UIterminal is a (dumb) terminal
      //
      UI->ApplyCommand("/control/execute vis.mac");    
#if defined(G4UI_USE_XM) || defined(G4UI_USE_WIN32)
      // Customize the G4UIXm,Win32 menubar with a macro file
      //
      UI->ApplyCommand("/control/execute visTutor/gui.mac");
#endif
      session->SessionStart();
      delete session;
    }
  else           // Batch mode
    { 
#ifdef G4VIS_USE
      visManager->SetVerboseLevel("quiet");
#endif
      G4String command = "/control/execute ";
      G4String fileName = argv[1];
      UI->ApplyCommand(command+fileName);
    }

  // Job termination
  // Free the store: user actions, physics_list and detector_description are
  //                 owned and deleted by the run manager, so they should not
  //                 be deleted in the main() program !

#ifdef G4VIS_USE
  delete visManager;
#endif
  delete runManager;
  delete verbosity;

  return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
