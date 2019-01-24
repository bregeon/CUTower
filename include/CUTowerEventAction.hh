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
// $Id: CUTowerEventAction.hh,v 1.1 2007/08/08 10:07:53 bregeon Exp $
// GEANT4 tag $Name:  $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef CUTowerEventAction_h
#define CUTowerEventAction_h 1

#include "CUTowerRootTree.hh"

#include "G4UserEventAction.hh"
#include "globals.hh"
#include "G4ThreeVector.hh"

class CUTowerRunAction;
class CUTowerEventActionMessenger;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class CUTowerEventAction : public G4UserEventAction
{
 public:
   CUTowerEventAction(CUTowerRunAction*);
  ~CUTowerEventAction();

 public:
   void  BeginOfEventAction(const G4Event*);
   void    EndOfEventAction(const G4Event*);
    
   void AddLayer(G4double enedep, G4int ilayer, G4int icylinder);                     
   void SetPrintModulo(G4int    val)  {printModulo = val;};
   
   void SetInitialDir(const G4ThreeVector& aValue) {initialDir = aValue;}; 
   G4ThreeVector GetInitialDir()		   {return initialDir;}; 
   void SetFinalDir(const G4ThreeVector& aValue)   {finalDir = aValue;}; 
   G4ThreeVector GetFinalDir()			   {return finalDir;}; 

    
 private:
   CUTowerRunAction*  runAct;
   
   G4double totalenergy;
   G4double edepLayer[100];
   G4double edepColumn[100];
   G4double edepLog[100][100];

   G4int nHitLayer[36];
  
                  
   G4int     printModulo;

   G4ThreeVector initialDir;
   G4ThreeVector finalDir;
   G4int trackerCollID;                         
                                                    
   CUTowerEventActionMessenger*  eventMessenger;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

    
