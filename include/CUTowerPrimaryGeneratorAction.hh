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
// $Id: CUTowerPrimaryGeneratorAction.hh,v 1.2 2007/08/20 10:19:49 bregeon Exp $
// GEANT4 tag $Name:  $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef CUTowerPrimaryGeneratorAction_h
#define CUTowerPrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "globals.hh"
#include "G4ThreeVector.hh"

class G4ParticleGun;
class G4Event;
class CUTowerDetectorConstruction;
class CUTowerPrimaryGeneratorMessenger;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class CUTowerPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
  public:
    CUTowerPrimaryGeneratorAction(CUTowerDetectorConstruction*);    
   ~CUTowerPrimaryGeneratorAction();

  public:
    void GeneratePrimaries(G4Event*);
    void SetRndmFlag(G4String val) 	   { rndmFlag = val;}
    void SetBeamFlag(G4String val) 	   { beamFlag = val;} 
    void SetGunWidth(G4double val) 	   { gunWidth = val;}
    void SetGunEnergy(G4double val)	   { gunEnergy = val;}
    void SetGunPosition(G4ThreeVector val) { gunPosition = val;}
    void SetGunParticle(G4String val)      { gunParticleName = val;}
    void SetGunAngle(G4double val)         { gunAngle = val;}
    void SetGunDivergence(G4double val)    { gunDiv = val;}
    void SetRunNumber(G4int val)	   { runNo = val;}

  private:
    G4ParticleGun*			particleGun;	  //pointer a to G4  class
    CUTowerDetectorConstruction*	CUTowerDetector;  //pointer to the geometry
    
    CUTowerPrimaryGeneratorMessenger* gunMessenger; //messenger of this class
    G4String		rndmFlag;	  //flag for a rndm impact point
    G4String		beamFlag;	  //flag for reading from beamtest06
    G4ThreeVector	gunPosition;	  //Gun Position
    G4double		gunWidth;	  //Gun Width
    G4double		gunEnergy;	  //Gun Energy
    G4String		gunParticleName; //Gun particle name
    G4double		gunAngle;	  //Gun Width
    G4double		gunDiv;	  //Gun Width
    G4int		runNo;
    G4int		seedNo;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif


