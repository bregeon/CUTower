//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: TrackerSD.cc,v 1.1 2007/08/08 10:07:55 bregeon Exp $
// GEANT4 tag $Name:  $
// ------------------------------------------------------------
//      GEANT 4 class implementation file
//      CERN Geneva Switzerland
//
//
//      ------------ TrackerSD  ------
//           by  R.Giannitrapani, F.Longo & G.Santin (13 nov 2000)
//
// ************************************************************

#include "TrackerSD.hh"
#include "TrackerHit.hh"

#include "G4VPhysicalVolume.hh"
#include "G4Step.hh"
#include "G4VTouchable.hh"
#include "G4TouchableHistory.hh"
#include "G4SDManager.hh"

#include "G4ios.hh"

#define NSTRIPS 1500+10000*36  // correspond to 1500 strips/plane * 36 planes

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

TrackerSD::TrackerSD(G4String name)
  :G4VSensitiveDetector(name)
{
  collectionName.insert("TrackerCollection");
  HitID = new G4int[NSTRIPS];
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

TrackerSD::~TrackerSD()
{
  delete [] HitID;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void TrackerSD::Initialize(G4HCofThisEvent*HCE)
{
  TrackerCollection = new TrackerHitsCollection
    (SensitiveDetectorName,collectionName[0]);
  for (G4int j=0;j<NSTRIPS;j++) {
    HitID[j] = -1;   
  };
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4bool TrackerSD::ProcessHits(G4Step* aStep,G4TouchableHistory* ROhist)
{ 
  G4double edep = aStep->GetTotalEnergyDeposit();
  if ((edep/keV == 0.)) return false;      

  G4TouchableHistory* theTouchable
    = (G4TouchableHistory*)(aStep->GetPreStepPoint()->GetTouchable());
  G4VPhysicalVolume* phys_strip = theTouchable->GetVolume();
  G4int StripNumber = 0;
  StripNumber=phys_strip->GetCopyNo();
  //G4cout << StripNumber << G4endl;

  if (HitID[StripNumber]==-1)
    { 
      //LayerHit* layerHit = new LayerHit();
      //layerHit->AddSil(edep,stepl);
      TrackerHit* trackerHit = new TrackerHit();
      //      trackerHit->SetEnergy(edep);
      // G4cout << "sono qui " << edep << G4endl;
      
      trackerHit->AddEnergy(edep);
      trackerHit->SetNumber(StripNumber);
      trackerHit->SetIPosition(aStep->GetPreStepPoint()->GetPosition());
      trackerHit->SetFPosition(aStep->GetPostStepPoint()->GetPosition());
      //TrackerCollection->insert(trackerHit) -1;
      //G4cout << "sono qua " << edep << G4endl;
      HitID[StripNumber] = TrackerCollection->insert(trackerHit) -1;
      //G4cout << " New Layer Hit "<< StripNumber <<  HitID[StripNumber] << G4endl;
    }

  else
    { 
      (*TrackerCollection)[HitID[StripNumber]]->AddEnergy(edep);
      //G4cout << " Energy added to Layer" << G4endl; 
    }
  

//   TrackerHit* trackerHit = new TrackerHit();
//   trackerHit->SetEnergy(edep);
//   trackerHit->SetIPosition(aStep->GetPreStepPoint()->GetPosition());
//   trackerHit->SetFPosition(aStep->GetPostStepPoint()->GetPosition());
//   TrackerCollection->insert(trackerHit) -1;
	  
  return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void TrackerSD::EndOfEvent(G4HCofThisEvent* HCE)
{
  static G4int HCID = -1;
  if(HCID<0)
    { 
      HCID = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]);
    }
  HCE->AddHitsCollection(HCID,TrackerCollection);
  
  for (G4int i=0;i<NSTRIPS;i++) 
    {
      HitID[i] = -1;
    };
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void TrackerSD::clear()
{} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void TrackerSD::DrawAll()
{} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void TrackerSD::PrintAll()
{} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

























