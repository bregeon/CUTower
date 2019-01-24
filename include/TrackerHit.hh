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
// $Id: TrackerHit.hh,v 1.1 2007/08/08 10:07:54 bregeon Exp $
// GEANT4 tag $Name:  $
// ------------------------------------------------------------
//      GEANT 4 class header file
//      CERN Geneva Switzerland
//
//
//      ------------ TrackerHit  ------
//           by R.Giannitrapani, F.Longo & G.Santin (13 nov 2000)
//
// ************************************************************
// This Class describe the hits on the Tracker

#ifndef TrackerHit_h
#define TrackerHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class TrackerHit : public G4VHit
{
public:
  
  TrackerHit();
  ~TrackerHit();
  TrackerHit(const TrackerHit&);
  const TrackerHit& operator=(const TrackerHit&);
  int operator==(const TrackerHit&) const;
  
  inline void* operator new(size_t);
  inline void  operator delete(void *aHit);

  void Draw();
  void Print();

private:
  
  G4double EdepSil;  // Energy deposited on the silicon strip
  G4int Number;  // Strip address
  G4ThreeVector Ipos; // Position of the hit 
  G4ThreeVector Fpos; // Position of the hit 

public:
  
  inline void SetEnergy(G4double de) {EdepSil = de;};
  inline void SetNumber(G4int N) {Number = N;};
  inline void AddEnergy(G4double de) {EdepSil += de;};
  inline void SetIPosition(G4ThreeVector xyz){ Ipos = xyz; }
  inline void SetFPosition(G4ThreeVector xyz){ Fpos = xyz; }
  
  inline G4double GetEnergy()     { return EdepSil; };
  inline G4int GetNumber()     { return Number; };
  inline G4ThreeVector GetIPosition() { return Ipos; };
  inline G4ThreeVector GetFPosition() { return Fpos; };
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

typedef G4THitsCollection<TrackerHit> TrackerHitsCollection;
extern G4Allocator<TrackerHit> TrackerHitAllocator;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void* TrackerHit::operator new(size_t)
{
  void* aHit;
  aHit = (void*) TrackerHitAllocator.MallocSingle();
  return aHit;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void TrackerHit::operator delete(void* aHit)
{
  TrackerHitAllocator.FreeSingle((TrackerHit*) aHit);
}

#endif









