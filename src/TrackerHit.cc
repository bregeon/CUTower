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
// $Id: TrackerHit.cc,v 1.1 2007/08/08 10:07:55 bregeon Exp $
// GEANT4 tag $Name:  $
// ------------------------------------------------------------
//      GEANT 4 class implementation file
//      CERN Geneva Switzerland
//
//
//      ------------ TrackerHit  ------
//           by R.Giannitrapani, F.Longo & G.Santin (13 nov 2000)
//
// ************************************************************

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "TrackerHit.hh"

G4Allocator<TrackerHit> TrackerHitAllocator;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

TrackerHit::TrackerHit()
{
  EdepSil = 0.; 
  Ipos = 0.;
  Fpos = 0.;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

TrackerHit::~TrackerHit()
{;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

TrackerHit::TrackerHit(const TrackerHit& right)
{
  EdepSil = right.EdepSil; 
  Ipos = right.Ipos;
  Fpos = right.Fpos;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

const TrackerHit& TrackerHit::operator=(const TrackerHit& right)
{
  EdepSil = right.EdepSil; 
  Ipos =right.Ipos;
  Fpos =right.Fpos;
  return *this;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

int TrackerHit::operator==(const TrackerHit& right) const
{
  return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void TrackerHit::Draw()
{;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void TrackerHit::Print()
{;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....











