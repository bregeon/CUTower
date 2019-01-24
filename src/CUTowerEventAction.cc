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
// $Id: CUTowerEventAction.cc,v 1.1 2007/08/08 10:07:55 bregeon Exp $
// GEANT4 tag $Name:  $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "CUTowerEventAction.hh"

#include "CUTowerRunAction.hh"
#include "CUTowerEventActionMessenger.hh"

#include "G4Event.hh"
#include "G4TrajectoryContainer.hh"
#include "G4VTrajectory.hh"
#include "G4VVisManager.hh"
#include "G4UnitsTable.hh"

#include "Randomize.hh"
#include <iomanip>
#include <cmath>

#include "TrackerHit.hh"
#include "G4HCofThisEvent.hh"
#include "G4VHitsCollection.hh"
#include "G4SDManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

CUTowerEventAction::CUTowerEventAction(CUTowerRunAction* run)
:runAct(run),printModulo(1),eventMessenger(0)
{
  eventMessenger = new CUTowerEventActionMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

CUTowerEventAction::~CUTowerEventAction()
{
  delete eventMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void CUTowerEventAction::BeginOfEventAction(const G4Event* evt)
{  
 G4int evtNb = evt->GetEventID();
 if (evtNb%printModulo == 0) { 
   G4cout << "\n---> Begin of event: " << evtNb << G4endl;
   CLHEP::HepRandom::showEngineStatus();
 }
 
 // initialisation per event
  totalenergy = 0;
  for(G4int i=0;i<100;++i)
    {
      edepLayer[i] = 0;
      edepColumn[i] = 0;
      for(G4int j=0;j<100;++j)
	edepLog[i][j] = 0;
    }
  // reset num hits per tkr layer
  for(G4int l=0;l<36;l++) nHitLayer[l]=0;
  //get the Root Tree
  CUTowerRootTree *theRootTree;
  theRootTree = runAct->getTheRootTree();
  theRootTree->ResetTkrHitsVariables();

  //Create Tracker Hit Collection
 G4SDManager * SDman = G4SDManager::GetSDMpointer();  
 if (trackerCollID==-1)
   {
     trackerCollID = SDman->GetCollectionID("TrackerCollection");
     G4cout << "COLLECTION TRK" << trackerCollID  << G4endl;
   }
 


}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void CUTowerEventAction::EndOfEventAction(const G4Event* evt)
{

  // Local variables
  G4double mass = 0;
  G4double px=0, py=0, pz=0;
  G4double x0=0, y0=0, z0=0;

  //print per event (modulo n)
  //
  G4int evtNb = evt->GetEventID();
  if (evtNb%printModulo == 0) {
    G4cout << "Energy: " << totalenergy << G4endl;	
    for(G4int l=0; l<12; l++){
    G4cout <<"Layer "<<l<<": "<<edepLayer[l];
    G4cout <<"MeV Column "<<l<<": "<<edepColumn[l];
    G4cout <<"MeV"<<G4endl;
    }
    G4cout << "---> End of event: " << evtNb << G4endl;	
  }
  
  //Save data in the tree
  //get the Root Tree
  CUTowerRootTree *theRootTree;
  theRootTree = runAct->getTheRootTree();
  
  //Retrieve Event Number
  theRootTree->EventNumber = evtNb;
  
  // Get the first Primary Vertex to fill MC info
  // G4PrimaryVertex is a linked list of G4PrimaryParticle
  G4PrimaryVertex* primVtx;
  primVtx = evt->GetPrimaryVertex(); 
  //G4PrimaryParticle
  G4PrimaryParticle* primPart;
  primPart = primVtx->GetPrimary();
  x0 = primVtx->GetX0();
  y0 = primVtx->GetY0();
  z0 = primVtx->GetZ0();

  theRootTree->ParticleType = primPart->GetPDGcode();
  px = primPart->GetPx();
  py = primPart->GetPy();
  pz = primPart->GetPz();
  mass = primPart->GetMass();

  // Fill MC Variables
  theRootTree->McPosition[0] = x0;
  theRootTree->McPosition[1] = y0;
  theRootTree->McPosition[2] = z0;
  
  theRootTree->McMomentum[0] = px;
  theRootTree->McMomentum[1] = py;
  theRootTree->McMomentum[2] = pz;
  theRootTree->McEnergy = sqrt((px*px + py*py + pz*pz) + mass*mass);
 
  //Copy Calorimeter variables to Tree branch variables
  theRootTree->CalEnergy = totalenergy;
  for(G4int i=0;i<100;++i)
    {
      theRootTree->ELayer[i]    = edepLayer[i];
      theRootTree->EColumn[i] = edepColumn[i];
      for(G4int j=0;j<100;++j)
	theRootTree->ELayerColumn[i][j] = edepLog[i][j];
    }

  //Now look at Tracker Hits
  G4HCofThisEvent* HCE = evt->GetHCofThisEvent();
  TrackerHitsCollection* THC = 0;
  if (HCE)
    {
      THC = (TrackerHitsCollection*)(HCE->GetHC(trackerCollID));
    }
  if (THC)
    {
      
      int n_hit = THC->entries();
      theRootTree->TkrTotalNumHits =n_hit;
      //G4cout << "Number of tracker hits in this event =  " << n_hit << G4endl;
      G4double ESil=0;
      G4int Number = -1;
      G4int layer  = -1;
      //G4ThreeVector pos=0;
      // This is a cycle on all the tracker hits of this event

      for (int i=0;i<n_hit;i++)
	{
	  if(i>=500) break;
	  // Here we put the hit data in a an ASCII file for 
	  // later analysis 
	  ESil = (*THC)[i]->GetEnergy();
	  Number = (*THC)[i]->GetNumber();
	  //  pos = (*THC)[i]->Getpos();
	  layer = (G4int) Number/10000;
	  if (ESil/keV >28)  nHitLayer[layer] +=1; //cut at 28 keV correspondig to ~1/4 of MIP

	  theRootTree->TkrLayerHit[i]=layer;
	  theRootTree->TkrStripHit[i]=(G4int) Number%10000;
	  theRootTree->TkrEnergyHit[i]=ESil/keV;

	  if(layer == 0||layer == 3||layer == 4||layer == 7||layer == 8||layer == 11||layer == 12||layer == 15||layer == 16||layer ==19 || layer == 20 || layer ==23 ||  layer == 24||  layer == 27||  layer == 28||  layer == 31|| layer == 32||  layer == 35){ 
	    theRootTree->TkrLayerView[i]= 0  ;
	  }
	  else{
	    theRootTree->TkrLayerView[i]= 1  ;
	  }
	  /*
	    G4cout << std::setw(7) << evtNb << " " <<  
	    ESil/keV << " " <<
	    Number << " " << (G4int) Number/10000 << " " <<
 	    (*THC)[i]->GetIPosition().x()/mm <<" "<<
	    (*THC)[i]->GetIPosition().y()/mm <<" "<<
	    (*THC)[i]->GetIPosition().z()/mm <<" "<<
	    (*THC)[i]->GetFPosition().x()/mm <<" "<<
	    (*THC)[i]->GetFPosition().y()/mm <<" "<<
	    (*THC)[i]->GetFPosition().z()/mm <<" "<<
	    G4endl;
	  */
	}

    }
  
  //G4cout << "Tkr Hits per Layer" ;
  for(G4int i=0;i<36;++i)
    {
      theRootTree->TkrNumLayerHits[i] = nHitLayer[i];
      //G4cout << " " <<  nHitLayer[i];
      //G4cout << i << " " <<  nHitLayer[i] << G4endl;
      
    }
  //G4cout << G4endl;
  
  //Fill the Tree
  theRootTree->FillTree();
    
}  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void CUTowerEventAction::AddLayer(G4double enedep, G4int ilayer, G4int icolumn)
{
  totalenergy += enedep;
  if(ilayer<100)
    edepLayer[ilayer] += enedep;
  if(icolumn<100)
    edepColumn[icolumn] += enedep;
  if(ilayer<100 && icolumn<100)
    edepLog[ilayer][icolumn] += enedep;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
