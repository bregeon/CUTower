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
//---------------------------------------------------------------------------
//
// ClassName:   Histo - Generic histogram/ntuple manager class
//
//
// Author:      V.Ivanchenko 30.10.03
//
//----------------------------------------------------------------------------
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
#include <time.h>
#include "CUTowerRootTree.hh"


CUTowerRootTree::CUTowerRootTree(){
				verbose = 0;
				rootFileName = "output.root";
				rootTreeName = "CalEvent";
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

CUTowerRootTree::~CUTowerRootTree()
{
  G4cout<<"RootTree Destructor"<<G4endl;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void CUTowerRootTree::SetFileNameByTime()
{
  char toto[100];
  time_t temps;
  tm * ptm;
  time ( &temps );
  ptm = gmtime ( &temps );
  sprintf (toto,"G4v8.2_CUTower_%2dh%02d%02d.root", ptm->tm_hour, ptm->tm_min, ptm->tm_sec);
  rootFileName = toto;
  G4cout<<"Output File is "<<rootFileName<<G4endl;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void CUTowerRootTree::SetFileNameByEnv()
{
  rootFileName = "$G4OUTFILENAME.root";
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

TFile* CUTowerRootTree::OpenFile()
{
  rootOutFile = new TFile(rootFileName, "RECREATE");
  return rootOutFile;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

TTree* CUTowerRootTree::CreateTree()
{
  rootOutFile->cd();
  calTree = new TTree(rootTreeName, rootTreeName);
  calTree->Branch("EventNumber",  &EventNumber,   "EventNumber/I");
  calTree->Branch("ParticleType", &ParticleType,  "ParticleType/I");
  calTree->Branch("McEnergy",     &McEnergy,      "McEnergy/D");
  calTree->Branch("McMomentum",    McMomentum,    "McMomentum[3]/D");
  calTree->Branch("McPosition",    McPosition,    "McPosition[3]/D");
  calTree->Branch("CalEnergy",    &CalEnergy,     "CalEnergy/D");
  calTree->Branch("ELayer",        ELayer,        "ELayer[100]/D");
  calTree->Branch("EColumn",       EColumn,       "EColumn[100]/D");
  calTree->Branch("ELayerColumn",ELayerColumn,"ELayerColumn[100][100]/D");  
  calTree->Branch("TkrNumLayerHits",TkrNumLayerHits,"TkrNumLayerHits[36]/I");
  calTree->Branch("TkrTotalNumHits",&TkrTotalNumHits,"TkrTotalNumHits/I");
  calTree->Branch("TkrLayerHit",TkrLayerHit,"TkrLayerHit[500]/I");
  calTree->Branch("TkrLayerView",TkrLayerView,"TkrLayerView[500]/I");
  calTree->Branch("TkrStripHit",TkrStripHit,"TkrStripHit[500]/I");
  calTree->Branch("TkrEnergyHit",TkrEnergyHit,"TkrEnergyHit[500]/D");

  return calTree;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void CUTowerRootTree::FillTree()
{
  calTree->Fill();
  //G4cout<<"Tree Filled"<<G4endl;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void CUTowerRootTree::Save()
{
  rootOutFile->Write();
  G4cout<<"File Saved"<<G4endl;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
void CUTowerRootTree:: ResetTkrHitsVariables()
{
  for (G4int i=0;i<500;i++){
    TkrLayerHit[i]  = -1;
    TkrStripHit[i]  = -1;
    TkrLayerView[i] = -1;
    TkrEnergyHit[i] = -1;
  }


}



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
