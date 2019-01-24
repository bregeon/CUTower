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
// $Id: CUTowerRootTree.hh,v 1.1 2007/08/08 10:07:53 bregeon Exp $
// GEANT4 tag $Name:  $

#ifndef CUTowerRootTree_h
#define CUTowerRootTree_h 1

//ROOT includes
#include "TFile.h"
#include "TTree.h"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
#include "globals.hh"
#include <vector>
#include "G4DynamicParticle.hh"
#include "G4VPhysicalVolume.hh"
#include "G4DataVector.hh"
#include "G4Track.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class CUTowerRootTree
  {
	public:
		CUTowerRootTree();
		~CUTowerRootTree();
		
		//TFile *GetRootFile(){return rootOutFile;};
		TFile* OpenFile();
		TTree* CreateTree();
		void  FillTree();
		void  Save();
		void  SetFileNameByTime();
		void  SetFileNameByEnv();
		
				
		G4int EventNumber;
		G4int ParticleType;		
		G4double McEnergy;
		G4double McMomentum[3];
		G4double McPosition[3];
		G4double CalEnergy;
   		G4double ELayer[100];
 		G4double EColumn[100];
 		G4double ELayerColumn[100][100];

                G4int TkrNumLayerHits[36];
                G4int TkrTotalNumHits;
                G4int TkrLayerView[500];
                G4int TkrLayerHit[500];
                G4int TkrStripHit[500];
                G4double TkrEnergyHit[500];
                void ResetTkrHitsVariables();
    

				
	private:	
		G4int verbose;
		G4String rootFileName;
		G4String rootTreeName;
		TTree *calTree;
		TFile *rootOutFile;
};
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
#endif

