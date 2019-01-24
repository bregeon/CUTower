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
// $Id: CUTowerDetectorConstruction.hh,v 1.7 2006/06/29 17:48:32 gunter Exp $
// GEANT4 tag $Name: geant4-08-02 $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef CUTowerDetectorConstruction_h
#define CUTowerDetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"
#include "G4AssemblyVolume.hh"

class G4Box;
class G4Tubs;
class G4LogicalVolume;
class G4VPhysicalVolume;
class G4Material;
class CUTowerDetectorMessenger;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class CUTowerDetectorConstruction : public G4VUserDetectorConstruction
{
  public:
  
    CUTowerDetectorConstruction();
   ~CUTowerDetectorConstruction();

  public:
          
     void SetNbOfLayers (G4int);         
     void SetNbOfColumns  (G4int);        
     void SetLayerGap   (G4double);   
     void SetColumnGap   (G4double);   
     
     void SetAlTkrWidthX0   (G4double);   
     void SetHoneyCombType  (G4int);   

     G4VPhysicalVolume* Construct();

     void UpdateGeometry();
     
  public:
  
     void PrintCalorParameters(); 
                    
     G4double GetWorldSizeZ()        {return WorldSizeZ;}; 
     G4double GetWorldSizeXY()       {return WorldSizeXY;};
     
     // Calorimeter
     G4double GetCalorThickness()    {return CalorThickness;}; 
     G4int    GetNbOfLayers()        {return NbOfLayers;}; 
     G4int    GetNbOfColumns()     {return NbOfColumns;}; 
     
     G4Material* GetCalMaterial()    {return CalMaterial;};
     G4double    GetLayerThickness() {return LayerThickness;};

     // Tracker
     G4double GetTrayCopyNumber()	 {return TrayCopyNumber;};
     G4double SetTrayCopyNumber(G4int n) {TrayCopyNumber = n;
     					  return TrayCopyNumber;};


     const G4VPhysicalVolume* GetphysiWorld() {return physiWorld;};           
                 
  private:
     
     G4Material*        defaultMaterial;
     G4Material*        CalMaterial;
     G4Material*        CalGapMaterial;
     G4Material*        AirMaterial;
     G4Material*        AlMaterial;
     G4Material*        SiMaterial;
     G4Material*        TungstenMaterial;
     G4Material*        HCThinAlMaterial;
     G4Material*        HCThickAlMaterial;
  G4Material*        ColMaterial;     

     // Calorimeter
     G4int              NbOfLayers;
     G4int              NbOfColumns;
     G4double           LayerThickness;          
     G4double           CalorThickness;
     G4double		LayerGap;
     G4double		ColumnGap;
     G4double           LogWidth;          
     G4double           LogLength;          

     // Tracker
     G4double 		TrackerSideLength;
     G4double 		TrackerHeight ;  

     G4double 		TraySideLength; 
     G4double 		TrayThickness; 
     G4double 		ThinHoneyComb; 
     G4double 		ThickHoneyComb; 
     G4double 		ThinHoneyCombHeight; 
     G4double 		SiliconThickness;
     G4double		SiliconPitch;
     G4double 		ThinTungsten; 
     G4double 		ThickTungsten; 
     G4double		TrayCopyNumber;

     G4double		AlTkrWidthX0;
     G4int		HoneyCombType;
     
     //World
     G4double           WorldSizeXY;
     G4double           WorldSizeZ;
            
     G4Box*             solidWorld;    //pointer to the solid World 
     G4LogicalVolume*   logicWorld;    //pointer to the logical World
     G4VPhysicalVolume* physiWorld;    //pointer to the physical World
         
     // Calorimeter
     G4Box*            solidCalor;    //pointer to the solid Calor 
     G4LogicalVolume*   logicCalor;    //pointer to the logical Calor
     G4VPhysicalVolume* physiCalor;    //pointer to the physical Calor

     G4Box*		solidLog;    //pointer to the solid Layer 
     G4LogicalVolume*	logicLog;    //pointer to the logical Layer
     G4VPhysicalVolume *physiLog;    //pointer to the physical Layer

     // Tracker
     G4LogicalVolume*   TrayLogical;    //pointer to the logical Tray
     G4VPhysicalVolume* TrayPhysical;    //pointer to the physical Tray

     G4LogicalVolume*   TrackerLogical;    //pointer to the logical Tracker 
     G4VPhysicalVolume* TrackerPhysical;    //pointer to the physical Tracker

     G4Box*             solidTKRStripX;
     G4LogicalVolume*   logicTKRStripX;    //pointer to the logical Tracker Strips
     G4VPhysicalVolume* physiTKRStripX;    //pointer to the physical Tracker Strips
     G4Box*             solidTKRStripY;
     G4LogicalVolume*   logicTKRStripY;    //pointer to the logical Tracker Strips
     G4VPhysicalVolume* physiTKRStripY;    //pointer to the physical Tracker Strips     
     CUTowerDetectorMessenger* detectorMessenger;  //pointer to the Messenger
      
  private:
    
     void DefineMaterials();
     void ComputeCalorParameters();
     G4VPhysicalVolume* ConstructCalorimeter();     
     G4AssemblyVolume*  ConstructHoneyComb(G4double HoneyCombWidth);     
     G4LogicalVolume*   BuildTrackerTray(G4int trayCopyNumber, G4int SiPlanesFlag, G4double HoneyCombWidth, G4double TungstenThickness);
     G4LogicalVolume*	ConstructTracker(G4int nTopTrays,G4int nMidTrays, G4int nHeavyTrays,
					 G4int nLightTrays, G4int nBottomTrays);
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline void CUTowerDetectorConstruction::ComputeCalorParameters()
{
  if(NbOfLayers==0)NbOfColumns=0;  
  // Compute derived parameters of the calorimeter
  CalorThickness = NbOfLayers*LayerThickness;   
  G4cout<<"\n\n-----------------------------------------------------"<<G4endl;
  G4cout<<"Starting Detector Construction"<<G4endl;
  G4cout<<"Calorimeter : "<<NbOfLayers<<" Layers and ";
  G4cout<<NbOfColumns<<" Columns"<<G4endl;

  solidLog = 0;  logicLog = 0; physiLog = 0;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

