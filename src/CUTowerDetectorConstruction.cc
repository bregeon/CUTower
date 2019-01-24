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
// $Id: CUTowerDetectorConstruction.cc,v 1.23 2006/06/29 17:48:54 gunter Exp $
// GEANT4 tag $Name: geant4-08-02 $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "CUTowerDetectorConstruction.hh"
#include "CUTowerDetectorMessenger.hh"

#include "G4Material.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4AssemblyVolume.hh"

#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "Randomize.hh"

#include "G4SDManager.hh"
#include "TrackerSD.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

CUTowerDetectorConstruction::CUTowerDetectorConstruction()
:defaultMaterial(0), CalMaterial(0), CalGapMaterial(0), AirMaterial(0), AlMaterial(0), SiMaterial(0),
 solidWorld(0),logicWorld(0),physiWorld(0),
 solidCalor(0),logicCalor(0),physiCalor(0),
 solidLog(0),  logicLog(0), physiLog(0),
 TrayLogical(0), TrayPhysical(0), TrackerLogical(0), TrackerPhysical(0),
 solidTKRStripX(0), logicTKRStripX(0), physiTKRStripX(0),
 solidTKRStripY(0), logicTKRStripY(0), physiTKRStripY(0)
{
  // default parameter values of the world
  WorldSizeZ  = 2.*m;
  WorldSizeXY = 2.*m;

  // default parameter values of the calorimeter
  LayerThickness    = 19.9*mm;
  LogWidth	    = 26.7*mm;
  LogLength	    = 326*mm;
  NbOfLayers        = 8;       // Can be specified in macro
  NbOfColumns       = 12;      //G4int
  LayerGap	    = 1.45*mm;
  ColumnGap	    = 1.14*mm;
  
  ComputeCalorParameters();
  
  // default parameter values of the Trays
  TrackerSideLength   = 36.*cm;
  TrackerHeight	      = 63.*cm;
  TraySideLength      = 35.1*cm;
  TrayThickness       = 3.1*cm;
  ThinHoneyComb       = 0.020*mm;
  ThickHoneyComb      = 0.060*mm;
  ThinHoneyCombHeight = 2.8*cm;
  SiliconThickness    = 0.400*mm;
  SiliconPitch	      = 0.228*mm;
  ThinTungsten        = 0.097*mm; //Mid   Bottom Side PayLoad
  ThickTungsten       = 0.630*mm; //Heavy Bottom Side PayLoad
  TrayCopyNumber      = 0;
  
  
  AlTkrWidthX0      = 1.4;	// Equivalent Al Tkr in X0
  HoneyCombType     = 1;

  // materials
  DefineMaterials();
  
  // create commands for interactive definition of the calorimeter
  detectorMessenger = new CUTowerDetectorMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

CUTowerDetectorConstruction::~CUTowerDetectorConstruction()
{ delete detectorMessenger;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* CUTowerDetectorConstruction::Construct()
{
  return ConstructCalorimeter();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void CUTowerDetectorConstruction::DefineMaterials()
{ 
 //This function illustrates the possible ways to define materials
 
G4String symbol;             //a=mass of a mole;
G4double a, z, density;      //z=mean number of protons;  
G4int ncomponents;
G4double fractionmass;

//
// define Elements
//

new G4Element("Hydrogen",symbol="H" , z= 1., a= 1.01*g/mole);
new G4Element("Carbon"  ,symbol="C" , z= 6., a= 12.01*g/mole);
G4Element* N   = new G4Element("Nitrogen",symbol="N" , z= 7., a= 14.01*g/mole);
G4Element* O   = new G4Element("Oxygen"  ,symbol="O" , z= 8., a= 16.00*g/mole);

//
// define simple materials
//
G4Material *Carbon  = new G4Material("Carbon"  ,z= 6., a= 12.01*g/mole, density= 2.265*g/cm3);
G4Material *Al   = new G4Material("Aluminium", z=13., a=26.98*g/mole, density=2.700*g/cm3);
G4Material *Si   = new G4Material("Silicon"  ,z= 14., a= 28.0*g/mole, density= 2.33*g/cm3);
G4Material *W    = new G4Material("Tungsten"  ,z= 74., a= 183.8*g/mole, density= 19.3*g/cm3);
new G4Material("liquidArgon", z=18., a= 39.95*g/mole, density= 1.390*g/cm3);
new G4Material("Lead"     , z=82., a= 207.19*g/mole, density= 11.35*g/cm3);
  // special Al for homogeneous HoneyComb
G4Material *HCThinAl  = new G4Material("Aluminium", z=13., a=26.98*g/mole, density=0.016*g/cm3); 
G4Material *HCThickAl = new G4Material("Aluminium", z=13., a=26.98*g/mole, density=0.048*g/cm3); 

 G4Material *Pb = new G4Material("Lead", z=82., a=207.2*g/mole,density=11.350*g/cm3);



//
// CsI
//
G4Element* elI = new G4Element("Iodine", "I",  53., 126.9*g/mole);
G4Element* elCs= new G4Element("Cesium", "Cs", 55., 132.9*g/mole);

G4Material *CsIMaterial = new G4Material("CsI", 4.51*g/cm3, 2);
CsIMaterial->AddElement(elI, .5);
CsIMaterial->AddElement(elCs,.5);

//
// define a material from elements.   case 2: mixture by fractional mass
//

G4Material* Air = new G4Material("Air"  , density= 1.290*mg/cm3, ncomponents=2);
Air->AddElement(N, fractionmass=0.7);
Air->AddElement(O, fractionmass=0.3);

//
// examples of vacuum
//
G4Material* Vacuum =
new G4Material("Galactic", z=1., a=1.01*g/mole,density= universe_mean_density,
                           kStateGas, 2.73*kelvin, 3.e-18*pascal);

G4Material* beam = 
new G4Material("Beam", density= 1.e-5*g/cm3, ncomponents=1,
                       kStateGas, STP_Temperature, 2.e-2*bar);
beam->AddMaterial(Air, fractionmass=1.);

G4cout << *(G4Material::GetMaterialTable()) << G4endl;

//default materials of the World
defaultMaterial  = Vacuum;
AirMaterial = Air;
AlMaterial = Al;
SiMaterial = Si;
CalMaterial = CsIMaterial;
CalGapMaterial = Vacuum;
TungstenMaterial = W;
HCThinAlMaterial  = HCThinAl;
HCThickAlMaterial = HCThickAl;

 ColMaterial = Al;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4AssemblyVolume* CUTowerDetectorConstruction::ConstructHoneyComb(G4double HoneyCombWidth)
{
   // Make the elementary assembly of the whole structure
   // Use mother volume instead of assembly as Geant4 does not 
   // support assembly of assemblies
 
   G4double xplate = 5.5*mm;      	  // Honeycomb Cell Wall is 5.5mm long
   G4double yplate = ThinHoneyCombHeight; // Honeycomb Cell Wall 2cm high
   G4double zplate = HoneyCombWidth;      // Honeycomb Cell Wall Width
   G4double dshift = sqrt(3)*xplate/2.;
 
   G4AssemblyVolume* tplate = new G4AssemblyVolume();
 
   // plate volume : Honeycomb Cell Wall 
   G4Box* plateS           = new G4Box("plateS", xplate/2., yplate/2., zplate/2);
   G4LogicalVolume* plateV = new G4LogicalVolume(plateS, AlMaterial, "PLATE");
   
   // compose assembly 
   G4ThreeVector pos0(0.,0., 0.);
   tplate->AddPlacedVolume(plateV, pos0, 0);
   
   G4RotationMatrix* rot1 = new G4RotationMatrix();
   rot1->rotateX(90.*deg);
   G4RotationMatrix *rot;
 
   // Make a hexagone cell out of 6 toothplates. These can zip togeather
   // without generating overlaps (they are self-contained)
   G4AssemblyVolume* cell = new G4AssemblyVolume();
   for (G4int i2=0; i2<3; i2++) {
     G4double phi =  60.*i2 * deg;
     G4double xp = dshift*std::sin(phi);
     G4double yp = -dshift*std::cos(phi);
     rot = new G4RotationMatrix(*rot1);
     rot->rotateZ(phi); 
     G4ThreeVector pos(xp, yp, 0.);
     cell->AddPlacedAssembly(tplate, pos, rot);
   }   
 
   // Make a row as an assembly of cells, then combine rows in a honeycomb
   // structure. This again works without any need to define rows as
   // "overlapping"
   G4AssemblyVolume* row = new G4AssemblyVolume();
   G4int ncells = 16; 
   for (G4int i3=0; i3<ncells; i3++) {
     G4double ycell = (2*i3+1)*dshift;
     G4ThreeVector pos1(0., ycell, 0.);
     row->AddPlacedAssembly(cell, pos1, 0);
     G4ThreeVector pos2(0., -ycell, 0.);
     row->AddPlacedAssembly(cell, pos2, 0);
   }
 
   // Make a tray as an assembly of rows
   G4AssemblyVolume* honeyComb = new G4AssemblyVolume();
   G4double dxrow = 3.*(dshift)*std::tan(30.*deg);
   G4double dyrow = dshift;
   G4int nrows  = 20; 
   
   for (G4int i4=0; i4<nrows; i4++) {
     G4double xrow = 0.5*(2*i4+1)*dxrow;
     G4double yrow = 0.5*dyrow;
     if ((i4%2)==0) yrow = -yrow;
     G4ThreeVector pos1(xrow, yrow, 0);
     honeyComb->AddPlacedAssembly(row, pos1, 0);
     G4ThreeVector pos2(-xrow, -yrow, 0);
     honeyComb->AddPlacedAssembly(row, pos2, 0);
     }

  return honeyComb;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4LogicalVolume* CUTowerDetectorConstruction::BuildTrackerTray(G4int trayCopyNumber, G4int SiPlanesFlag, G4double HoneyCombWidth, G4double TungstenThickness)
{
  //Tray Volume
  G4Box*           TraySolid;
  G4LogicalVolume* TrayLogical;
   
  TraySolid    = new G4Box("TraySolid", TraySideLength/2., TraySideLength/2., TrayThickness/2.);
  TrayLogical  = new G4LogicalVolume(TraySolid, AirMaterial, "TrayLogical");

  if(HoneyCombType == 2) // realistic
    {
     //Build a HoneyComb Tray
     G4AssemblyVolume *honeyComb;
     honeyComb = ConstructHoneyComb(HoneyCombWidth);
     // Make a tray with the imprint of honeyComb
     G4double randomizeX = 1.*cm*(G4UniformRand()-0.5);
     G4double randomizeY = 1.*cm*(G4UniformRand()-0.5);  
     G4ThreeVector posTray(randomizeX, randomizeY, 0);
     honeyComb->MakeImprint(TrayLogical, posTray, 0, 0);  
     }
  else if(HoneyCombType == 1) // homogeneous
    {
    // Define Material depending whether we have Thin or Thick Honey comb
    // HoneyCombHeight should also be different
    G4Material *HoneyCombAlMaterial;
    HoneyCombAlMaterial = 0;
    if(HoneyCombWidth == ThinHoneyComb)
      HoneyCombAlMaterial = HCThinAlMaterial;
    else if(HoneyCombWidth == ThickHoneyComb)
      HoneyCombAlMaterial = HCThickAlMaterial;
    
    G4Box* HoneyCombSolid = new G4Box("HoneyCombSolid", TraySideLength/2., TraySideLength/2., ThinHoneyCombHeight/2.);
    G4LogicalVolume* HoneyCombLogical  = new G4LogicalVolume(HoneyCombSolid, HoneyCombAlMaterial, "HoneyCombLogical");
    G4PVPlacement *HoneyCombPhys;
    G4ThreeVector posHoneyComb(0, 0, 0); // Center of the tray
    HoneyCombPhys = new G4PVPlacement(0, posHoneyComb, HoneyCombLogical,
    					      "HoneyCombPhys", TrayLogical, false, 0); 
    }
  else if(HoneyCombType == 0) // nothing
    {
    G4cout<<"No honey comb layer"<<G4endl;
    }    

  /* --------------------------------------------------------old ------------------------------------------------*/
  //Silicon Plane - 400mu
  //   G4Box* SiPlaneSolid    = new G4Box("SiPlaneSolid", TraySideLength/2., TraySideLength/2., SiliconThickness/2.);
  //   G4LogicalVolume* SiPlaneLogical  = new G4LogicalVolume(SiPlaneSolid, SiMaterial, "SiPlaneLogical");
  //   G4PVPlacement *SiPlanePhysUp, *SiPlanePhysDown;
  //   G4ThreeVector posSiliconUp(0, 0, (2.8/2.+ 0.02+0.1+0.02)*cm); //HalfHoneyComb + Carbon Sheet + W Sheet +  HalfSiliconPlane
  //   G4ThreeVector posSiliconDown(0, 0, -(2.8/2.+ 0.02+0.02)*cm);  //HalfHoneyComb + Carbon Sheet + HalfSiliconPlane
  //   SiPlanePhysUp   = new G4PVPlacement(0, posSiliconUp, SiPlaneLogical,
  //   					    "SiPlanePhysicalup", TrayLogical, false, 0); 
  //   SiPlanePhysDown = new G4PVPlacement(0, posSiliconDown, SiPlaneLogical,
  //   					    "SiPlanePhysicalDown", TrayLogical, false, 1); 
  /* --------------------------------------------------------old ------------------------------------------------*/
  
  // Check the SiPlane flag
  //G4cout << "Checking Si Flag for tray " << trayCopyNumber  << G4endl;
  //if (SiPlanesFlag&0x1) G4cout << "Si Flag set to enable layer 1 " << G4endl;
  //if ((SiPlanesFlag>>1)&0x1) G4cout << "Si Flag set to enable layer 2 " << G4endl;

  // Check the SiPlane orientation:
  //if (trayCopyNumber%2){G4cout << " tray " <<trayCopyNumber << " with X strips" <<G4endl;}
  //else{G4cout << " tray " <<trayCopyNumber << " with Y strips" <<G4endl;}

  G4int nStripsperplane = 1500; // number of strips (PVPlacements) in a plane

  if (trayCopyNumber%2){// Strip in 'vertical' direction, reading X
     
    //New Silicon Plane made of strips - logicTKRStripY is created once in ConstructTracker so that it can be added to SDManager
    if (SiPlanesFlag&0x1){
      // Creating TOP Si Plane; WARNING if defined for tray 0 (top) it define negative copynumber for strips 
      for (G4int i=0; i<nStripsperplane; i++)
	{      
	  G4double posX = -TraySideLength/2.+SiliconPitch/2. + i*SiliconPitch;
	  G4ThreeVector posSiliconDown(posX, 0, -(2.8/2.+ 0.02+0.02)*cm);  //HalfHoneyComb + Carbon Sheet + HalfSiliconPlane
	  physiTKRStripY =  new G4PVPlacement(0,posSiliconDown, logicTKRStripY, "TKRStripY",
					      TrayLogical, false, i+10000*(2*trayCopyNumber-1));					 
	}
    }
    if ((SiPlanesFlag>>1)&0x1){
      // Creating BOTTOM Si Plane 
      for (G4int i=0; i<nStripsperplane; i++)
	{      
	  G4double posX = -TraySideLength/2.+SiliconPitch/2. + i*SiliconPitch;
	  G4ThreeVector posSiliconUp(posX, 0, (2.8/2.+ 0.02+0.1+0.01)*cm); //HalfHoneyComb + Carbon Sheet + W Sheet +  HalfSiliconPlane
	  physiTKRStripY =  new G4PVPlacement(0,posSiliconUp, logicTKRStripY, "TKRStripY",
					      TrayLogical, false, i+10000*2*trayCopyNumber);					 
	}
    }
    
  }
  else{  // Strip in 'horizontal' direction, reading Y
  
  //New Silicon Plane made of strips - logicTKRStripX is created once in ConstructTracker so that it can be added to SDManager
  if (SiPlanesFlag&0x1){
    // Creating TOP Si Plane; WARNING if defined for tray 0 (top) it define negative copynumber for strips 
    for (G4int i=0; i<nStripsperplane; i++)
      {      
	G4double posY = -TraySideLength/2.+SiliconPitch/2. + i*SiliconPitch;
	G4ThreeVector posSiliconDown(0, posY, -(2.8/2.+ 0.02+0.02)*cm);  //HalfHoneyComb + Carbon Sheet + HalfSiliconPlane
	physiTKRStripX =  new G4PVPlacement(0,posSiliconDown, logicTKRStripX, "TKRStripX",
					    TrayLogical, false, i+10000*(2*trayCopyNumber-1));					 
      }
  }
  if ((SiPlanesFlag>>1)&0x1){
    // Creating BOTTOM Si Plane 
    for (G4int i=0; i<nStripsperplane; i++)
      {      
	G4double posY = -TraySideLength/2.+SiliconPitch/2. + i*SiliconPitch;
	G4ThreeVector posSiliconUp(0, posY, (2.8/2.+ 0.02+0.1+0.01)*cm); //HalfHoneyComb + Carbon Sheet + W Sheet +  HalfSiliconPlane
	physiTKRStripX =  new G4PVPlacement(0,posSiliconUp, logicTKRStripX, "TKRStripX",
					    TrayLogical, false, i+10000*2*trayCopyNumber);					 
      }
  }
  }

  //Tungsten Plane 
  if(TungstenThickness>0)
    {
    G4Box* TungstenPlaneSolid    = new G4Box("TungstenPlaneSolid", TraySideLength/2., TraySideLength/2., TungstenThickness/2.);
    G4LogicalVolume* TungstenPlaneLogical  = new G4LogicalVolume(TungstenPlaneSolid, TungstenMaterial, "TungstenPlaneLogical");
    G4PVPlacement *TungstenPlanePhysDown;
    G4ThreeVector posTungsten(0, 0, (2.8/2.+ 0.02+0.05)*cm); //HalfHoneyComb + Carbon Sheet + Half W Sheet   
    TungstenPlanePhysDown = new G4PVPlacement(0, posTungsten, TungstenPlaneLogical,
    					      "TungstenPlanePhyTungstencalDown", TrayLogical, false, 0); 
    }


  return TrayLogical;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4LogicalVolume* CUTowerDetectorConstruction::ConstructTracker(G4int nTopTrays,G4int nMidTrays, G4int nHeavyTrays,
							     G4int nLightTrays, G4int nBottomTrays)
{
  //Tracker Volume
  G4Box*           TrackerSolid;
  G4LogicalVolume* TrackerLogical;
  
  TrackerSolid    = new G4Box("TrackerSolid", TrackerSideLength/2., TrackerSideLength/2., TrackerHeight/2.);
  TrackerLogical  = new G4LogicalVolume(TrackerSolid, defaultMaterial, "TrackerLogical");

  //Strip Volume X
  solidTKRStripX = new G4Box("TKRStripX", TraySideLength/2., SiliconPitch/2., SiliconThickness/2.); 
  logicTKRStripX = new G4LogicalVolume(solidTKRStripX, SiMaterial, "TKRStripX");
  //Strip Volume Y
  solidTKRStripY = new G4Box("TKRStripY", SiliconPitch/2., TraySideLength/2., SiliconThickness/2.); 
  logicTKRStripY = new G4LogicalVolume(solidTKRStripY, SiMaterial, "TKRStripY");
    
    //Make a Tracker Tower with many trays
  G4double Zstart  = -TrackerHeight/2. +  TrayThickness/2. + 1*cm;
  G4double dZtrays = 3.2*cm;
  G4int trayCopyNumber = 0; 

  // Top tray        
  for(G4int i=0; i<nTopTrays; i++)
    {
     //Create a new Tray 
     TrayLogical = BuildTrackerTray(trayCopyNumber, 2, ThinHoneyComb, ThinTungsten);  //ThinTungsten
     char trayName[20];
     sprintf(trayName, "TrayPhysical_%i",i);
  
     G4ThreeVector trayPosition(0, 0, Zstart + trayCopyNumber*dZtrays);
     TrayPhysical = new G4PVPlacement(0, trayPosition, TrayLogical,
  				      trayName, TrackerLogical, false, trayCopyNumber); 
     trayCopyNumber++;				      
    }
  
  // Mid Trays + a Top tray        
  for(G4int i=0; i<nMidTrays; i++)
    {
     //Create a new Tray 
     TrayLogical = BuildTrackerTray(trayCopyNumber, 3, ThinHoneyComb, ThinTungsten);  //ThinTungsten
     char trayName[20];
     sprintf(trayName, "TrayPhysical_%i",i);
  
     G4ThreeVector trayPosition(0, 0, Zstart + trayCopyNumber*dZtrays);
     TrayPhysical = new G4PVPlacement(0, trayPosition, TrayLogical,
  				      trayName, TrackerLogical, false, trayCopyNumber); 
     trayCopyNumber++;				      
    }
    
  // Heavy Trays								       				       
  for(G4int i=0; i<nHeavyTrays; i++)						       
    {										       
     //Create a new Tray 							       
     TrayLogical = BuildTrackerTray(trayCopyNumber, 3, ThickHoneyComb, ThickTungsten);  //ThinTungsten    
     char trayName[20]; 							       
     sprintf(trayName, "TrayPhysical_%i",i);					       
  
     G4ThreeVector trayPosition(0, 0, Zstart + trayCopyNumber*dZtrays);				       
     TrayPhysical = new G4PVPlacement(0, trayPosition, TrayLogical,		       
  				      trayName, TrackerLogical, false, trayCopyNumber);  
     trayCopyNumber++;  			      				       
    }										       
  							       
  // Light Trays								       				       
  for(G4int i=0; i<nLightTrays; i++)						       
    {										       
     //Create a new Tray 							       
     TrayLogical = BuildTrackerTray(trayCopyNumber, 3, ThinHoneyComb, 0);  //ThinTungsten    
     char trayName[20]; 							       
     sprintf(trayName, "TrayPhysical_%i",i);					       
  
     G4ThreeVector trayPosition(0, 0, Zstart + trayCopyNumber*dZtrays);				       
     TrayPhysical = new G4PVPlacement(0, trayPosition, TrayLogical,		       
  				      trayName, TrackerLogical, false, trayCopyNumber);  
     trayCopyNumber++;  			      				       
    }										       

  // Bottom Trays - To Be Reviewed								       							       
  for(G4int i=0; i<nBottomTrays; i++)						       
    {										       
     //Create a new Tray 							       
     TrayLogical = BuildTrackerTray(trayCopyNumber, 1, ThickHoneyComb, 0);  //ThinTungsten    
     char trayName[20]; 							       
     sprintf(trayName, "TrayPhysical_%i",i);					       
  
     G4ThreeVector trayPosition(0, 0, Zstart + trayCopyNumber*dZtrays);				       
     TrayPhysical = new G4PVPlacement(0, trayPosition, TrayLogical,		       
  				      trayName, TrackerLogical, false, trayCopyNumber);  
     trayCopyNumber++;  			      				       
    }										       
  SetTrayCopyNumber(trayCopyNumber);
  return TrackerLogical;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* CUTowerDetectorConstruction::ConstructCalorimeter()
{
  
  // Clean old geometry, if any
  //
  G4cout<<"Cleaning Up Geometry"<<G4endl;
  G4GeometryManager::GetInstance()->OpenGeometry();
  G4PhysicalVolumeStore::GetInstance()->Clean();
  G4LogicalVolumeStore::GetInstance()->Clean();
  G4SolidStore::GetInstance()->Clean();

  // complete the Calor parameters definition
  ComputeCalorParameters();
   
  //     
  // World
  //
  solidWorld = new G4Box("World",				//its name
                   WorldSizeXY/2,WorldSizeXY/2,WorldSizeZ/2);	//its size
                         
  logicWorld = new G4LogicalVolume(solidWorld,		//its solid
                                   defaultMaterial,	//its material
                                   "World");		//its name
                                   
  physiWorld = new G4PVPlacement(0,			//no rotation
  				 G4ThreeVector(),	//at (0,0,0)
                                 logicWorld,		//its logical volume				 
                                 "World",		//its name
                                 0,			//its mother  volume
                                 false,			//no boolean operation
                                 0);			//copy number

  //                               
  // Calorimeter
  //  
  solidCalor=0; logicCalor=0; physiCalor=0;
  solidLog = 0;  logicLog = 0; physiLog = 0;

  if (CalorThickness > 0.)
    {
    G4cout<<"Building a Calorimeter : CalorThickness = "<<CalorThickness<<"mm"<<G4endl;
    
    solidCalor = new G4Box("Calorimeter", (LogLength+ColumnGap*12+10*cm)/2., (LogLength+ColumnGap*12+10*cm)/2., (CalorThickness+LayerGap*8+4.*cm)/2.); 												 
  
    logicCalor = new G4LogicalVolume(solidCalor,  //its solid															 
    				     CalGapMaterial,	  //its material														 
    				     "Calorimeter");	  //its name														 
    					 																	 
    physiCalor = new G4PVPlacement(0,			  //no rotation 													 
    				   G4ThreeVector(0,0,0),  //at (0,0,0)														 
    				   logicCalor,    //its logical volume														 
    				   "Calorimeter", //its name															 
    				   logicWorld,    //its mother  volume														 
    				   false,	  //no boolean operation													 
    				   0);  	  //copy number 														 
  
    //  			        																	 
    // Layers Columns
    //																						 
   G4RotationMatrix* rotateMatrix = new G4RotationMatrix();
   rotateMatrix -> rotateZ(90.0*deg);
   
   double Xpos=0., Ypos=0., Zpos=0.;   
   
   solidLog = new G4Box("Calorimeter", LogWidth/2., LogLength/2., LayerThickness/2.);	   
   logicLog = new G4LogicalVolume(solidLog, CalMaterial, CalMaterial->GetName());
 
   for(int i=0;i<NbOfLayers;i+=2)
     {													 
     for(int j=0;j<NbOfColumns;++j)																		 
       {
	Xpos = -((NbOfColumns-1)/2)*(LogWidth+ColumnGap) + j*(LogWidth+ColumnGap);
	Ypos = 0.;
	Zpos = -((NbOfLayers-1)/2)*(LayerThickness+LayerGap) + i*(LayerThickness+LayerGap);
    	physiLog = new G4PVPlacement(0, G4ThreeVector(Xpos, Ypos, Zpos),
    				     logicLog, "Layer", logicCalor, false, 100*i+j);
     
	Xpos = 0.;
	Ypos = -((NbOfColumns-1)/2)*(LogWidth+ColumnGap) + j*(LogWidth+ColumnGap);
	Zpos = -((NbOfLayers-1)/2)*(LayerThickness+LayerGap) + (i+1)*(LayerThickness+LayerGap);
    	physiLog = new G4PVPlacement(rotateMatrix, G4ThreeVector(Xpos, Ypos, Zpos),
    				     logicLog, "Layer", logicCalor, false, 100*(i+1)+j);
	}
      }
    } //End of Calorimeter                                    
  else{G4cout<<"CalorThickness == 0. : No Calorimeter."<<G4endl;}
  
  
  //Construct Tracker and place it in the world logical					               
  G4double AlTkrWidth_cm = AlTkrWidthX0*8.9*cm;
  G4cout<<"Al Tkr Plate is now "<<AlTkrWidth_cm<<" cm."<<G4endl; 		               
  
    // Realistic LAT Tracker Tower
  if(AlTkrWidth_cm < 0.01)
    {
     // ConstructTracker(G4int nTopTrays, G4int nMidTrays, G4int nHeavyTrays, G4int nLightTrays, G4int nBottomTrays)            
     TrackerLogical = ConstructTracker(1,11, 4, 2, 1);						                 
     G4cout<<"Built a Tracker Tower with "<<GetTrayCopyNumber()<<" trays."<<G4endl; 		               
    }
  else    // Al Plate of Equivalent Radiation Length
    {
    G4Box *AlTkrPlateSolid    = new G4Box("AlTkrPlateSolid",		
                   	 TrackerSideLength/2,TrackerSideLength/2,AlTkrWidth_cm/2);
    TrackerLogical  = new G4LogicalVolume(AlTkrPlateSolid,
    					    AlMaterial,
    					    "TrackerLogical");
    G4cout<<"Built a fake Tracker Tower of "<<AlTkrWidthX0<<" X0."<<G4endl; 		               
    }
   
   G4double ZTrackerPos = -CalorThickness/2. - TrackerHeight/2. - 5.*cm;				     
   G4ThreeVector trackerPosition(0, 0, ZTrackerPos);							     
   TrackerPhysical = new G4PVPlacement(0, trackerPosition, TrackerLogical,				     
   				       "TrackerPhysical", logicWorld, false, 0);			     

//    Collimator simulation
//    G4double Ext_x = 5.*mm; //7.5 * mm;
//    G4double Ext_y = 5.*mm;//7.5 *mm;
//    G4double Ext_z = 0.2 * 8.9*cm; //1.6 * 8.9*cm; // X0 = 8.9*cm
//    G4double MedA_x = 4.1*mm;//5. * mm;
//    G4double MedA_y = 4.1*mm;//5. *mm;
//    G4double MedA_z = 0.2 * 8.9*cm; //1.6 * 8.9*cm; // X0 = 8.9*cm
//    G4double Med_x = 4.1*mm;//5. * mm;
//    G4double Med_y = 4.1*mm;//5. *mm;
//    G4double Med_z = 0.1 * 8.9*cm; //1.4 * 8.9*cm; // X0 = 8.9*cm
//    G4double IntA_x = 2.8*mm;//2.5 * mm;
//    G4double IntA_y = 2.8*mm;//2.5 *mm;
//    G4double IntA_z = 0.1 * 8.9*cm; //1.4 * 8.9*cm; // X0 = 8.9*cm
//    G4double Int_x = 2.8*mm;//2.5 * mm;
//    G4double Int_y = 2.8*mm;//2.5 *mm;
//    G4double Int_z = 0.1 * 8.9*cm; //1.2 * 8.9*cm; // X0 = 8.9*cm
// 
//    // Collimator volumes
// 
//    G4Box* ExternalCollimator = new G4Box("ExternalCollim", Ext_x/2., Ext_y/2., Ext_z/2.);
//    G4Box* MediumCollimatorAir = new G4Box("MediumCollimAir", MedA_x/2., MedA_y/2., MedA_z/2.);
//    G4Box* MediumCollimator = new G4Box("MediumCollim", Med_x/2., Med_y/2., Med_z/2.);
//    G4Box* InternalCollimatorAir = new G4Box("InternalCollimAir", IntA_x/2., IntA_y/2., IntA_z/2.);
//    G4Box* InternalCollimator = new G4Box("InternalCollim", Int_x/2., Int_y/2., Int_z/2.);
// 
//    // Logical volume
// 
//     G4LogicalVolume* ExtCollLogical  = new G4LogicalVolume(ExternalCollimator, ColMaterial, "ExtCollim");
//     G4LogicalVolume* MedCollAirLogical  = new G4LogicalVolume(MediumCollimatorAir, AirMaterial, "MedCollimAir");
//     G4LogicalVolume* MedCollLogical  = new G4LogicalVolume(MediumCollimator, ColMaterial, "MedCollim");
//     G4LogicalVolume* IntCollAirLogical  = new G4LogicalVolume(InternalCollimatorAir, AirMaterial, "IntCollimAir");
//     G4LogicalVolume* IntCollLogical  = new G4LogicalVolume(InternalCollimator, ColMaterial, "IntCollim");
// 
//    // Physical Volume
// 
//    G4double CollPos_x = 0. * mm;
//    G4double CollPos_y = 0. *mm;
//    G4double CollPos_z = -60. * cm; // TBD
// 
//    G4ThreeVector ExtCollPosition(CollPos_x, CollPos_y, CollPos_z);
//    G4ThreeVector zeroPosition(CollPos_x, CollPos_y, CollPos_z);
//    
//    G4VPhysicalVolume*  ExtCollPhysical    = new G4PVPlacement(0, ExtCollPosition, ExtCollLogical, "ExternalCollimatorPhys", logicWorld, false, 0); 
//    G4VPhysicalVolume*  MedCollAirPhysical = new G4PVPlacement(0, zeroPosition, MedCollAirLogical, "MedCollimatorAirPhys", ExtCollLogical, false, 0); 
//    G4VPhysicalVolume*  MedCollPhysical    = new G4PVPlacement(0, zeroPosition, MedCollLogical, "MedCollimatorPhys", MedCollAirLogical, false, 0); 
//    G4VPhysicalVolume*  IntCollAirPhysical = new G4PVPlacement(0, zeroPosition, IntCollAirLogical, "IntCollimatorAirPhys", MedCollLogical, false, 0); 
//    G4VPhysicalVolume*  IntCollPhysical    = new G4PVPlacement(0, zeroPosition, IntCollLogical, "IntCollimatorPhys", IntCollAirLogical, false, 0); 


  //                                        
  // Visualization attributes
  //
  logicWorld->SetVisAttributes (G4VisAttributes::Invisible);

  G4VisAttributes* simpleBoxVisAtt= new G4VisAttributes(G4Colour(1.0,1.0,1.0));
  simpleBoxVisAtt->SetVisibility(true);
  logicCalor->SetVisAttributes(simpleBoxVisAtt);


 //------------------------------
 // Sensitive Detector Manager
 //------------------------------  
 G4SDManager* SDman = G4SDManager::GetSDMpointer();
 TrackerSD* trackerSD = new TrackerSD("TrackerSD");
 SDman->AddNewDetector( trackerSD );
 if (logicTKRStripX)  
   logicTKRStripX->SetSensitiveDetector(trackerSD);
 if (logicTKRStripY)  
   logicTKRStripY->SetSensitiveDetector(trackerSD);
 G4cout << "SDManager Created" << G4endl;
  
  //
  //always return the physical World
  //
  G4cout<<"Detector Construction Done."<<G4endl;
  G4cout<<"-----------------------------------------------------\n\n"<<G4endl;
  return physiWorld;
}



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void CUTowerDetectorConstruction::SetNbOfLayers(G4int val)
{
  NbOfLayers = val;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void CUTowerDetectorConstruction::SetNbOfColumns(G4int val)
{
  NbOfColumns = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void CUTowerDetectorConstruction::SetLayerGap(G4double val)
{
  LayerGap = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void CUTowerDetectorConstruction::SetColumnGap(G4double val)
{
  ColumnGap = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void CUTowerDetectorConstruction::SetAlTkrWidthX0(G4double val)
{
  AlTkrWidthX0 = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void CUTowerDetectorConstruction::SetHoneyCombType(G4int val)
{
  HoneyCombType = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4RunManager.hh"

void CUTowerDetectorConstruction::UpdateGeometry()
{
  G4RunManager::GetRunManager()->DefineWorldVolume(ConstructCalorimeter());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
