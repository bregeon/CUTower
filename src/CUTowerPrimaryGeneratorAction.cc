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
// $Id: CUTowerPrimaryGeneratorAction.cc,v 1.2 2007/08/20 10:19:49 bregeon Exp $
// GEANT4 tag $Name:  $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "CUTowerPrimaryGeneratorAction.hh"

#include "CUTowerDetectorConstruction.hh"
#include "CUTowerPrimaryGeneratorMessenger.hh"

#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "Randomize.hh"
#include "CLHEP/Random/RandFlat.h"

extern std::ifstream inFile;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

CUTowerPrimaryGeneratorAction::CUTowerPrimaryGeneratorAction(
                                             CUTowerDetectorConstruction* CUTowerDC)
  :CUTowerDetector(CUTowerDC),rndmFlag("off"), beamFlag("off"),
   runNo(0), seedNo(0)
{
  G4int n_particle = 1;
  particleGun  = new G4ParticleGun(n_particle);
  
  //create a messenger for this class
  gunMessenger = new CUTowerPrimaryGeneratorMessenger(this);

  // default particle kinematic

  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4String particleName;
  G4ParticleDefinition* particle
    = particleTable->FindParticle(particleName="e-");
  particleGun->SetParticleDefinition(particle);
  particleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0., 1.));
  particleGun->SetParticleEnergy(10.*GeV);
  
  G4double position = -0.5*(CUTowerDetector->GetWorldSizeZ());
  particleGun->SetParticlePosition(G4ThreeVector(0.*cm,0.*cm, position));  

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

CUTowerPrimaryGeneratorAction::~CUTowerPrimaryGeneratorAction()
{
  //delete particleGun;
  delete gunMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void CUTowerPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  //this function is called at the begining of event
  //first set the seed according to run number and event number
  //inspired from beamtest06 code     
  seedNo = 100000 * ((runNo+1) % 20000) + anEvent->GetEventID();
  // G4cout << "Random Seed : " << seedNo << "\tRun : " << runNo << G4endl;
  CLHEP::HepRandom::setTheSeed(seedNo); 
  
  // particle name
  G4String particleName = gunParticleName;
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition* particle = particleTable->FindParticle(particleName);

  // particle energy
  G4double energy = gunEnergy;

  // position
  G4double x0 = gunPosition[0];
  G4double y0 = gunPosition[1];
  G4double z0 = gunPosition[2];  
  
  // beam width
  G4double beamWidth = gunWidth;

//  G4cout << beamFlag << G4endl;
  if (beamFlag == "off") 
    {
      
      // random position
      if (rndmFlag == "on")
	{
	  x0 += beamWidth*(G4UniformRand()-0.5);
	  y0 += beamWidth*(G4UniformRand()-0.5);
	} 
      
      // Beam Divergence
      G4double a_0=0, a_1=0, a_2=1, theta=0, phi=0;
      if (gunDiv != 0)
	{
	  theta = G4RandGauss::shoot(0.,gunDiv);
	  phi = G4UniformRand()*2*M_PI;
	  a_0 = sin(theta)*sin(phi);
	  a_1 = sin(theta)*cos(phi);
	  a_2 = cos(theta);
	}
      
      // Beam Angle
      G4double na_0=0, na_1=0, na_2=0;
      if(gunAngle != 0)
	{
	  theta = gunAngle;
	  na_0 = a_0;
	  na_1 =  a_1*cos(theta) + a_2*sin(theta);
	  na_2 = -a_1*sin(theta ) +a_2*cos(theta);
	  a_0 = na_0;
	  a_1 = na_1;
	  a_2 = na_2;
	}
      
      // set parameters
      particleGun->SetParticleEnergy(energy);
      particleGun->SetParticleDefinition(particle);
      particleGun->SetParticlePosition(G4ThreeVector(x0,y0,z0));
      particleGun->SetParticleMomentumDirection(G4ThreeVector(a_0,a_1,a_2));
      particleGun->GeneratePrimaryVertex(anEvent);
      
    }
  else
    {
      G4int evtNb;
      G4double a, b, c, d, e, f, g, h;
      
      evtNb = anEvent->GetEventID(); 

      G4PrimaryVertex* v0;
      G4PrimaryVertex* v1;
      G4ParticleDefinition* particle1;
      G4PrimaryParticle* p0;
      G4PrimaryParticle* p1;

      //Jump over first event and any bad event
      if (!inFile.eof()) inFile >> a >> b >> c >> d >> e >> f >> g >> h;       
      if (evtNb==0)  inFile >> a >> b >> c >> d >> e >> f >> g >> h;         
      while( (abs((int)c)>50) || (abs((int)d)>50) )
         {
         G4cout<<"Found a Bad Primary : Jumping Event"<<G4endl;
	 G4cout << a<<" "<<b<<" "<<c<<" "<< d <<" "<< e <<" "<< f <<" "<< g <<" "<< h << G4endl;
     	 while( (a!=-1) )
	     inFile >> a >> b >> c >> d >> e >> f >> g >> h; 
	 inFile >> a >> b >> c >> d >> e >> f >> g >> h; 
	 G4cout << a<<" "<<b<<" "<<c<<" "<< d <<" "<< e <<" "<< f <<" "<< g <<" "<< h << G4endl;
	 }

      //Now have a good event, let's generate the particles
       v0 = new G4PrimaryVertex(G4ThreeVector(c*mm,d*mm,z0),0.);
       particle1 = G4ParticleTable::GetParticleTable()->FindParticle((G4int)a);
       p0 = new G4PrimaryParticle(particle1);
       p0->SetMomentum(f*MeV, g*MeV, e*MeV);
       v0->SetPrimary(p0);
       anEvent->AddPrimaryVertex(v0);
       
//       while((a!=-1)||(!inFile.eof())) // Does not work, don't know why !
       while((a!=-1))
         {
           inFile >> a >> b >> c >> d >> e >> f >> g >> h; 
	   if( (abs((int)c)<50) && (abs((int)d)<50) && (e>0) ) //Create vertex only if particle is good
	     {
             v1 = new G4PrimaryVertex(G4ThreeVector(c*mm,d*mm,z0),0.);
     	     particle1 = G4ParticleTable::GetParticleTable()->FindParticle((G4int)a);
             p1 = new G4PrimaryParticle(particle1);
             p1->SetMomentum(f*MeV, g*MeV, e*MeV);
             v1->SetPrimary(p1);
             anEvent->AddPrimaryVertex(v1);
	     }
	   else if(a!=-1)
	     {
	     G4cout<<"Found a Bad Secondary"<<G4endl;
	     G4cout << a<<" "<<b<<" "<<c<<" "<< d <<" "<< e <<" "<< f <<" "<< g <<" "<< h << G4endl;
	     }
         }
      }

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

