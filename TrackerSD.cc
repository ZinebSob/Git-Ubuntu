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
// This example is provided by the Geant4-DNA collaboration
// Any report or published results obtained using the Geant4-DNA software 
// shall cite the following Geant4-DNA collaboration publications:
// Phys. Med. 31 (2015) 861-874
// Med. Phys. 37 (2010) 4692-4708
// The Geant4-DNA web site is available at http://geant4-dna.org
//
/// \file TrackerSD.cc
/// \brief Implementation of the TrackerSD class

#include "TrackerSD.hh"
#include "Randomize.hh"
#include "G4AnalysisManager.hh"
#include "G4SDManager.hh"

#include "G4RandomDirection.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TrackerSD::TrackerSD(const G4String& name,
                     const G4String& hitsCollectionName) 
:G4VSensitiveDetector(name),
fHitsCollection(NULL)
{
  collectionName.insert(hitsCollectionName);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TrackerSD::~TrackerSD() 
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TrackerSD::Initialize(G4HCofThisEvent* hce)
{
  // Create hits collection
  fHitsCollection 
    = new TrackerHitsCollection(SensitiveDetectorName, collectionName[0]); 

  // Add this collection in hce

  G4int hcID 
    = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]);
  
  hce->AddHitsCollection( hcID, fHitsCollection ); 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool TrackerSD::ProcessHits(G4Step* aStep, 
                                     G4TouchableHistory*)
{  
  // energy deposit
  G4double edep = aStep->GetTotalEnergyDeposit();
  // Step Length
  G4double steplength = aStep->GetStepLength();

  if (edep==0.) return false;

  TrackerHit* newHit = new TrackerHit();

  newHit->SetTrackID  (aStep->GetTrack()->GetTrackID());
  newHit->SetEdep(edep);
  newHit->SetPos (aStep->GetPostStepPoint()->GetPosition());
  newHit->SetLength(steplength);

  if (aStep->GetTrack()->GetTrackID()==1&&aStep->GetTrack()->GetParentID()==0)
    newHit->SetIncidentEnergy(aStep->GetTrack()->GetVertexKineticEnergy());

  fHitsCollection->insert( newHit );

  //newHit->Print();

  return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TrackerSD::EndOfEvent(G4HCofThisEvent*)
{
  G4int nofHits = fHitsCollection->entries();
//  G4int eventID = aStep->GetTrack()->GetTrackID();	
  G4double Einc=0;
  G4double xpos=0;
  G4double ypos=0;
  G4double zpos=0;
  G4double edepose =0.;
  //G4double vedepose =0.;
  G4double range = 0.;
  G4double length = 0.;
  G4double spower = 0.;
  /*
  G4cout << G4endl
  << "-------->Hits Collection: in this event they are " 
  << nofHits 
  << " hits in the target volume " << G4endl;
  */
  
  // PROCESSING OF MICRODOSIMETRY Y & Z SPECTRA
  
  // *************************************
  // Please select herebelow :
  // the radius of the target sphere:
  // variable name = radius
  // it is set to 5 nm by default)
  //
  
  G4double radius = 1.*nm;

  //
 
  //***************
  // y and z
  //***************
  
  // select random hit
  G4int randHit=0; // Runs from 0 to number of hits - 1
  randHit = static_cast<G4int>( G4UniformRand()*nofHits );
  
   
  /*
  G4cout 
  << "======> random selection of hit number randHit =" 
  << randHit << G4endl;
  */
  
  // get selected random hit position
  G4ThreeVector hitPos =  (*fHitsCollection)[randHit]->GetPos();
//G4cout << "======> random hit position x/nm =" << hitPos.x()/nm << G4endl; 
//G4cout << "======> random hit position y/nm =" << hitPos.y()/nm << G4endl; 
//G4cout << "======> random hit position z/nm =" << hitPos.z()/nm << G4endl; 
  
  // set random position of center of sphere within radius
/*  G4double chord = 4.*radius/3;
  G4double density = 1 * g/cm3;
  G4double mass = (4./3)*CLHEP::pi*radius*radius*radius*density;*/
  
  // random placement of sphere: method 1
  /*  
  G4ThreeVector randDir = G4RandomDirection();
  G4double randRadius = G4UniformRand()*radius;
  G4ThreeVector randCenterPos = randRadius*randDir + hitPos;
  */  

  // random placement of sphere: method 2

  G4double xRand = 1.01*radius;
  G4double yRand = 1.01*radius;
  G4double zRand = 1.01*radius;
  G4double randRad = 1.01*radius;
  do
  {
    xRand = (2*G4UniformRand()-1)*radius;
    yRand = (2*G4UniformRand()-1)*radius;
    zRand = (2*G4UniformRand()-1)*radius;
    randRad = std::sqrt( xRand*xRand+yRand*yRand+zRand*zRand );
  }
  while (randRad>radius);

  G4ThreeVector 
    randCenterPos(xRand+hitPos.x(),yRand+hitPos.y(),zRand+hitPos.z());
    
      // added by zineb to take random set position of center of sphere within radius
  randRad = radius;
  G4double chord = 4.*randRad/3;
  G4double density = 1 * g/cm3;
  G4double mass = (4./3)*CLHEP::pi*randRad*randRad*randRad*density;
     // end by zineb
  // search for neighbouring hits in the sphere and cumulate deposited energy 
  //  in epsilon
  G4double epsilon = 0;
  G4int nbEdep = 0;
  G4double LocalRange =0.;
  G4double LocalSP =0.;
  
  for ( G4int i=0; i<nofHits; i++ ) 
  { 
	G4double Cedepose = (*fHitsCollection)[i]->GetEdep();
	if(Cedepose !=0.)
	 { 
    if ((*fHitsCollection)[i]->GetIncidentEnergy()>0)
      Einc = (*fHitsCollection)[i]->GetIncidentEnergy();
      
   	G4ThreeVector localPos = (*fHitsCollection)[i]->GetPos();
         xpos =  localPos.x()/nm;
         ypos =  localPos.y()/nm;
         zpos =  localPos.z()/nm;
    
    	
    	if(((*fHitsCollection)[i]->GetLength()) !=0.){ 
    	edepose += (*fHitsCollection)[i]->GetEdep();
    	length = (*fHitsCollection)[i]->GetLength();
  	range += length;
  	spower = edepose/range; 
  	}
  	
    
    if ( 
        (localPos.x()-randCenterPos.x()) * (localPos.x()-randCenterPos.x()) +
        (localPos.y()-randCenterPos.y()) * (localPos.y()-randCenterPos.y()) +
        (localPos.z()-randCenterPos.z()) * (localPos.z()-randCenterPos.z()) 
         <= radius*radius
       ) 
       
    { 
      epsilon += (*fHitsCollection)[i]->GetEdep() ;
      nbEdep = nbEdep+1;
      LocalRange += (*fHitsCollection)[i]->GetLength();
      LocalSP = epsilon/LocalRange;
    }
      }  
  }

  // for testing only
  /*
  G4cout << "======> for hit number #" << randHit <<
  ", we collect " 
  << nbEdep << " energy depositions in a sphere of radius " 
  << radius/nm << " nm and mass " 
  << mass/kg << " kg for a total of " 
  << epsilon/eV << " eV or " 
  << (epsilon/joule)/(mass/kg) << " Gy" << G4endl;
  G4cout << "-" << G4endl;
  
  
  
  FILE* myFile;
  myFile=fopen("yz.txt","a");
  fprintf(myFile,"%e %e %e %e %e \n",radius/nm,(epsilon/eV)/(chord/nm),(epsilon/joule)/(mass/kg),spower/(eV/nm),range/micrometer);
  fclose(myFile);
 */

  // get analysis manager
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();

  // fill ntuple including weighting
 // analysisManager->FillNtupleDColumn(0, eventID);
  analysisManager->FillNtupleDColumn(1, randRad/nm);
  analysisManager->FillNtupleDColumn(2, xpos/nm);
  analysisManager->FillNtupleDColumn(3, ypos/nm);
  analysisManager->FillNtupleDColumn(4, zpos/nm);
  analysisManager->FillNtupleDColumn(5, length/nm);
  analysisManager->FillNtupleDColumn(6, range/nm);
  analysisManager->FillNtupleDColumn(7, nofHits);
  analysisManager->FillNtupleDColumn(8, nbEdep);
  analysisManager->FillNtupleDColumn(9, (epsilon/eV)/(chord/nm));
  analysisManager->FillNtupleDColumn(10, (epsilon/mass)/gray);
  analysisManager->FillNtupleDColumn(11, edepose/eV);
  analysisManager->FillNtupleDColumn(12, LocalSP/(eV/nm));
  analysisManager->FillNtupleDColumn(13, spower/(eV/nm));
  analysisManager->FillNtupleDColumn(14, Einc/eV);
  analysisManager->AddNtupleRow();

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
