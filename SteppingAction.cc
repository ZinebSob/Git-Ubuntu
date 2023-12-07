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
// shall cite the following Geant4-DNA collaboration publication:
// Med. Phys. 37 (2010) 4692-4708
// The Geant4-DNA web site is available at http://geant4-dna.org
//
/// \file SteppingAction.cc
/// \brief Implementation of the SteppingAction class

#include "SteppingAction.hh"
#include "RunAction.hh"
#include "DetectorConstruction.hh"
#include "PrimaryGeneratorAction.hh"

#include "G4EmCalculator.hh"
#include "G4AnalysisManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4SteppingManager.hh"
#include "G4VTouchable.hh"
#include "G4VPhysicalVolume.hh"
#include "CommandLineParser.hh"

#include "G4Electron.hh"
#include "G4Proton.hh"
#include "G4Gamma.hh"
#include "G4Alpha.hh"
#include "G4DNAGenericIonsManager.hh"
#include "G4EventManager.hh"
#include "G4Event.hh"

#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

using namespace G4DNAPARSER;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::SteppingAction() : G4UserSteppingAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::~SteppingAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//Modification
std::ofstream  file1 ("resultat.dat");	
  G4int NumberOfStep;
  //G4double i;
  G4int CompteurIoni;		// Compteur des ionisations secondaires
  G4int CompteurExci;		// Compteur des excitations secondaires
  G4int CompteurSolvation;	// Compteur électron solvation
  G4int CompteurAtt;		// Compteur électron attachement
  G4int CompteurVibExci;	// Compteur Vibration Excitation
  G4int CompteurElastic;	// Compteur des processus élastiques pour un électron
  //G4int CompteurCaptureElectronique;	// Compteur de capture électronique 
  //G4int CompteurMSC;			// Compteur Multiple Scattering
  //G4int CompteureIonisationStandard;	// Compteur Ionisation Standard
  G4int CompteurpElastic;		// Compteur des processus elastiques du proton
  G4int CompteurpExci;			// Compteur d'excitation par proton
  G4int CompteurpIoni;			// Compteur d'ionisation par proton
  G4int CompteurpCapture;		// Compteur de capture électronique
  G4int CompteurpIonisationStandard;   // Compteur d'ionisation standard par proton
  G4int Compteurpmsc;
  G4int CompteurCoulombScat;
  G4int CompteurnuclearStopping;
  G4int CompteurhElastic;		// Compteur des processus élastiques après transmutation du proton en atome d'hydrogène
  G4int CompteurhExci;			// Compteur d'excitation par hydrogène
  G4int CompteurhIoni;  		// Compteur d'ionisation par hydrogène
  G4int CompteurhPerte1electron;	// Compteur de perte d'un électron
  G4double SumOfStepLength;	// Longueur total sommé sur le nombre de step	
  G4double SumOfTrackLength;	
  G4double SumMsteplength; 	// Somme des longueurs des steps 
  G4double SumOfEnergy1; 	// Somme d'énergie déposée pour chaque step
  G4double SumOfEnergy2;	// Somme d'énergie perdue pour chaque step
  //G4double Range1;
  G4double StoppingPower1;	// Pouvoir d'arret = Energie déposée totale /Somme des longueurs des steps
  G4double StoppingPower2;	// Pouvoir d'arret = Energie perdue totale /Somme Moyenne des longueurs des steps
  G4double StoppingPower3;	// Pouvoir d'arret de chaque step
  G4double StoppingPower4;	// Sommation du pouvoir d'arret de chaque step
  G4double StoppingPower5;	// Pouvoir d'arret de chaque step
  G4double StoppingPower6;	// Sommation du pouvoir d'arret de chaque step
  G4double MeanFreePath;	// Libre Parcours Moyen
  G4double steplength ;	// Longueur de Step instantanée
  G4double tracklength ;
  G4double Msteplength ;	// Moyen de Longueur de Step 
  G4double edepstep ;          // Energie déposée pour chaque step
  G4double destep ;  		// Energie perdue pour chaque step
  G4double Ikinetic ;  	// Energie incidente pour chaque step
  G4double Ipkinetic ;		// Energie post step 
  G4double deltaEnergy ;
  G4ThreeVector deltaPosition;
  G4ThreeVector deltaMomentum ; 
  G4double SumOfEnergy ;
  G4double ProjectedRange ;
  G4double SumProjectedRange ;
  G4double SP ;
  G4double SPIntrinseque;
  G4double Range;
  G4double R;
  //G4double Range ;
  //G4double SP ;

  
void SteppingAction::UserSteppingAction(const G4Step* step)
{ 

/* Instanciate EmCalculator
 G4EmCalculator emCal;
 Get material
 const G4Material* material = fDetector->GetMaterial();
  G4NistManager * pMan = G4NistManager::Instance();
  const G4Material *pWater = pMan->FindOrBuildMaterial("G4_WATER");
*/ 
  
  G4double flagParticle=-1;
  G4double flagProcess=-1;
  G4int TID = 0;			// Track ID
  G4int PID = 0;			// Parent ID
  G4int CSN = 0;		        // Current Step Number
  G4double x,y,z;			// PreStep Positions
  G4double xp,yp,zp;			// PostStep Positions
  G4int EventIn = 0;			// Evenement Initial
  G4int EventID = 0;			// Nouveau Evenement

  // Particle identification

  // The following method avoids the usage of string comparison 

  G4ParticleDefinition * partDef = 
   step->GetTrack()->GetDynamicParticle()->GetDefinition();

  if (partDef == G4Gamma::GammaDefinition())
     flagParticle = 0; 

  if (partDef == G4Electron::ElectronDefinition())
     flagParticle = 1; 

  if (partDef == G4Proton::ProtonDefinition())
     flagParticle = 2; 

  if (partDef == G4Alpha::AlphaDefinition())
     flagParticle = 4; 

  G4DNAGenericIonsManager *instance;
  instance = G4DNAGenericIonsManager::Instance();
  
  // Usage example
  /*
  G4ParticleDefinition* protonDef = G4Proton::ProtonDefinition();
  G4ParticleDefinition* hydrogenDef = instance->GetIon("hydrogen");
  G4ParticleDefinition* alphaPlusPlusDef = instance->GetIon("alpha++");
  G4ParticleDefinition* alphaPlusDef = instance->GetIon("alpha+");
  G4ParticleDefinition* heliumDef = instance->GetIon("helium");
  */

  if (partDef == instance->GetIon("hydrogen"))
     flagParticle = 3; 

  if (partDef == instance->GetIon("alpha+"))
     flagParticle = 5; 

  if (partDef == instance->GetIon("helium"))
     flagParticle = 6; 
     
  // Alternative method (based on string comparison)
  
  /*
  const G4String& particleName = step->GetTrack()->GetDynamicParticle()->
      GetDefinition()->GetParticleName();
      
  if (particleName == "gamma")         flagParticle = 0;
  else if (particleName == "e-")       flagParticle = 1;
  else if (particleName == "proton")   flagParticle = 2;
  else if (particleName == "hydrogen") flagParticle = 3;
  else if (particleName == "alpha")    flagParticle = 4;
  else if (particleName == "alpha+")   flagParticle = 5;
  else if (particleName == "helium")   flagParticle = 6;
  */

  // Process identification

  // Process sub-types are listed in G4PhysicsListHelper.cc
  // or in Geant4-DNA process class implementation files (*.cc)

  G4StepPoint * preStep = step->GetPreStepPoint();
  G4StepPoint * postStep = step->GetPostStepPoint();
 // G4int procID = postStep->GetProcessDefinedStep()->GetProcessSubType();

  const G4String& processName = postStep->
      GetProcessDefinedStep()->GetProcessName();
 
  if (processName=="eCapture") flagProcess =1;
  // (no subType and procID exists at the moment for this process)
/*
  else if (flagParticle == 1)

  {
    if      (procID==58) flagProcess =10;
    else if (procID==51) flagProcess =11;
    else if (procID==52) flagProcess =12;
    else if (procID==53) flagProcess =13;
    else if (procID==55) flagProcess =14;
    else if (procID==54) flagProcess =15;
    else if (procID==10) flagProcess =110;
    else if (procID==1)  flagProcess =120;
    else if (procID==2)  flagProcess =130;
  }
        
  else if (flagParticle == 2)

  {
    if      (procID==51) flagProcess =21;
    else if (procID==52) flagProcess =22;
    else if (procID==53) flagProcess =23;
    else if (procID==56) flagProcess =24;
    else if (procID==10) flagProcess =210;
    else if (procID==1)  flagProcess =220;
    else if (procID==2)  flagProcess =230;
    else if (procID==8)  flagProcess =240;
  }
        
  else if (flagParticle == 3)

  {
    if      (procID==51) flagProcess =31;
    else if (procID==52) flagProcess =32;
    else if (procID==53) flagProcess =33;
    else if (procID==57) flagProcess =35;
  }
        
  else if (flagParticle == 4)

  {
    if      (procID==51) flagProcess =41;
    else if (procID==52) flagProcess =42;
    else if (procID==53) flagProcess =43;
    else if (procID==56) flagProcess =44;
    else if (procID==10) flagProcess =410;
    else if (procID==1)  flagProcess =420;
    else if (procID==2)  flagProcess =430;
    else if (procID==8)  flagProcess =440;
  }
        
  else if (flagParticle == 5)

  {
    if      (procID==51) flagProcess =51;
    else if (procID==52) flagProcess =52;
    else if (procID==53) flagProcess =53;
    else if (procID==56) flagProcess =54;
    else if (procID==57) flagProcess =55;
    else if (procID==10) flagProcess =510;
    else if (procID==1)  flagProcess =520;
    else if (procID==2)  flagProcess =530;
    else if (procID==8)  flagProcess =540;
  }
        
  else if (flagParticle == 6)

  {
    if      (procID==51) flagProcess =61;
    else if (procID==52) flagProcess =62;
    else if (procID==53) flagProcess =63;
    else if (procID==57) flagProcess =65;
  }

  else if (processName=="GenericIon_G4DNAIonisation") flagProcess =73;
  else if (processName=="msc")                        flagProcess =710;
  else if (processName=="CoulombScat")                flagProcess =720;
  else if (processName=="ionIoni")                    flagProcess =730;
  else if (processName=="nuclearStopping")            flagProcess =740;
  // (for all GenericIons)
  */
  // Alternatively, using process names

  
  if (processName=="e-_G4DNAElectronSolvation")         flagProcess =10;
  else if (processName=="e-_G4DNAElastic")              flagProcess =11;
  else if (processName=="e-_G4DNAExcitation")           flagProcess =12;
  else if (processName=="e-_G4DNAIonisation")           flagProcess =13;
  else if (processName=="e-_G4DNAAttachment")           flagProcess =14;
  else if (processName=="e-_G4DNAVibExcitation")        flagProcess =15;
  else if (processName=="eCapture")                     flagProcess =16;
  else if (processName=="msc")                          flagProcess =17;
  else if (processName=="eIoni")                        flagProcess =130;

  else if (processName=="proton_G4DNAElastic")          flagProcess =21;
  else if (processName=="proton_G4DNAExcitation")       flagProcess =22;
  else if (processName=="proton_G4DNAIonisation")       flagProcess =23;
  else if (processName=="proton_G4DNAChargeDecrease")   flagProcess =24;
  else if (processName=="hIoni")                        flagProcess =230;

  else if (processName=="hydrogen_G4DNAElastic")        flagProcess =31;
  else if (processName=="hydrogen_G4DNAExcitation")     flagProcess =32;
  else if (processName=="hydrogen_G4DNAIonisation")     flagProcess =33;
  else if (processName=="hydrogen_G4DNAChargeIncrease") flagProcess =35;
  
  else if (processName=="GenericIon_G4DNAIonisation") flagProcess =73;
  else if (processName=="msc")                        flagProcess =710;
  else if (processName=="CoulombScat")                flagProcess =720;
  else if (processName=="ionIoni")                    flagProcess =730;
  else if (processName=="nuclearStopping")            flagProcess =740;
  
/*
  else if (processName=="alpha_G4DNAElastic")           flagProcess =41;
  else if (processName=="alpha_G4DNAExcitation")        flagProcess =42;
  else if (processName=="alpha_G4DNAIonisation")        flagProcess =43;
  else if (processName=="alpha_G4DNAChargeDecrease")    flagProcess =44;

  else if (processName=="alpha+_G4DNAElastic")          flagProcess =51;
  else if (processName=="alpha+_G4DNAExcitation")       flagProcess =52;
  else if (processName=="alpha+_G4DNAIonisation")       flagProcess =53;
  else if (processName=="alpha+_G4DNAChargeDecrease")   flagProcess =54;
  else if (processName=="alpha+_G4DNAChargeIncrease")   flagProcess =55;

  else if (processName=="helium_G4DNAElastic")          flagProcess =61;
  else if (processName=="helium_G4DNAExcitation")       flagProcess =62;
  else if (processName=="helium_G4DNAIonisation")       flagProcess =63;
  else if (processName=="helium_G4DNAChargeIncrease")   flagProcess =65;

  else if (processName=="GenericIon_G4DNAIonisation")   flagProcess =73;
*/
  
  if (processName!="Transportation")
{  
  
// Initialization for each new event 
 EventID =  G4EventManager::GetEventManager()->
                         GetConstCurrentEvent()->GetEventID();
 
  if (EventID != EventIn )
   	{
   		EventIn = EventID;
   		SumOfStepLength = 0.;
   		SumMsteplength = 0.;
   		x = 0;
   		y = 0;
   		z = 0;
   		xp = 0;
   		yp = 0;
   		zp = 0;
   		steplength = 0.;
   		Msteplength = 0;
   		edepstep = 0;
   		destep = 0;
   		SumOfEnergy1 = 0;
   		SumOfEnergy2 = 0.;
   		Ikinetic = 0;
   		Ipkinetic =0.; 
   		TID = 0;
   		PID = 0;
   		CSN = 0;
   		//psquare = 0;
   		//Range1 =0.;
               StoppingPower1 = 0.;
               StoppingPower2 = 0.;   
		StoppingPower3 = 0.;
		StoppingPower4 = 0.; 
		StoppingPower5 = 0.; 
		StoppingPower6 = 0.; 
		deltaEnergy = 0.;
    		SumOfEnergy =0.;
    		ProjectedRange =0.;
    		SumProjectedRange =0.;
    		SP = 0.;
    		SPIntrinseque = 0.;
    		//Range = 0.;
    		//SP = 0.;	
	} 
	// End Initialization for each new event 
	
	
	 if (flagParticle == 2){

// G4cout<< " -----------> Initialisation Event : "<<EventIn<< G4endl;

    //NumberOfStep ++;  
    x=preStep->GetPosition().x()/nanometer;
    y=preStep->GetPosition().y()/nanometer;
    z=preStep->GetPosition().z()/nanometer;

    xp=postStep->GetPosition().x()/nanometer;
    yp=postStep->GetPosition().y()/nanometer;
    zp=postStep->GetPosition().z()/nanometer;
    /*
    Mx=preStep->GetMomentumDirection().x();
    My=preStep->GetMomentumDirection().y();
    Mz=preStep->GetMomentumDirection().z();
    
    Mxp=postStep->GetMomentumDirection().x();
    Myp=postStep->GetMomentumDirection().y();
    Mzp=postStep->GetMomentumDirection().z();
    */
    PID = step->GetTrack()->GetParentID();
    TID = step->GetTrack()->GetTrackID();
    CSN = step->GetTrack()->GetCurrentStepNumber();

 // Les variables utilisées en opérations 
 
    Ikinetic = preStep->GetKineticEnergy();
    Ipkinetic = postStep->GetKineticEnergy();
    destep =  std::fabs(Ikinetic - Ipkinetic);
    deltaEnergy = step->GetDeltaEnergy();
   // destep = std::fabs(deltaEnergy)
   // psquare = preStep->GetMomentumDirection()*postStep->GetMomentumDirection();
   // MeanFreePath = (SumOfStepLength / NumberOfStep);
   // Range1 = emCal.GetCSDARange(1.*MeV,G4Proton::Proton(),pWater);
   
    Msteplength = std::sqrt((x-xp)*(x-xp) + (y-yp)*(y-yp) + (z-zp)*(z-zp));
    steplength = step->GetStepLength();
    tracklength = step->GetTrack()->GetTrackLength();
    edepstep = step->GetTotalEnergyDeposit(); 
    
    
  //if(std::fabs(deltaEnergy) > 0.*CLHEP::eV)					// (8.22*eV) Seuil d'énergie pour calculer les cassures simples et doubles brins
  
  //Les opérations pour chaque step
  if(edepstep!=0)	
    {   
    	NumberOfStep ++;
    	//deltaEnergy = step->GetDeltaEnergy();
    	deltaPosition = step->GetDeltaPosition();
    	deltaMomentum = step->GetDeltaMomentum();
    	ProjectedRange = std::fabs(deltaPosition * deltaMomentum);
    	SumProjectedRange += ProjectedRange;
        SumOfStepLength += steplength;
        SumOfTrackLength += tracklength;
        Range += SumOfStepLength;
        R=Range/1000;
        SumMsteplength += Msteplength;
        SumOfEnergy1 += edepstep;
        SumOfEnergy2 += destep;
        SumOfEnergy += deltaEnergy;					 // Sommation des énergies perdues avant et après chaque step
        //Range += deltaPosition;					 // Sommation des différences des positions avant et après chaque step
        SP = -(SumOfEnergy/SumOfStepLength);				 // Pouvoir d'arret
        SPIntrinseque = -(deltaEnergy/steplength);			 // Pouvoir d'arret intrinseque
    	StoppingPower1 = (edepstep/steplength);			 // Pouvoir d'arret total instantané
    	StoppingPower2 += StoppingPower1;				 // Somme du pouvoir d'arret instantané
    	StoppingPower3 = (SumOfEnergy1/SumOfStepLength); 		 // Pour vérifier Pouvoir d'arret = somme d'énergie perdue/somme des longueurs des steps
    	//SP1 = (StoppingPower2/NumberOfStep);
    	StoppingPower4 = (destep/ Msteplength);
    	StoppingPower5 += StoppingPower4;
    	StoppingPower6 = (SumOfEnergy2/SumOfStepLength); 		// Vérification
    	//SP2 = (StoppingPower6/NumberOfStep); 
	    
 if (processName=="e-_G4DNAElectronSolvation")        CompteurSolvation++;
  else if (processName=="e-_G4DNAElastic")            CompteurElastic++;
  else if (processName=="e-_G4DNAExcitation")         CompteurExci++;
  else if (processName=="e-_G4DNAIonisation")         CompteurIoni++;
  else if (processName=="e-_G4DNAAttachment")         CompteurAtt++;
  else if (processName=="e-_G4DNAVibExcitation")      CompteurVibExci++;
  
  else if (processName=="proton_G4DNAElastic")        CompteurpElastic++;
  else if (processName=="proton_G4DNAExcitation")     CompteurpExci++;
  else if (processName=="proton_G4DNAIonisation")     CompteurpIoni++;
  else if (processName=="proton_G4DNAChargeDecrease") CompteurpCapture++;
  else if (processName=="hIoni")                      CompteurpIonisationStandard++;
  else if (processName=="msc")                        Compteurpmsc++;
  else if (processName=="CoulombScat")                CompteurCoulombScat++;
  else if (processName=="nuclearStopping")            CompteurnuclearStopping++;

  else if (processName=="hydrogen_G4DNAElastic")       CompteurhElastic++;
  else if (processName=="hydrogen_G4DNAExcitation")    CompteurhExci++;
  else if (processName=="hydrogen_G4DNAIonisation")    CompteurhIoni++;
  else if (processName=="hydrogen_G4DNAChargeIncrease") CompteurhPerte1electron++;
 
 file1 << "     " << SumOfStepLength/cm << "     " << SumOfEnergy1/MeV << "     " << StoppingPower3/(MeV/cm)<< G4endl;
  
 /*file1 << "     " << CompteurSolvation << "     " << CompteurElastic << "     " << CompteurExci <<"     " << CompteurIoni << "     " << CompteurAtt << "     " << CompteurVibExci << "     " << CompteurpElastic << "     " << CompteurpExci << "     " << CompteurpIoni << "     " << CompteurpCapture << "     " << CompteurpIonisationStandard << "     " << Compteurpmsc << "     " << CompteurCoulombScat << "     " << CompteurnuclearStopping << "     " << CompteurhElastic << "     " << CompteurhExci << "     " << CompteurhIoni << "     " << CompteurhPerte1electron << "     " << SumOfStepLength/cm << G4endl;*/
  
      }
 
  /*
  else if (processName=="eCapture")                   CompteurCaptureElectronique;
  else if (processName=="msc")                        CompteurMSC;
  else if (processName=="eIoni")                      CompteureIonisationStandard;
*/

   
// get analysis manager
    CommandLineParser* parser = CommandLineParser::GetParser();
    Command* command(0);
    if((command = parser->GetCommandIfActive("-out"))==0) return;

    G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();

    // fill ntuple
    analysisManager->FillNtupleDColumn(0, flagParticle);
    analysisManager->FillNtupleDColumn(1, flagProcess);
    analysisManager->FillNtupleDColumn(2, x);
    analysisManager->FillNtupleDColumn(3, y);
    analysisManager->FillNtupleDColumn(4, z);
    analysisManager->FillNtupleDColumn(5, SumOfEnergy);
    analysisManager->FillNtupleDColumn(6, steplength/micrometer);
    analysisManager->FillNtupleDColumn(7, destep/eV );
    analysisManager->FillNtupleDColumn(8, Ikinetic/eV);
    //analysisManager->FillNtupleDColumn(9, psquare);        
    analysisManager->FillNtupleIColumn(9, EventID);
    analysisManager->FillNtupleIColumn(10, TID);
    analysisManager->FillNtupleIColumn(11, PID);
    analysisManager->FillNtupleIColumn(12, CSN);                                   
    analysisManager->FillNtupleDColumn(13, SumOfStepLength/micrometer);
    analysisManager->FillNtupleDColumn(14, SumMsteplength);
    analysisManager->FillNtupleDColumn(15, SP);
    analysisManager->FillNtupleDColumn(16, StoppingPower2);
    analysisManager->FillNtupleDColumn(17, StoppingPower3);
    analysisManager->FillNtupleDColumn(18, StoppingPower4);
    analysisManager->FillNtupleDColumn(19, StoppingPower5);
    analysisManager->FillNtupleDColumn(20, StoppingPower6);
    //analysisManager->FillNtupleDColumn(20, SigmaMFP);
    //analysisManager->FillNtupleDColumn(21, Sigma);
    //analysisManager->FillNtupleDColumn(22, meanfreepath);
    //analysisManager->FillNtupleDColumn(23, Summeanfreepath);
    analysisManager->AddNtupleRow();
    //}
    edepstep=0.; // Initialisation pour chaque step
    steplength=0.; // Initialisation pour chaque step
    tracklength=0.;
}
   }
 
  //Output Data
  //G4cout << "############################################################################################"  << G4endl;
  //G4cout << "  Event number " << EventID << " ---------> Range1 = " << G4BestUnit( Range1,"Length") << G4endl;
  /*
  G4cout << "  Event number " << EventID << " ---------> steplength = " << G4BestUnit( steplength,"Length") << G4endl;
  G4cout << "  Event number " << EventID << " ---------> edepstep = " << G4BestUnit( edepstep,"Energy") << G4endl;
  G4cout << "  Event number " << EventID << " ---------> destep = " << G4BestUnit( destep,"Energy") << G4endl;
  G4cout << "  Event number " << EventID << " ---------> SumOfStepLength = " << G4BestUnit( SumOfStepLength,"Length") << G4endl;
  //G4cout << "  Event number " << EventID << " ---------> Msteplength = " << G4BestUnit( Msteplength,"Length") << G4endl;
  //G4cout << "  Event number " << EventID << " ---------> SumMStepLength = " << G4BestUnit( SumMStepLength,"Length") << G4endl;
  //G4cout << "  Event number " << EventID << " ---------> Ikinetic  = " << G4BestUnit( Ikinetic,"Energy") << G4endl;
  //G4cout << "  Event number " << EventID << " ---------> Ipkinetic  = " << G4BestUnit( Ipkinetic,"Energy") << G4endl;
  //G4cout << "  Event number " << EventID << " ---------> Deposited Energy  = " << G4BestUnit( edepstep,"Energy") << G4endl;
  //G4cout << "  Event number " << EventID << " ---------> Difference Kinetic Energy  = " << G4BestUnit( destep,"Energy") << G4endl;
  G4cout << "  Event number " << EventID << " ---------> Sum of Deposited Energy  = " << G4BestUnit( SumOfEnergy1,"Energy") << G4endl;
  //G4cout << "  Event number " << EventID << " ---------> Sum of Difference Kinetic Energy  = " << G4BestUnit( SumOfEnergy2,"Energy") << G4endl;
  G4cout << "  Event number " << EventID << " ---------> StoppingPower1 = " << G4BestUnit(StoppingPower1,"Energy/Length") << G4endl;
  G4cout << "  Event number " << EventID << " ---------> StoppingPower2 = " << G4BestUnit(StoppingPower2,"Energy/Length") << G4endl;
  G4cout << "  Event number " << EventID << " ---------> StoppingPower3 = " << G4BestUnit(StoppingPower3,"Energy/Length") << G4endl;
  G4cout << "  Event number " << EventID << " ---------> StoppingPower4 = " << G4BestUnit(StoppingPower4,"Energy/Length") << G4endl;
  G4cout << "  Event number " << EventID << " ---------> StoppingPower5 = " << G4BestUnit(StoppingPower5,"Energy/Length") << G4endl;
  G4cout << "  Event number " << EventID << " ---------> StoppingPower6 = " << G4BestUnit(StoppingPower6,"Energy/Length") << G4endl;
  G4cout << "  Event number " << EventID << " ---------> NumberOfStep = " << NumberOfStep << G4endl;
 //G4cout << "  Event number " << EventID << " ---------> CSN = " << CSN << G4endl;
 //G4cout << "  Event number " << EventID << " ---------> CompteurIoni = " << CompteurIoni << G4endl;
 // G4cout << "  Event number " << EventID << " ---------> CompteurExci = " << CompteurExci << G4endl;
 // G4cout << "  Event number " << EventID << " ---------> CompteurSolvation = " << CompteurSolvation << G4endl;
 // G4cout << "  Event number " << EventID << " ---------> CompteurAtt = " << CompteurAtt << G4endl;
 // G4cout << "  Event number " << EventID << " ---------> CompteurVibExci = " << CompteurVibExci << G4endl;
 // G4cout << "  Event number " << EventID << " ---------> Range1 = " << Range1 << G4endl;
 // G4cout << "  Event number " << EventID << " ---------> MeanFreePath = " << G4BestUnit(MeanFreePath,"Length") << G4endl;
 // G4cout << "  Event number " << EventID << " ---------> AverageMFP = " << G4BestUnit(averagemfp,"Length") << G4endl;
 // G4cout << "  Event number " << EventID << " ---------> SigmaMFP = " << G4BestUnit(SigmaMFP,"Length") << G4endl;
 // G4cout << "  Event number " << EventID << " ---------> Sigma = " << G4BestUnit(Sigma,"Surface") << G4endl;
 // G4cout << "  Event number " << EventID << " ---------> Summeanfreepath = " << G4BestUnit(Summeanfreepath,"Length") << G4endl;
  
 */
    }


  
	
