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
/// \file TrackingAction.cc
/// \brief Implementation of the TrackingAction class
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "TrackingAction.hh"

#include "DetectorConstruction.hh"
#include "Run.hh"
#include "HistoManager.hh"

#include "G4RunManager.hh"
#include "G4Track.hh"
#include "G4HadronicProcessType.hh"

#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TrackingAction::TrackingAction(DetectorConstruction* det, EventAction* evt)
 :G4UserTrackingAction(), fDetector(det), fEventAction(evt)
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TrackingAction::PreUserTrackingAction(const G4Track* track)
{
  
  Run* run = static_cast<Run*>(
        G4RunManager::GetRunManager()->GetNonConstCurrentRun());    
  
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();

  //which volume ?
  //
  G4LogicalVolume* lVolume = track->GetVolume()->GetLogicalVolume();
  G4int iVol = 0;
//  if (lVolume == fDetector->GetLogicTarget())   iVol = 1;
 // if (lVolume == fDetector->GetLogicDetector()) iVol = 2;
    
  //secondary particles only. EC wants all.
  //  if (track->GetTrackID() == 1) return;
  G4int event;
  event = G4EventManager::GetEventManager()->GetConstCurrentEvent()->GetEventID();
  const G4ParticleDefinition* particle = track->GetParticleDefinition();  
  G4String name   = particle->GetParticleName();
  G4int pid       = particle->GetPDGEncoding();
  G4int Z         = particle->GetAtomicNumber();
  G4int A         = particle->GetAtomicMass();
  G4double charge = particle->GetPDGCharge();    
  G4double energy = track->GetKineticEnergy();
  G4double time   = track->GetGlobalTime();
  G4double tID = track->GetTrackID();
  G4double weight = track->GetWeight();
  const G4ThreeVector& vtx = track->GetVertexPosition();
  const G4double length = track->GetTrackLength();
  const G4ThreeVector pvtx =  track->GetMomentumDirection () ;
  
  const G4VProcess* process   = track->GetCreatorProcess();
  G4String processName("null") ;
  if (process)
    {
      processName = process->GetProcessName();
    }  
  run->ParticleCount(name,energy,iVol);
  
  std::vector<std::string> procOfInterest({"capt","beta","radioactive"}); // to catch ncapt, Rn222 chain, Ar39, Ar42 chain betas. And primaries with "null".
  // Add PairProduction and ComptonScattering to get the gammas I launch from Bi214,Tl208 from outside. EC, 5-Aug-2021.
  procOfInterest.push_back("pair");
  procOfInterest.push_back("comp");
  std::for_each(processName.begin(), processName.end(), [](char & c) {
      c = ::tolower(c);
    });



  /*
  const std::vector<double> fidv {
      GetEvtAct()->GetPrimGenAct()->GetParticleGun()->GetCurrentSource()->GetPosDist()->GetHalfX(),
      GetEvtAct()->GetPrimGenAct()->GetParticleGun()->GetCurrentSource()->GetPosDist()->GetHalfY(),
      GetEvtAct()->GetPrimGenAct()->GetParticleGun()->GetCurrentSource()->GetPosDist()->GetHalfZ()
      } ;
  */

  // no worky when launching n's, gammas from outside the fidv. Hard code it. EC, 4-Aug-2021.
  const std::vector<double> fidv {3000,4500,20000};

  if (abs(vtx[0])<fidv.at(0) && abs(vtx[1])<fidv.at(1) && abs(vtx[2])<fidv.at(2) && !fEventAction->GetFiducial())
    {
      for (const auto& proc : procOfInterest) {
	if (processName.find(proc) != std::string::npos) {
	  fEventAction->SetProcVtx(vtx);
	  //  	  std::cout << "TrackingAction: Interesting Process is " << proc << std::endl;
	  fEventAction->SetFiducial(true);
	  break;
	}
      }
    }

  //  std::cout << "TrackingAction():  Track at pos " << vtx[0] << ", " << vtx[1]  << ", " << vtx[2] << std::endl;
  //Radioactive decay products
  //G4int procaessType = track->GetCreatorProcess()->GetProcessSubType();
  //  if (processType == fRadioactiveDecay) {
    //fill ntuple id = 3

  //    std::cout << "TrackingAction::PreUserTrackingAction()..." << std::endl ;
  // std::cout << "\t pid, energy, length, " << pid << ", " << energy << ", "  << length << std::endl;
    G4int id = 3;
    analysisManager->FillNtupleDColumn(id,0, double(pid));
    analysisManager->FillNtupleDColumn(id,1, double(Z));
    analysisManager->FillNtupleDColumn(id,2, double(A));
    analysisManager->FillNtupleDColumn(id,3, energy);
    analysisManager->FillNtupleDColumn(id,4, time/s);
    analysisManager->FillNtupleDColumn(id,5, tID);
    analysisManager->FillNtupleDColumn(id,6, vtx[0]);
    analysisManager->FillNtupleDColumn(id,7, vtx[1]);
    analysisManager->FillNtupleDColumn(id,8, vtx[2]);
    analysisManager->FillNtupleDColumn(id,9,  pvtx[0]/sqrt(pvtx[0]*pvtx[0]+pvtx[1]*pvtx[1]+pvtx[2]*pvtx[2]));
    analysisManager->FillNtupleDColumn(id,10, pvtx[1]/sqrt(pvtx[0]*pvtx[0]+pvtx[1]*pvtx[1]+pvtx[2]*pvtx[2]));
    analysisManager->FillNtupleDColumn(id,11, pvtx[2]/sqrt(pvtx[0]*pvtx[0]+pvtx[1]*pvtx[1]+pvtx[2]*pvtx[2]));
    analysisManager->FillNtupleDColumn(id,12, length);
    analysisManager->FillNtupleDColumn(id,13, event);
    analysisManager->FillNtupleSColumn(id,14, processName);
    analysisManager->AddNtupleRow(id);

    if (tID == 1) {   // primaries
      //fill ntuple id = 0
      id = 0;
      analysisManager->FillNtupleDColumn(id,0, double(pid));
      analysisManager->FillNtupleDColumn(id,1, energy);
      analysisManager->FillNtupleDColumn(id,2, time/s);
      analysisManager->FillNtupleDColumn(id,3, event);
      analysisManager->AddNtupleRow(id);
    
      analysisManager->FillH1(7, energy, weight);
      analysisManager->FillH1(8, energy, weight);
    } 

    //  }
  
  //all unstable ions produced in target
  G4bool unstableIon = ((charge > 2.) && !(particle->GetPDGStable()));
  if ((unstableIon) && (iVol == 1)) {
    //fill ntuple id = 1
    id = 1;
    analysisManager->FillNtupleDColumn(id,0, double(pid));
    analysisManager->FillNtupleDColumn(id,1, time/s);
    analysisManager->FillNtupleDColumn(id,2, weight);
    analysisManager->AddNtupleRow(id);  
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TrackingAction::PostUserTrackingAction(const G4Track* )
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

