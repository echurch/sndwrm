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
/// \file SteppingAction.cc
/// \brief Implementation of the SteppingAction class
//
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "SteppingAction.hh"

#include "DetectorConstruction.hh"
#include "Run.hh"
#include "EventAction.hh"
#include "HistoManager.hh"

#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"

#include "G4OpticalPhoton.hh" 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::SteppingAction(DetectorConstruction* det, EventAction* event)
: G4UserSteppingAction(), fDetector(det), fEventAction(event)
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::~SteppingAction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SteppingAction::UserSteppingAction(const G4Step* aStep)
{
  Run* run = static_cast<Run*>(
        G4RunManager::GetRunManager()->GetNonConstCurrentRun());    
  G4int  event = G4EventManager::GetEventManager()->GetConstCurrentEvent()->GetEventID();

  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();

  //which volume ?
  //
  G4LogicalVolume* lVolume = aStep->GetPreStepPoint()->GetTouchableHandle()
                             ->GetVolume()->GetLogicalVolume();
  G4int iVol = 0;
  G4Track* track = aStep->GetTrack();
  G4double tID = track->GetTrackID();



//  if (lVolume == fDetector->GetLogicTarget())   iVol = 1;
 // if (lVolume == fDetector->GetLogicDetector()) iVol = 2;


  // count processes
  // 
  const G4StepPoint* endPoint = aStep->GetPostStepPoint();
  const G4VProcess* tprocess   = endPoint->GetProcessDefinedStep();
  const G4VProcess* sprocess   = aStep->GetPreStepPoint()->GetProcessDefinedStep();
  run->CountProcesses(tprocess, iVol);

  const G4ParticleDefinition* particle = track->GetParticleDefinition();  
  G4int pID       = particle->GetPDGEncoding();
  G4TouchableHandle touch = endPoint->GetTouchableHandle();
  G4VPhysicalVolume* eVolume = touch->GetVolume();
  G4String eVname("null");
  if (eVolume)
    {
      eVname = eVolume->GetName();  
      //            std::cout << "SteppingAction eVname: " << eVname << std::endl;
/*      if (lVolume == fDetector->GetLogicSiPM() || eVolume->GetLogicalVolume() == fDetector->GetLogicSiPM() || 
	  eVname.find("SiPM")!=std::string::npos || (lVolume->GetName()).find("SiPM")!=std::string::npos)     iVol = 3;
*/
      if (( eVname.find("Arapuca")!=std::string::npos || (lVolume->GetName()).find("Arapuca")!=std::string::npos ) and (pID==0 || pID==-22))
	//	if ( (eVname.find("Arapuca")!=std::string::npos)  and (pID==0 || pID==-22))
	{
	  iVol = 4;
	  //	  std::cout << "SteppingAction: Optical photon hit an Arapuca: " << eVname << std::endl;
	}
    }


  
  // energy deposit
  //
  G4double edepStep = aStep->GetTotalEnergyDeposit();

  tprocess = aStep->GetPostStepPoint()->GetProcessDefinedStep();




  if (edepStep <= 0. && (pID!=0 && pID!=-22) ) return; // the deception version of G4 uses -22 for optical photons; my Mac's uses 0.

  G4double time   = aStep->GetPreStepPoint()->GetGlobalTime();
  G4double weight = aStep->GetPreStepPoint()->GetWeight();   

  //fill ntuple id = 2
  G4int id = 4;   
  const G4double length = aStep->GetStepLength();
  const G4ThreeVector pos(aStep->GetPreStepPoint()->GetPosition());
  const G4ThreeVector tpos(aStep->GetPostStepPoint()->GetPosition());

  std::string startp("null");
  std::string endp("null");


  /*
  const std::vector<double> fidv {
      GetEvtAct()->GetPrimGenAct()->GetParticleGun()->GetCurrentSource()->GetPosDist()->GetHalfX(),
      GetEvtAct()->GetPrimGenAct()->GetParticleGun()->GetCurrentSource()->GetPosDist()->GetHalfY(),
      GetEvtAct()->GetPrimGenAct()->GetParticleGun()->GetCurrentSource()->GetPosDist()->GetHalfZ()
      } ;
  */

  // Above does not work when launching n's, gammas from outside the fidV. G'arr! Hard code it for now. EC, 4-Aug-2021.

  const std::vector<double> fidv {3000,4500,20000};
  if (sprocess)
      startp = sprocess->GetProcessName();

  // If an optical photon is born outside fiducialvolume let's kill it. EC, 4-Aug-2021.
  // Idea being that we'd reco this vtx outside our fidv and cut the event.
  /*
  if ((pID == 0 or pID == -22 ) and
      ( ( abs(pos[0]) > fidv.at(0) ) or ( abs(pos[1]) > fidv.at(1) ) or ( abs(pos[2]) > fidv.at(2) ) ) 
      //      and track->GetCurrentStepNumber() <= 1  // seems ok looking at steps, but concerned it's biasing. EC, 5-Aug-2021.
      and track->GetCurrentStepNumber() == 0 and iVol!=4 // don't do this for Arapuca hit counting mode
      )
    {
      //      std::cout << "SteppingAction(): Killing optical photon Track at pos " << pos[0] << ", " << pos[1]  << ", " << pos[2] << std::endl;
      track->SetTrackStatus(fStopAndKill);
      return;
    }
  */

  if (iVol!=3)
    fEventAction->AddEdep(iVol, edepStep, time, weight);

  if (eVname=="SiPM") {
	  fEventAction->AddEdep(3, 1.0, time, weight);	  
  }

  //  if (iVol==4)
    //    std::cout << "Hit Arapuca, evolume/copyNo, particle: " << eVname << "/" << eVolume->GetCopyNo() << ", " << pID << std::endl;
//  if (iVol!=4 and fDetector->GetAPEX()) return; // do not fill the steps TTree if we haven't stepped into an Arapuca
  
  const G4ThreeVector tspos(track->GetVertexPosition()); // Let's grab the opt photon's point of origin, not the step's, which may not be the same thing. EC, 16-Aug-2023.
  analysisManager->FillNtupleDColumn(id,0, edepStep);
  analysisManager->FillNtupleDColumn(id,1, time/s);
  analysisManager->FillNtupleDColumn(id,2, weight);
  analysisManager->FillNtupleDColumn(id,3, tspos[0]/mm);
  analysisManager->FillNtupleDColumn(id,4, tspos[1]/mm);
  analysisManager->FillNtupleDColumn(id,5, tspos[2]/mm);
  analysisManager->FillNtupleDColumn(id,6, length/mm);
  analysisManager->FillNtupleDColumn(id,7, event);
  analysisManager->FillNtupleDColumn(id,8, pID);
  analysisManager->FillNtupleDColumn(id,9, tID);

  analysisManager->FillNtupleSColumn(id,10, track->GetVolume()->GetLogicalVolume()->GetName());
  if (eVolume)
    analysisManager->FillNtupleDColumn(id,11, eVolume->GetCopyNo());
  else
    analysisManager->FillNtupleDColumn(id,11, 0);
  analysisManager->FillNtupleSColumn(id,12, track->GetVolume()->GetLogicalVolume()->GetMaterial()->GetName());
  analysisManager->FillNtupleDColumn(id,13, track->GetCurrentStepNumber());
  if (sprocess)
      startp = sprocess->GetProcessName();
  if (tprocess)
    endp = tprocess->GetProcessName();
  analysisManager->FillNtupleSColumn(id,14, startp);
  analysisManager->FillNtupleSColumn(id,15, endp);
  analysisManager->FillNtupleDColumn(id,16, tpos[0]/mm);
  analysisManager->FillNtupleDColumn(id,17, tpos[1]/mm);
  analysisManager->FillNtupleDColumn(id,18, tpos[2]/mm);
  analysisManager->FillNtupleSColumn(id,19,eVname);

  analysisManager->AddNtupleRow(id);      



}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
