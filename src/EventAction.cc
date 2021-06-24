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
/// \file EventAction.cc
/// \brief Implementation of the EventAction class
//
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "EventAction.hh"


#include "HistoManager.hh"

#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"

#include "G4LogicalVolumeStore.hh"
#include "G4Material.hh"
#include "G4MaterialPropertiesTable.hh"
#include "G4MaterialPropertyVector.hh"
#include "G4ThreeVector.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::EventAction()
:G4UserEventAction(),
 fEdep1(0.), fEdep2(0.), fWeight1(0.), fWeight2(0.),
 fTime0(-1*s)
{ } 

EventAction::EventAction(PrimaryGeneratorAction* prim)
 :G4UserEventAction(), fPGA(prim), 
 fEdep1(0.), fEdep2(0.), fWeight1(0.), fWeight2(0.),
 fTime0(-1*s)
{ } 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::~EventAction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::BeginOfEventAction(const G4Event*)
{
  fEdep1 = fEdep2 = fWeight1 = fWeight2 = 0.;
  fTime0 = -1*s;
  fFiducial = false;
  procVtx[0] = procVtx[1] = procVtx[2] = 0.0;

  G4LogicalVolumeStore * lvs =   G4LogicalVolumeStore::GetInstance();

  for ( G4LogicalVolumeStore::iterator i = lvs->begin(); i != lvs->end(); ++i ){

    G4LogicalVolume* volume = (*i);
    G4Material* TheMaterial = volume->GetMaterial();
    std::string Material = TheMaterial->GetName();
    std::string volLower = volume->GetName();
    G4MaterialPropertiesTable* aMaterialPropertiesTable = TheMaterial->GetMaterialPropertiesTable();

    if (aMaterialPropertiesTable)
      {
	/*
	G4MaterialPropertyVector* Fast_Intensity =
	  aMaterialPropertiesTable->GetProperty(kFASTCOMPONENT);
	G4MaterialPropertyVector* Slow_Intensity =
	  aMaterialPropertiesTable->GetProperty(kSLOWCOMPONENT);

	std::cout <<  Material << " Fast Intensity" << std::endl;
	Fast_Intensity->DumpValues();
	std::cout << Material << " Slow Intensity" << std::endl;
	Slow_Intensity->DumpValues();
	*/
      }
  }


}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::AddEdep(G4int iVol, G4double edep,
                                      G4double time, G4double weight)
{
  // initialize t0
  if (fTime0 < 0.) fTime0 = time;
  
  // out of time window ?
  const G4double TimeWindow (1*microsecond);
  //  if (std::fabs(time - fTime0) > TimeWindow) return;
  
  G4ThreeVector vtx(fPGA->GetPrimaryGenerator()->GetParticlePosition());
  if (iVol == 1) { fEdep1 += edep; fWeight1 += edep*weight;}

  // Use GetFiducial for the physics processes that follow in the particle/radioactive decay samples 
  // If we ran Marley (vtx!=0) just assume they were run in the fiducial volume. EC, 6-May-2021.
  if (iVol == 3 and (GetFiducial() || vtx[0]!=0.0)) { fEdep2 += edep; fWeight2 += edep*weight; }

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::EndOfEventAction(const G4Event*)
{
 G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
 
 // pulse height in target
 //
 if (fEdep1 > 0.) {
   fWeight1 /= fEdep1;

   //   std::cout << "fEdep1:: " << fEdep1  << std::endl;
   analysisManager->FillH1(0, fEdep1, fWeight1);   
 }
 
 // pulse height in SiPM
 //   
 if (fEdep2 > 0.) {
   fWeight2 /= fEdep2;
   std::cout << "EndofEvtAction: Total SipM hits: " << fEdep2 << std::endl;
   analysisManager->FillH1(1, fEdep2, 1.0);    //fWeight2);
   if (fPGA->GetPrimaryGenerator()->GetFSNeutrino())
     analysisManager->FillH1(2, fEdep2, 1.0);
   else 
     analysisManager->FillH1(3, fEdep2, 1.0);

 }
   

 // Now get PhotonsToMeV data and fill calibrated-by-Primary-Launch-Location Histograms of Energy.
 Run* run = static_cast<Run*>(
            G4RunManager::GetRunManager()->GetNonConstCurrentRun());

 if (fEdep2 > 0.) {
   // Fill 3 histograms according after converting to MeV
   double Energy = EnergyCalc(fEdep2,run);

   std::cout << " EventAction::EndOfEventAction(): Total SipM hits: " << fEdep2 << std::endl;

   analysisManager->FillH1(4, Energy, 1.0);    //fWeight2);
   if (fPGA->GetPrimaryGenerator()->GetFSNeutrino())
     analysisManager->FillH1(5, Energy, 1.0);
   else 
     analysisManager->FillH1(6, Energy, 1.0);

   G4ThreeVector vtx(fPGA->GetPrimaryGenerator()->GetParticlePosition());
   // std::cout << "EventAction::EndOfEventAction(): primaryVtx: " << vtx[0]/m <<"," << vtx[1]/m << "," << vtx[2]/m  << std::endl;
   analysisManager->FillP1(0, fabs(vtx[0]/m) , fEdep2); 
   analysisManager->FillP1(1, fabs(vtx[1]/m) , fEdep2); 
   analysisManager->FillP1(2, fabs(vtx[2]/m) , fEdep2); 
   analysisManager->FillP2(0, fabs(vtx[0]/m), fabs(vtx[1]/m) , fEdep2); 
   analysisManager->FillP2(1, fabs(vtx[0]/m), fabs(vtx[2]/m) , fEdep2); 
   analysisManager->FillP2(2, fabs(vtx[1]/m), fabs(vtx[2]/m) , fEdep2); 

 }



 // run->AddEdep (fEdep1, fEdep2);             

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double EventAction::EnergyCalc(G4double E, Run* run)
{

 std::vector<double> *data( run->GetPhotonsToMeVData());
 std::vector<double> *bins( run->GetPhotonsToMeVBins());
 G4double energyConv(0.);

 G4ThreeVector vtx(fPGA->GetPrimaryGenerator()->GetParticlePosition());
 std::cout << "EnergyCalc(): vtx is " << vtx[0] << "," << vtx[1] << "," << vtx[2] << std::endl;

 // If vtx=0, it's because this event was initiatiated by PrimaryGeneratorAction, not PrimaryGenerator. That means
 // we want to use the vtx of the interesting  process we found and saved in TrackingAction::PreUserTrackingAction().
 if (vtx[0]==0. && vtx[1]==0. && vtx[2]==0.) {
   vtx = GetProcVtx();
   std::cout << "EnergyCalc(): Interesting fiducial process's vtx is " << vtx[0] << "," << vtx[1] << "," << vtx[2] << std::endl;
 }

 int nbinedges(bins->size()/3);  
 int nbins = nbinedges-1;

 // Below takes advantage of the x8 symmetry of our detector
 for (int ii=0; ii<int(nbinedges*3/nbins-1) ;ii++)
   {
     if (!(abs(vtx[0]/1000.) >= bins->at(ii) && abs(vtx[0]/1000.)< bins->at(ii+1))) continue;
     for (int jj=nbinedges*3/nbins; jj<int(2*nbinedges*3/nbins-1) ;jj++)
       {
	 if (!(abs(vtx[1]/1000.) >= bins->at(jj) && abs(vtx[1]/1000.)<bins->at(jj+1))) continue;
	 for (int kk=2*nbinedges*3/nbins; kk<int(3*nbinedges*3/nbins-1) ;kk++)
	   {
	     if (!(abs(vtx[2]/1000.) >= bins->at(kk) && abs(vtx[2]/1000.)<bins->at(kk+1))) continue;
	     std::cout << "EnergyCalc using nbins " << ii <<"," << jj << "," << kk << std::endl;
	     int dataind (ii*9 + (jj-nbinedges)*3 + (kk-nbinedges*2)*1); 
	     energyConv = data->at(dataind);
	     std::cout << "EnergyCalc returning ph/MeV,MeV " << energyConv << ", " << E/energyConv << std::endl;
	     break;
	   }
       }
   }


  return E/energyConv;

}

