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
  fFiducial = true;

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
  if (std::fabs(time - fTime0) > TimeWindow) return;
  
  if (iVol == 1) { fEdep1 += edep; fWeight1 += edep*weight;}
  if (iVol == 3 and GetFiducial()) { fEdep2 += edep; fWeight2 += edep*weight;}  

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

   std::cout << "EndofEvtAction: Total SipM hits: " << fEdep2 << std::endl;

   analysisManager->FillH1(4, Energy, 1.0);    //fWeight2);
   if (fPGA->GetPrimaryGenerator()->GetFSNeutrino())
     analysisManager->FillH1(5, Energy, 1.0);
   else 
     analysisManager->FillH1(6, Energy, 1.0);
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
  
 // Below takes advantage of the x8 symmetry of our detector
 for (int ii=0; ii<int(bins->size()/3-1) ;ii++)
   {
     if (!(abs(vtx[0]/1000.) >= bins->at(ii) && abs(vtx[0]/1000.)< bins->at(ii+1))) continue;
     for (int jj=bins->size()/3; jj<int(2*bins->size()/3-1) ;jj++)
       {
	 if (!(abs(vtx[1]/1000.) >= bins->at(jj) && abs(vtx[1]/1000.)<bins->at(jj+1))) continue;
	 for (int kk=2*bins->size()/3; kk<int(3*bins->size()/3-1) ;kk++)
	   {
	     if (!(abs(vtx[2]/1000.) >= bins->at(kk) && abs(vtx[2]/1000.)<bins->at(kk+1))) continue;
	     std::cout << "EnergyCalc using bins " << ii <<"," << jj << "," << kk << std::endl;
	     int dataind (ii + jj + kk);
	     energyConv = data->at(dataind);
	     std::cout << "EnergyCalc returning ph/MeV,MeV " << energyConv << ", " << E/energyConv << std::endl;
	     break;
	   }
       }
   }


  return E/energyConv;

}
