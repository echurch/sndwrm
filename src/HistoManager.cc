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
/// \file HistoManager.cc
/// \brief Implementation of the HistoManager class
//
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "HistoManager.hh"
#include "G4UnitsTable.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

HistoManager::HistoManager()
  : fFileName("rdecay02")
{
  Book();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

HistoManager::~HistoManager()
{
  delete G4AnalysisManager::Instance();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HistoManager::Book()
{
  // Create or get analysis manager
  // The choice of analysis technology is done via selection of a namespace
  // in HistoManager.hh
  //
  G4AnalysisManager* analysis = G4AnalysisManager::Instance();
  
  analysis->SetFileName(fFileName);
  analysis->SetVerboseLevel(1);
  analysis->SetActivation(true);     //enable inactivation of histos, nTuples
    
  // Default values (to be reset via /analysis/h1/set command)               
  G4int nbins = 100;
  G4double vmin = 0.;
  G4double vmax = 100.;

  // Create all histograms as inactivated 
  // as we have not yet set nbins, vmin, vmax
  //
  ////analysis->SetHistoDirectoryName("histo");  
  ////analysis->SetFirstHistoId(1);
    
  G4int id = analysis->CreateH1("H10","Energy deposit (MeV) in the target",
                       nbins, vmin, vmax);
  analysis->SetH1Activation(id, false);
    
  id = analysis->CreateH1("H11"," total deposits (hits) in the SiPMs",
			  nbins, vmin, 24000*0.8*20.); //12 MeV max
  analysis->SetH1Activation(id, false);

  id = analysis->CreateH1("H12"," total deposits (hits) in the SiPMs - ES",
			  nbins, vmin, 24000*0.8*20.); //12 MeV max
  analysis->SetH1Activation(id, false);

  id = analysis->CreateH1("H13"," total deposits (hits) in the SiPMs - non-ES",
			  nbins, vmin, 24000*0.8*20.); //12 MeV max
  analysis->SetH1Activation(id, false);  

  id = analysis->CreateH1("H21"," total Energy in the SiPMs",
			  nbins, vmin, 0.8*24.); //12 MeV max
  analysis->SetH1Activation(id, false);

  id = analysis->CreateH1("H22"," total Energy in the SiPMs - ES",
			  nbins, vmin, 0.8*24.); //12 MeV max
  analysis->SetH1Activation(id, false);

  id = analysis->CreateH1("H23"," total Energy in the SiPMs - non-ES",
			  nbins, vmin, 0.8*24.); //12 MeV max
  analysis->SetH1Activation(id, false);  
  std::cout << "HistoManager: H1 id: " << id << std::endl;

  const int nbinsx(10);
  const int nbinsy(10);
  const int nbinsz(20);
  double xmin, ymin, zmin, xmax, ymax, zmax;
  xmin = ymin = zmin = 0.0;
  xmax = 3.0;
  ymax = 4.5;
  zmax = 25.0;
  const double maxphotons(1300.); // for photon bombs this is true

  id = analysis->CreateP1("H31"," total deposits (hits) in the SiPMs in xbins",
			  nbinsx, xmin, xmax, 0., maxphotons);
  analysis->SetH1Activation(id, false);
  std::cout << "HistoManager: P1 id: " << id << std::endl;
  id = analysis->CreateP1("H32"," total deposits (hits) in the SiPMs in ybins",
			  nbinsy, ymin, ymax, 0., maxphotons);
  analysis->SetH1Activation(id, false);
  std::cout << "HistoManager: P1 id: " << id << std::endl;
  id = analysis->CreateP1("H33"," total deposits (hits) in the SiPMs in zbins",
			  nbinsz, zmin, zmax, 0., maxphotons);
  analysis->SetH1Activation(id, false);  
  std::cout << "HistoManager: P1 id: " << id << std::endl;

  id = analysis->CreateP2("H34"," total deposits (hits) in the SiPMs in x-y bins",
			  nbinsx, xmin, xmax, nbinsy, ymin, ymax, 0., maxphotons);
  analysis->SetH2Activation(id, false);
  std::cout << "HistoManager: P2 id: " << id << std::endl;
  id = analysis->CreateP2("H35"," total deposits (hits) in the SiPMs in x-z bins",
			  nbinsx, xmin, xmax, nbinsz, zmin, zmax, 0., maxphotons);
  analysis->SetH2Activation(id, false);
  std::cout << "HistoManager: P2 id: " << id << std::endl;
  id = analysis->CreateP2("H36"," total deposits (hits) in the SiPMs in y-z bins",
			  nbinsy, ymin, ymax, nbinsz, zmin, zmax, 0., maxphotons);
  analysis->SetH2Activation(id, false);
  std::cout << "HistoManager: P2 id: " << id << std::endl;

  id = analysis->CreateH1("H14",
                "Anti-coincidence spectrum (MeV) in the traget",
                 nbins, vmin, vmax);
  analysis->SetH1Activation(id, false);

  id = analysis->CreateH1("H15",
                "Anti-coincidence spectrum (MeV) in the detector",
                 nbins, vmin, vmax);
  analysis->SetH1Activation(id, false);  

  id = analysis->CreateH1("H16","Decay emission spectrum (0 - 10 MeV)",
                 nbins, vmin, vmax);
  analysis->SetH1Activation(id, false);  
  
  id = analysis->CreateH1("H17","Decay emission spectrum (0 - 1 MeV)",
                 nbins, vmin, vmax);
  analysis->SetH1Activation(id, false);

  id = analysis->CreateH1("H18","Decay emission spectrum (0 - 0.1 MeV)",
                 nbins, vmin, vmax);
  analysis->SetH1Activation(id, false);
  
  // nTuples
  //
  ////analysis->SetNtupleDirectoryName("ntuple");
  ////analysis->SetFirstNtupleId(1);
  //       
  analysis->CreateNtuple("T1", "Emitted Particles");
  analysis->CreateNtupleDColumn("PID");       //column 0
  analysis->CreateNtupleDColumn("KEnergy");    //column 1
  analysis->CreateNtupleDColumn("Time");      //column 2
  analysis->CreateNtupleDColumn("Weight");    //column 3
  analysis->FinishNtuple();
  
  analysis->CreateNtuple("T2", "RadioIsotopes");
  analysis->CreateNtupleDColumn("PID");       //column 0
  analysis->CreateNtupleDColumn("Time");      //column 1
  analysis->CreateNtupleDColumn("Weight");    //column 2
  analysis->FinishNtuple();
  
  analysis->CreateNtuple("T3", "Energy depositions");
  analysis->CreateNtupleDColumn("EnergyDep");    //column 0
  analysis->CreateNtupleDColumn("Time");      //column 1
  analysis->CreateNtupleDColumn("Weight");    //column 2
  analysis->FinishNtuple();
  
  analysis->CreateNtuple("Tracks", "Track Summaries");
  analysis->CreateNtupleDColumn("PID");       //column 0
  analysis->CreateNtupleDColumn("Z");         //column 1
  analysis->CreateNtupleDColumn("A");         //column 2    
  analysis->CreateNtupleDColumn("KEnergy");    //column 3
  analysis->CreateNtupleDColumn("Time");      //column 4
  analysis->CreateNtupleDColumn("TrkID");    //column 5
  analysis->CreateNtupleDColumn("Startx");         //column 6    
  analysis->CreateNtupleDColumn("Starty");         //column 7    
  analysis->CreateNtupleDColumn("Startz");         //column 8    
  analysis->CreateNtupleDColumn("Px");         //column 9    
  analysis->CreateNtupleDColumn("Py");         //column 10    
  analysis->CreateNtupleDColumn("Pz");         //column 11   
  analysis->CreateNtupleDColumn("Length");         //column 12
  analysis->CreateNtupleDColumn("Event");         //column 13
  analysis->CreateNtupleSColumn("Process");         //column 14
  analysis->FinishNtuple();

  analysis->CreateNtuple("Steps", "Step Summaries"); 
  analysis->CreateNtupleDColumn("Edep");       //column 0
  analysis->CreateNtupleDColumn("Time");         //column 1
  analysis->CreateNtupleDColumn("Weight");         //colum 2
  analysis->CreateNtupleDColumn("X");         //column 3
  analysis->CreateNtupleDColumn("Y");         //column 4
  analysis->CreateNtupleDColumn("Z");         //column 5
  analysis->CreateNtupleDColumn("Length");         //column 6
  analysis->CreateNtupleDColumn("Event");         //column 7
  analysis->CreateNtupleDColumn("PID");         //column 8
  analysis->CreateNtupleDColumn("TrkID");         //column 9
  analysis->CreateNtupleSColumn("Volume");         //column 10
  analysis->CreateNtupleDColumn("CopyNo");         //column 11
  analysis->CreateNtupleSColumn("Material");         //column 12
  analysis->CreateNtupleDColumn("StepNum");         //column 13
  analysis->CreateNtupleSColumn("SProcess");         //column 14
  analysis->CreateNtupleSColumn("TProcess");         //column 15
  analysis->CreateNtupleDColumn("TX");         //column 16
  analysis->CreateNtupleDColumn("TY");         //column 17
  analysis->CreateNtupleDColumn("TZ");         //column 18
  analysis->CreateNtupleSColumn("TVolume");         //column 19
  analysis->FinishNtuple();
  
  analysis->SetNtupleActivation(false);          
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
