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
/// \file rdecay02.cc
/// \brief Main program of the radioactivedecay/rdecay02 example
//
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4Types.hh"

#include "G4RunManagerFactory.hh"
#ifdef G4MULTITHREADED
#include "G4MTRunManager.hh"
#else
#include "G4RunManager.hh"
#endif

#include "G4UImanager.hh"
#include "Randomize.hh"

#include "DetectorConstruction.hh"
#include "PhysicsList.hh"
#include "ActionInitialization.hh"
#include "SteppingVerbose.hh"

#include "Shielding.hh"
#include "G4ThermalNeutrons.hh"

#include "G4UIExecutive.hh"
#include "G4VisExecutive.hh"

#include "G4OpticalPhysics.hh"
#include "G4EmStandardPhysics_SpacePhysics.hh"
#include "G4HadronInelasticQBBC.hh"
#include "XrayTESdetPhysicsList.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int main(int argc,char** argv) {

  std::cout << "Run interactively as ./sndwrm ../macros/mymacro.mac FloorShieldThicknessIntegerIncm Nthr." << std::endl;
  std::cout << "Can leave off last two arguments and use defaults." << std::endl;
  
  //detect interactive mode (if no arguments) and define UI session
  G4UIExecutive* ui = 0;
  if (argc == 1) ui = new G4UIExecutive(argc,argv);

  //choose the Random engine
  G4Random::setTheEngine(new CLHEP::RanecuEngine);

  // Construct the default run manager
#if defined(G4MULTITHREADED) && defined(USEG4MT)
  //G4MTRunManager* runManager = new G4MTRunManager;
  auto runManager = G4RunManagerFactory::CreateRunManager(G4RunManagerType::Default); // This, as opposed to above line, does not cause a million destructor complaints. EC, 6-May-2025.
  G4int nThreads = G4Threading::G4GetNumberOfCores();
  if (argc==4) nThreads = G4UIcommand::ConvertToInt(argv[3]);
  runManager->SetNumberOfThreads(nThreads);
  std::cout << "sndwrm: RUNNING IN G4MULTITHREADED MODE." << std::endl;
#else
  //my Verbose output class
  G4VSteppingVerbose::SetInstance(new SteppingVerbose);
  // To fix big dump-out error about leaks upon exit. EC, 30-Dec-2024.
  //G4RunManager* runManager = new G4RunManager
  std::cout << "sndwrm: RUNNING in single-threaded mode." << std::endl;
  auto* runManager =
    G4RunManagerFactory::CreateRunManager(G4RunManagerType::SerialOnly);  // SerialOnly);  //Default);
#endif

  //set mandatory initialization classes
  DetectorConstruction* det= new DetectorConstruction;
  runManager->SetUserInitialization(det);

  /*  
  PhysicsList* phys = new PhysicsList;
  runManager->SetUserInitialization(phys);
  */

  
  // EC, 30-Apr-2024. Replace longstanding use of crafting my own physics list. ... tacking on args that should enforce LIQMD_HPT, 23-June-2025  
  //  G4VModularPhysicsList* physlist = new Shielding(1,"HP","",true); // 1 for verbose.
  // EC, 20-Dec-2025. SpacePhysics
  auto *physlist = new XrayTESdetPhysicsList;
  runManager->SetUserInitialization(physlist);

  G4int verb(0);
  G4OpticalPhysics* opticalPhysics = new G4OpticalPhysics(verb);
  /* Comment out all below for successful compilation in G4 11, EC, 20-Dec-2024.
  opticalPhysics->SetWLSTimeProfile("delta");
  opticalPhysics->SetMaxBetaChangePerStep(10.0);
  opticalPhysics->SetTrackSecondariesFirst(kCerenkov,true);
  opticalPhysics->SetTrackSecondariesFirst(kScintillation,true);
  opticalPhysics->SetScintillationYieldFactor(1.);
  G4int fMaxNumPhotonStep(7000);
  opticalPhysics->SetMaxNumPhotonsPerStep(fMaxNumPhotonStep);
  */
  // Don't want OPs for space physics
  //physlist->RegisterPhysics( opticalPhysics);
  
  runManager->SetUserInitialization(new ActionInitialization(det));

  //initialize visualization
  G4VisManager* visManager = nullptr;

  // get the pointer to the User Interface manager
  G4UImanager* UImanager = G4UImanager::GetUIpointer();

  if (ui)  {
   //interactive mode
   visManager = new G4VisExecutive;
   visManager->Initialize();
   ui->SessionStart();
   delete ui;
  }
  else  {
   //batch mode
   G4String command = "/control/execute ";
   G4String fileName = argv[1];
   if (argc>=3)
     {
       std::stringstream alias;
       alias<<"ARG"<< 1 <<" "<< argv[2];	
       UImanager->SetAlias(alias.str().c_str());
     }
   UImanager->ApplyCommand(command+fileName);
  }

  std::cout << "sndwrm: About to call destructors, which might complain. Not calling 'em would dump a __ton__ of complaints." << std::endl;
  //job termination
  if (visManager)
    delete visManager;
  if (runManager)
    delete runManager;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
