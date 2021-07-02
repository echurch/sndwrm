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
/// \file PrimaryGeneratorAction.cc
/// \brief Implementation of the PrimaryGeneratorAction class
//
//
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "PrimaryGeneratorAction.hh"

#include "G4Event.hh"
#include "G4ParticleTable.hh"
#include "G4IonTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4Geantino.hh"
#include "G4NeutrinoE.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::PrimaryGeneratorAction()
  : G4VUserPrimaryGeneratorAction(),fParticleGun(0),
    fPrimaryGenerator(0), config_file_name(std::string(std::getenv("MARLEY"))+"/examples/config/annotated.js")
{
  //  G4int n_particle = 1;
  //  fParticleGun  = new G4ParticleGun(n_particle);
  fParticleGun  = new G4GeneralParticleSource();

  fPrimaryGenerator = new PrimaryGenerator();  
  /*
  fParticleGun->SetParticleEnergy(0*eV);
  fParticleGun->SetParticlePosition(G4ThreeVector(0.*cm,0.*cm,0.*cm));
  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,1.));
  */

  // now if user specifies a Marley config file in mac file, use that instead of annotated.js.
  fMarleyMessenger = new MarleyMessenger(this); 
  marley::JSONConfig config(config_file_name);
  std::cout << "PrimaryGenAction: Marley config file is " << config_file_name << std::endl;
  std::cout << "PrimaryGenAction: Unused if nu_e is not specified "  << std::endl;
  marley_generator_= config.create_generator();

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void PrimaryGeneratorAction::ReloadMarleyConfig()
{
  std::cout << "PrimaryGeneratorAction:ReloadMarleyConfig(): Loading new config " << config_file_name << std::endl;

  marley::JSONConfig config(config_file_name);
  std::cout << "PrimaryGenAction: Marley config file is " << config_file_name << std::endl;
  std::cout << "PrimaryGenAction: Unused if nu_e is not specified "  << std::endl;
  marley_generator_= config.create_generator();

}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
  if (fPrimaryGenerator)
    delete fPrimaryGenerator;
  if (fParticleGun)
    delete fParticleGun;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  // if there's a geantino in .mac file create a "photon bomb"
  std::vector<double> xyzbounds { 
      fParticleGun->GetCurrentSource()->GetPosDist()->GetHalfX(),
      fParticleGun->GetCurrentSource()->GetPosDist()->GetHalfY(),
      fParticleGun->GetCurrentSource()->GetPosDist()->GetHalfZ()};

  std::cout << "GeneratePrimaries() particle requested is " << fParticleGun->GetParticleDefinition() << std::endl;

  if (fParticleGun->GetParticleDefinition() == G4Geantino::Geantino()) {  
    std::cout << "GeneratePrimaries: detect that a 'photon bomb' is to be created." << std::endl;

    fPrimaryGenerator->GeneratePrimaryVertexOpt(anEvent,xyzbounds);
    return;
  }    
  //create vertex
  //   

  else if (fParticleGun->GetParticleDefinition() == G4NeutrinoE::NeutrinoE()) {
    std::cout << "GeneratePrimaries: detect that a Marley event is to be created." << std::endl;

    // Generate a new MARLEY event using the owned marley::Generator object  
    marley::Event ev = marley_generator_.create_event();
    fPrimaryGenerator->GeneratePrimaryVertexMarley(anEvent,xyzbounds,ev);

    return;
  }

  fParticleGun->GeneratePrimaryVertex(anEvent);


  // fPrimaryGenerator->GeneratePrimaryVertexExtra(event);  
  // see extended/eventgenerator/userPrimaryGeneratorSource/src/PrimaryGeneratorAction.cc::GeneratePrimaryVertex(), which should work to add other 
  // particles beyond those generated from our macro to the event. EC, 6-Feb-2021 

  //  ---------- OR  --------------------

  // do it perhaps from the mac file itself! 
  // see test31.g4mac at  http://hurel.hanyang.ac.kr/Geant4/Geant4_GPS/reat.space.qinetiq.com/gps/examples/examples.html

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


