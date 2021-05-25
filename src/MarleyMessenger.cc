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
/// \file DetectorMessenger.cc
/// \brief Implementation of the DetectorMessenger class
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "MarleyMessenger.hh"


#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithADouble.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

MarleyMessenger::MarleyMessenger(PrimaryGeneratorAction* PGA)
:G4UImessenger(), 
 fPrimaryGeneratorAction(PGA),
 fMarleyConfFileCmd(0)
{ 
         
  fRdecayDir = new G4UIdirectory("/rdecay02/");
  fRdecayDir->SetGuidance("commands specific to this example");

  G4bool broadcast = false;
  fMarleyDir = new G4UIdirectory("/rdecay02/marley/",broadcast);
  fMarleyDir->SetGuidance("detector construction commands");

  fMarleyConfFileCmd = new G4UIcmdWithAString("/rdecay02/marley/setConfigFile",this);
  fMarleyConfFileCmd->SetGuidance("Select Marley ConfigFile");
  fMarleyConfFileCmd->SetParameterName("choice",false);
  fMarleyConfFileCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

MarleyMessenger::~MarleyMessenger()
{
  delete fMarleyConfFileCmd;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void MarleyMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{ 
  if (command == fMarleyConfFileCmd )
   { fPrimaryGeneratorAction->SetMarleyConfFile(newValue);
     std::cout << "MarleyMessenger::SetNewValue() config file: " << fPrimaryGeneratorAction->GetMarleyConfFile() << std::endl;

     fPrimaryGeneratorAction->ReloadMarleyConfig();
   }

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
