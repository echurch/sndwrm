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

#include "MaterialMessenger.hh"

#include "MaterialPropertyLoader.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithADouble.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

MaterialMessenger::MaterialMessenger(MaterialPropertyLoader* MPL)
:G4UImessenger(), 
 fMaterialPropertyLoader(MPL),
 fMaterialG10SpecRefCmd(0),
 fMaterialAcrylicSpecRefCmd(0),
 fMaterialLArAbsLengthCmd(0)
{ 
         
  fRdecayDir = new G4UIdirectory("/rdecay02/");
  fRdecayDir->SetGuidance("commands specific to this example");

  std::cout << "MatMess:SetNewVal. ctor is " << std::endl;
  G4bool broadcast = false;
  fMaterialDir = new G4UIdirectory("/rdecay02/material/",broadcast);
  fMaterialDir->SetGuidance("material commands");


  fMaterialG10SpecRefCmd = new G4UIcmdWithADouble("/rdecay02/material/setG10SpecRef",this);
  fMaterialG10SpecRefCmd->SetGuidance("Set G10 Spectral Reflection");
  fMaterialG10SpecRefCmd->SetParameterName("choice",false);
  fMaterialG10SpecRefCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  fMaterialAcrylicSpecRefCmd = new G4UIcmdWithADouble("/rdecay02/material/setAcrylicSpecRef",this);
  fMaterialAcrylicSpecRefCmd->SetGuidance("Set Acrylic Spectral Reflection");
  fMaterialAcrylicSpecRefCmd->SetParameterName("choice",false);
  fMaterialAcrylicSpecRefCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  fMaterialLArAbsLengthCmd = new G4UIcmdWithADoubleAndUnit("/rdecay02/material/setLArAbsLength",this);
  fMaterialLArAbsLengthCmd->SetGuidance("Set LAr Absorption Length");
  fMaterialLArAbsLengthCmd->SetUnitCategory("Length");
  fMaterialLArAbsLengthCmd->SetParameterName("choice",false);
  fMaterialLArAbsLengthCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

MaterialMessenger::~MaterialMessenger()
{
  delete fMaterialG10SpecRefCmd;
  delete fMaterialAcrylicSpecRefCmd;
  delete fMaterialLArAbsLengthCmd;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void MaterialMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{ 

  if (command == fMaterialG10SpecRefCmd )
    { fMaterialPropertyLoader->SetMaterialG10SpecRef(fMaterialG10SpecRefCmd->GetNewDoubleValue(newValue));
     std::cout << "MaterialMessenger::SetNewValue() config file: " << fMaterialPropertyLoader->GetMaterialG10SpecRef() << std::endl;
   }
  else if (command == fMaterialAcrylicSpecRefCmd )
    { fMaterialPropertyLoader->SetMaterialAcrylicSpecRef(fMaterialAcrylicSpecRefCmd->GetNewDoubleValue(newValue));
     std::cout << "MaterialMessenger::SetNewValue() config file: " << fMaterialPropertyLoader->GetMaterialAcrylicSpecRef() << std::endl;
   }
  else if (command == fMaterialLArAbsLengthCmd )
    { fMaterialPropertyLoader->SetMaterialLArAbsLength(fMaterialLArAbsLengthCmd->GetNewDoubleValue(newValue));
     std::cout << "MaterialMessenger::SetNewValue() config file: " << fMaterialPropertyLoader->GetMaterialLArAbsLength() << std::endl;
   }

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
