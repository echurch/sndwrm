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

#include "DetectorMessenger.hh"

#include "DetectorConstruction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWith3VectorAndUnit.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithADouble.hh"



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorMessenger::DetectorMessenger(DetectorConstruction* Det)
:G4UImessenger(), 
 fDetector(Det), fRdecayDir(0), fDetDir(0),
 fTargMatCmd(0), fDetectMatCmd(0), fTargRadiusCmd(0), fInsetRadiusCmd(0),
 fDetectThicknessCmd(0), fTargLengthCmd(0), fDetectLengthCmd(0) 
{ 


  fRdecayDir = new G4UIdirectory("/rdecay02/");
  fRdecayDir->SetGuidance("commands specific to this example");
  
  G4bool broadcast = false;
  fDetDir = new G4UIdirectory("/rdecay02/det/",broadcast);
  fDetDir->SetGuidance("detector construction commands");
  fDetDirsw = new G4UIdirectory("/sndwrm/det/",broadcast);
  fDetDirsw->SetGuidance("detector construction commands");

  fFidVolumeCmd = new G4UIcmdWith3VectorAndUnit("/sndwrm/det/setHalfFidVol", this);
  fFidVolumeCmd->SetGuidance("Set the Fiducial Volume.");
  fFidVolumeCmd->SetParameterName("FidVolx","FidVoly","FidVolz",true);
  fFidVolumeCmd->SetDefaultUnit("m");
  fFidVolumeCmd->SetDefaultValue(fFidVolumeCmd->GetNew3VectorValue("6. 6. 30. m"));
  fFidVolumeCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  fDetector->SetFidVolume(fFidVolumeCmd->GetNew3VectorValue("6. 6. 30. m"));
  std::cout << "DetMess::pre-SetNewVal() FidVol " << fDetector->GetFidVolume() << std::endl;

  fFloorShieldCmd = new G4UIcmdWithADoubleAndUnit("/sndwrm/det/setFloorShield", this);
  fFloorShieldCmd->SetGuidance("Set the FloorShield thickiness.");
  fFloorShieldCmd->SetParameterName("FloorShield",true);
  fFloorShieldCmd->SetDefaultUnit("m");
  fFloorShieldCmd->SetDefaultValue(fFloorShieldCmd->GetNewDoubleValue("0. m"));
  fFloorShieldCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  fDetector->SetFloorShield(fFloorShieldCmd->GetNewDoubleValue("0. m"));
  std::cout << "DetMess::pre-SetNewVal() FloorShield " << fDetector->GetFloorShield() << std::endl;

  fGDMLCmd = new G4UIcmdWithAString("/sndwrm/det/setGDMLfile", this);
  fGDMLCmd->SetGuidance("Set the GDML file name.");
  fGDMLCmd->SetParameterName("GDMLfile",true);
  fGDMLCmd->SetDefaultValue("");
  fGDMLCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  fDetector->SetGDMLfile("");
  std::cout << "DetMess::pre-SetNewVal() GDML: " << fDetector->GetGDMLfile() << std::endl;

}

void DetectorMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{ 
  //if (command == fTargMatCmd )
  // { fDetector->SetTargetMaterial(newValue);}

  if (command == fFidVolumeCmd )
    {
      fDetector->SetFidVolume(fFidVolumeCmd->GetNew3VectorValue(newValue));
      std::cout << "DetMess::SNV() FidVolume " << fFidVolumeCmd->GetNew3VectorValue(newValue) << std::endl;
    }
  if (command == fFloorShieldCmd )
    {
      fDetector->SetFloorShield(fFloorShieldCmd->GetNewDoubleValue(newValue));
      std::cout << "DetMess::SNV() FloorShield " << fFloorShieldCmd->GetNewDoubleValue(newValue) << std::endl;
    }
  if (command == fGDMLCmd )
    {
      fDetector->SetGDMLfile(newValue);
      std::cout << "DetMess::SNV() GDML file " << newValue << std::endl;
    }
  
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorMessenger::~DetectorMessenger()
{
  delete fDetDirsw;
  delete fFidVolumeCmd;
  delete fDetDir;
  delete fRdecayDir;  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
