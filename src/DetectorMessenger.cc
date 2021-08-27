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

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorMessenger::DetectorMessenger(DetectorConstruction* Det)
:G4UImessenger(), 
 fDetector(Det), fRdecayDir(0), fDetDir(0),
 fTargMatCmd(0), fTargetPressureCmd(0),  fFlangeMatCmd(0), fFlangeThicknessCmd(0), fShieldRadiusCmd(0),
 fShieldThicknessCmd(0), fTubeThicknessCmd(0), fTubeLengthCmd(0), fTubeRadiusCmd(0)
{ 
  fRdecayDir = new G4UIdirectory("/rdecay02/");
  fRdecayDir->SetGuidance("commands specific to this example");
  
  G4bool broadcast = false;
  fDetDir = new G4UIdirectory("/rdecay02/det/",broadcast);
  fDetDir->SetGuidance("detector construction commands");
        
  fTargMatCmd = new G4UIcmdWithAString("/rdecay02/det/setTargetMate",this);
  fTargMatCmd->SetGuidance("Select material of the target");
  fTargMatCmd->SetParameterName("choice",false);
  fTargMatCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  fTargetPressureCmd =
       new G4UIcmdWithADoubleAndUnit("/rdecay02/det/setTargetPressure", this);
  fTargetPressureCmd->SetGuidance("Set the Target Pressure.");
  fTargetPressureCmd->SetUnitCategory("Pressure");
  fTargetPressureCmd->SetParameterName("choice",false);
  fTargetPressureCmd->AvailableForStates(G4State_PreInit);  
  
  fShieldThicknessCmd =
       new G4UIcmdWithADoubleAndUnit("/rdecay02/det/setShieldThickness", this);
  fShieldThicknessCmd->SetGuidance("Set the Shield Thickness.");
  fShieldThicknessCmd->SetUnitCategory("Length");
  fShieldThicknessCmd->SetParameterName("choice",false);
  fShieldThicknessCmd->AvailableForStates(G4State_PreInit);  

  fTubeRadiusCmd =
       new G4UIcmdWithADoubleAndUnit("/rdecay02/det/setTubeRadius", this);
  fTubeRadiusCmd->SetGuidance("Set the Tube Radius.");
  fTubeRadiusCmd->SetUnitCategory("Length");
  fTubeRadiusCmd->SetParameterName("choice",false);
  fTubeRadiusCmd->AvailableForStates(G4State_PreInit);  

  fTubeThicknessCmd =
       new G4UIcmdWithADoubleAndUnit("/rdecay02/det/setTubeThickness", this);
  fTubeThicknessCmd->SetGuidance("Set the Tube Radius.");
  fTubeThicknessCmd->SetUnitCategory("Length");
  fTubeThicknessCmd->SetParameterName("choice",false);
  fTubeThicknessCmd->AvailableForStates(G4State_PreInit);  
  
  fTubeLengthCmd =
       new G4UIcmdWithADoubleAndUnit("/rdecay02/det/setTubeLength", this);
  fTubeLengthCmd->SetGuidance("Set the Tube Length.");
  fTubeLengthCmd->SetUnitCategory("Length");
  fTubeLengthCmd->SetParameterName("choice",false);
  fTubeLengthCmd->AvailableForStates(G4State_PreInit);
  

  fFlangeMatCmd = new G4UIcmdWithAString("/rdecay02/det/setFlangeMate",this);
  fFlangeMatCmd->SetGuidance("Select Material of the Flange.");
  fFlangeMatCmd->SetParameterName("choice",false);
  fFlangeMatCmd->AvailableForStates(G4State_PreInit,G4State_Idle);  

  fFlangeThicknessCmd =
       new G4UIcmdWithADoubleAndUnit("/rdecay02/det/setFlangeThickness",this);
  fFlangeThicknessCmd->SetGuidance("Set the Flange Thickness.");
  fFlangeThicknessCmd->SetUnitCategory("Length");
  fFlangeThicknessCmd->SetParameterName("choice",false);
  fFlangeThicknessCmd->AvailableForStates(G4State_PreInit);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorMessenger::~DetectorMessenger()
{
  delete fTargMatCmd;
  delete fTargetPressureCmd;
  delete fTubeRadiusCmd;
  delete fFlangeMatCmd;
  delete fFlangeThicknessCmd;
  delete fShieldThicknessCmd;
  delete fShieldRadiusCmd;
  delete fTubeLengthCmd;

  delete fDetDir;
  delete fRdecayDir;  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{ 
  if (command == fTargMatCmd )
   { fDetector->SetTargetMaterial(newValue);}
  if (command == fTargetPressureCmd ) 
    {fDetector->SetTargetPressure(fTargetPressureCmd->GetNewDoubleValue(newValue));}

  if (command == fTubeThicknessCmd ) 
    {fDetector->SetTubeThickness(fTubeThicknessCmd->GetNewDoubleValue(newValue));}
  if (command == fTubeRadiusCmd ) 
    {fDetector->SetTubeRadius(fTubeRadiusCmd->GetNewDoubleValue(newValue));}
  if (command == fTubeLengthCmd ) 
    {fDetector->SetTubeLength(fTubeLengthCmd->GetNewDoubleValue(newValue));}

  if (command == fShieldThicknessCmd ) 
    {fDetector->SetShieldThickness(fShieldThicknessCmd->GetNewDoubleValue(newValue));}
   

  if (command == fShieldRadiusCmd ) 
    {fDetector->SetShieldRadius(fShieldRadiusCmd->GetNewDoubleValue(newValue));}
    
  if (command == fFlangeMatCmd )
    { fDetector->SetFlangeMaterial(newValue);}
  if (command == fFlangeThicknessCmd )
    { fDetector->SetFlangeThickness(fFlangeThicknessCmd->GetNewDoubleValue(newValue));}

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
