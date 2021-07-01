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
        
  fTargMatCmd = new G4UIcmdWithAString("/rdecay02/det/setTargetMate",this);
  fTargMatCmd->SetGuidance("Select material of the target");
  fTargMatCmd->SetParameterName("choice",false);
  fTargMatCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  fShieldThicknessCmd =
       new G4UIcmdWithADoubleAndUnit("/rdecay02/det/setShieldThickness", this);
  fShieldThicknessCmd->SetGuidance("Set the Shield Thickness.");
  fShieldThicknessCmd->SetUnitCategory("Length");
  fShieldThicknessCmd->SetParameterName("choice",false);
  fShieldThicknessCmd->AvailableForStates(G4State_PreInit);  

  fInsetRadiusCmd =
       new G4UIcmdWithADoubleAndUnit("/rdecay02/det/setInsetRadius", this);
  fInsetRadiusCmd->SetGuidance("Set the Inset Radius.");
  fInsetRadiusCmd->SetUnitCategory("Length");
  fInsetRadiusCmd->SetParameterName("choice",false);
  fInsetRadiusCmd->AvailableForStates(G4State_PreInit);  

  fTargRadiusCmd =
       new G4UIcmdWithADoubleAndUnit("/rdecay02/det/setTargetRadius", this);
  fTargRadiusCmd->SetGuidance("Set the Target Radius.");
  fTargRadiusCmd->SetUnitCategory("Length");
  fTargRadiusCmd->SetParameterName("choice",false);
  fTargRadiusCmd->AvailableForStates(G4State_PreInit);  
  
  fTargLengthCmd =
       new G4UIcmdWithADoubleAndUnit("/rdecay02/det/setTargetLength", this);
  fTargLengthCmd->SetGuidance("Set the Target Length.");
  fTargLengthCmd->SetUnitCategory("Length");
  fTargLengthCmd->SetParameterName("choice",false);
  fTargLengthCmd->AvailableForStates(G4State_PreInit);

  fAcrylicRadiusCmd =
       new G4UIcmdWithADoubleAndUnit("/rdecay02/det/setAcrylicRadius", this);
  fAcrylicRadiusCmd->SetGuidance("Set the Acrylic Radius.");
  fAcrylicRadiusCmd->SetUnitCategory("Length");
  fAcrylicRadiusCmd->SetParameterName("choice",false);
  fAcrylicRadiusCmd->AvailableForStates(G4State_PreInit);  
  
  fAcrylicLengthCmd =
       new G4UIcmdWithADoubleAndUnit("/rdecay02/det/setAcrylicLength", this);
  fAcrylicLengthCmd->SetGuidance("Set the Acrylic Length.");
  fAcrylicLengthCmd->SetUnitCategory("Length");
  fAcrylicLengthCmd->SetParameterName("choice",false);
  fAcrylicLengthCmd->AvailableForStates(G4State_PreInit);
  

  fDetectMatCmd = new G4UIcmdWithAString("/rdecay02/det/setDetectorMate",this);
  fDetectMatCmd->SetGuidance("Select Material of the Detector.");
  fDetectMatCmd->SetParameterName("choice",false);
  fDetectMatCmd->AvailableForStates(G4State_PreInit,G4State_Idle);  

  fDetectThicknessCmd =
       new G4UIcmdWithADoubleAndUnit("/rdecay02/det/setDetectorThickness",this);
  fDetectThicknessCmd->SetGuidance("Set the Detector Thickness.");
  fDetectThicknessCmd->SetUnitCategory("Length");
  fDetectThicknessCmd->SetParameterName("choice",false);
  fDetectThicknessCmd->AvailableForStates(G4State_PreInit);

  fDetectLengthCmd =
       new G4UIcmdWithADoubleAndUnit("/rdecay02/det/setDetectorLength",this);
  fDetectLengthCmd->SetGuidance("Set the Detector Length.");
  fDetectLengthCmd->SetUnitCategory("Length");
  fDetectLengthCmd->SetParameterName("choice",false);
  fDetectLengthCmd->AvailableForStates(G4State_PreInit);


  fSiPMsOnAcrylicCmd =
       new G4UIcmdWithABool("/rdecay02/det/setDetectorSiPMsOnAcrylic",this);
  fSiPMsOnAcrylicCmd->SetGuidance("Put SiPMs on side Acrylic walls");
  fSiPMsOnAcrylicCmd->SetParameterName("choice",false);
  fSiPMsOnAcrylicCmd->AvailableForStates(G4State_PreInit);
  fSiPMsOnCathodeCmd =
       new G4UIcmdWithABool("/rdecay02/det/setDetectorSiPMsOnCathode",this);
  fSiPMsOnCathodeCmd->SetGuidance("Put SiPMs on side Cathode");
  fSiPMsOnCathodeCmd->SetParameterName("choice",false);
  fSiPMsOnCathodeCmd->AvailableForStates(G4State_PreInit);
  fSiPMSizeCmd =
    new G4UIcmdWithADoubleAndUnit("/rdecay02/det/setDetectorSiPMSize",this);
  fSiPMSizeCmd->SetGuidance("Size of SiPM on side");
  fSiPMSizeCmd->SetUnitCategory("Length");
  fSiPMSizeCmd->SetParameterName("choice",false);
  fSiPMSizeCmd->AvailableForStates(G4State_PreInit);
  fSiPMThicknessCmd =
    new G4UIcmdWithADoubleAndUnit("/rdecay02/det/setDetectorSiPMThickness",this);
  fSiPMThicknessCmd->SetGuidance("Thickness of SiPM");
  fSiPMThicknessCmd->SetUnitCategory("Length");
  fSiPMThicknessCmd->SetParameterName("choice",false);
  fSiPMThicknessCmd->AvailableForStates(G4State_PreInit);
  fSiPMPhotoCathodeCoverageCmd =
    new G4UIcmdWithADouble("/rdecay02/det/setDetectorSiPMPhotoCathodeCoverage",this);
  fSiPMPhotoCathodeCoverageCmd->SetGuidance("Fractional PhotoCathodeCoverage of SiPM");
  fSiPMPhotoCathodeCoverageCmd->SetParameterName("choice",false);
  fSiPMPhotoCathodeCoverageCmd->AvailableForStates(G4State_PreInit);


}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorMessenger::~DetectorMessenger()
{
  delete fTargMatCmd;
  delete fDetectMatCmd;
  delete fAcrylicRadiusCmd;
  delete fAcrylicLengthCmd;
  delete fTargRadiusCmd;
  delete fShieldThicknessCmd;
  delete fInsetRadiusCmd;
  delete fDetectThicknessCmd;
  delete fTargLengthCmd;
  delete fDetectLengthCmd;
  delete fDetDir;
  delete fRdecayDir;  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{ 
  if (command == fTargMatCmd )
   { fDetector->SetTargetMaterial(newValue);}

    // I believe to get desired nesting of volumes, I must go in this order.
  if (command == fShieldThicknessCmd ) 
    {fDetector->SetShieldThickness(fShieldThicknessCmd->GetNewDoubleValue(newValue));}
   
  if (command == fTargLengthCmd ) 
    { fDetector->SetTargetLength(fTargLengthCmd->GetNewDoubleValue(newValue));}
    
  if (command == fTargRadiusCmd ) 
    {fDetector->SetTargetRadius(fTargLengthCmd->GetNewDoubleValue(newValue));}

  if (command == fAcrylicLengthCmd ) 
    { fDetector->SetAcrylicLength(fAcrylicLengthCmd->GetNewDoubleValue(newValue));}
    
  if (command == fAcrylicRadiusCmd ) 
    {fDetector->SetAcrylicRadius(fAcrylicLengthCmd->GetNewDoubleValue(newValue));}


  if (command == fInsetRadiusCmd ) 
    {fDetector->SetInsetRadius(fInsetRadiusCmd->GetNewDoubleValue(newValue));}
    
  if (command == fDetectMatCmd )
    { fDetector->SetDetectorMaterial(newValue);}
    
  if (command == fDetectLengthCmd ) 
    {fDetector->SetDetectorLength(
                     fDetectLengthCmd->GetNewDoubleValue(newValue));}

  if (command == fDetectThicknessCmd ) 
    {fDetector->SetDetectorThickness(
                     fDetectThicknessCmd->GetNewDoubleValue(newValue));}      



  if (command == fSiPMsOnAcrylicCmd ) 
    {fDetector->SetSiPMsOnAcrylic(
                     fSiPMsOnAcrylicCmd->GetNewBoolValue(newValue));}      
  if (command == fSiPMsOnCathodeCmd ) 
    {fDetector->SetSiPMsOnCathode(
                     fSiPMsOnCathodeCmd->GetNewBoolValue(newValue));}      
  if (command == fSiPMSizeCmd ) 
    {fDetector->SetSiPMSize(
                     fSiPMSizeCmd->GetNewDoubleValue(newValue));}      
  if (command == fSiPMThicknessCmd ) 
    {fDetector->SetSiPMThickness(
                     fSiPMThicknessCmd->GetNewDoubleValue(newValue));}      
  if (command == fSiPMPhotoCathodeCoverageCmd ) 
    {fDetector->SetSiPMPhotoCathodeCoverage(
                     fSiPMPhotoCathodeCoverageCmd->GetNewDoubleValue(newValue));}      


}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
