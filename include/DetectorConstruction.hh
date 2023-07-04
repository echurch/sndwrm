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
/// \file DetectorConstruction.hh
/// \brief Definition of the DetectorConstruction class
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

#include "MaterialPropertyLoader.hh"

class G4LogicalVolume;
class G4Material;
class DetectorMessenger;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class DetectorConstruction : public G4VUserDetectorConstruction
{
  public:
  
    DetectorConstruction();
   ~DetectorConstruction();

  public:
  
    virtual G4VPhysicalVolume* Construct();
    
    void SetTargetLength (G4double value);
    void SetTargetRadius (G4double value);
    void SetAcrylicLength (G4double value);
    void SetAcrylicRadius (G4double value);
    void SetShieldThickness (G4double value);
    void SetInsetRadius (G4double value);
    void SetTargetMaterial (G4String);
    
    void SetDetectorLength(G4double value);           
    void SetDetectorThickness(G4double value);  
    void SetDetectorRadius(G4double value);  
    void SetDetectorMaterial(G4String);               
         
  void SetSiPMsOnAcrylic(G4bool);
  void SetSiPMsOnCathode(G4bool);
  void SetSiPMSize(G4double);
  void SetSiPMThickness(G4double);
  void SetSiPMPhotoCathodeCoverage(G4double);
  void SetArapucasInCage(G4bool);
  
    void PrintParameters();
    
  public:
      
    G4double GetTargetLength();
    G4double GetTargetRadius();
    G4double GetAcrylicLength();
    G4double GetAcrylicRadius();
    G4double GetShieldThickness();
    G4double SetInsetRadius();
    G4Material* GetTargetMaterial();       
    G4Material* GetShieldMaterial(); 
    G4LogicalVolume* GetLogicTarget();
    G4LogicalVolume* GetLogicSiPM();
    bool GetAPEX() {return fAPEX;};
    G4double GetDetectorLength();
    G4double GetDetectorThickness();
    G4double GetDetectorRadius();
    G4Material* GetDetectorMaterial();                 
    G4LogicalVolume* GetLogicDetector();      
                       
  private:
  
    G4double           fTargetLength; 
    G4double           fTargetRadius;
    G4double           fAcrylicLength; 
    G4double           fAcrylicRadius;
    G4double           fShieldThickness;
    G4double           fG10Thickness;
    G4double           fWoodThickness;
    G4double           fAcrylicThickness;
    G4double           fTPBThickness; 
    G4double           fColdSkinThickness; 
    G4double           fInsetRadius;
    bool               fSiPMsOnAcrylic;
    bool               fSiPMsOnCathode;
    bool               fAPEX;
    G4double           fSiPMSize;
    G4double           fSiPMThickness;
    G4double           fSiPMPhotoCathodeCoverage;
    G4Material*        fTargetMater;
    G4Material*        fShieldMater;
    G4Material*        fG10Mater;
    G4Material*        fWoodMater;
    G4Material*        fColdSkinMater;
    G4Material*        fAcrylicMater;
    G4Material*        fSiPMMater;
    G4Material*        fTPBMater;
    G4Material*        fArapucaMater;
    G4Material*        fAluminumMater;
    G4LogicalVolume*   fLogicTarget;
    G4LogicalVolume*   fLogicShield;
    G4LogicalVolume*   fLogicG10;
    G4LogicalVolume*   fLogicWood;
    G4LogicalVolume*   fLogicColdSkin;
    G4LogicalVolume*   fLogicAcrylic;
    G4LogicalVolume*   fLogicTPB;
    G4LogicalVolume*   fLogicSiPM;
    G4LogicalVolume*   fLogicArapuca;
    G4LogicalVolume*   fLogicArapuca2; // shorter
  
    G4LogicalVolume*   fLogicHBack;
    G4LogicalVolume*   fLogicHIbar;
    G4LogicalVolume*   fLogicHFront;
    G4LogicalVolume*   fLogicHBackEnd;
    G4LogicalVolume*   fLogicHIbarEnd;
    G4LogicalVolume*   fLogicHFrontEnd;
                 
    G4double           fDetectorLength;
    G4double           fDetectorThickness;
    G4double           fDetectorRadius;
    G4Material*        fDetectorMater;
    G4LogicalVolume*   fLogicDetector;
               
    G4double           fWorldLength;
    G4double           fWorldRadius;
    G4Material*        fWorldMater;     
    G4VPhysicalVolume* fPhysiWorld;
                
    DetectorMessenger* fDetectorMessenger;

  private:
    
    void               DefineMaterials();
    G4VPhysicalVolume* ConstructVolumes();     
    void               DetSiPMs(G4String , G4LogicalVolume* );
    void               DetAPEX (G4String , G4LogicalVolume* );
    MaterialPropertyLoader* fMPL;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


#endif

