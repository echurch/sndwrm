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
    
    void SetTubeLength (G4double value);
    void SetTubeThickness (G4double value);
    void SetTubeRadius (G4double value);
    void SetShieldThickness (G4double value);
    void SetShieldRadius (G4double value);
    void SetFlangeThickness (G4double value);
    void SetInsetRadius (G4double value);
    void SetTargetPressure (G4double value);
    void SetTargetMaterial (G4String);
    void SetFlangeMaterial (G4String);
    
    void SetDetectorLength(G4double value);           
    void SetDetectorThickness(G4double value);  
    void SetDetectorRadius(G4double value);  
    void SetDetectorMaterial(G4String);               
                   
    void PrintParameters();
    
  public:
      
    G4double GetTargetLength();
    G4double GetTargetPressure();
    G4double GetTargetRadius();
    G4double GetShieldThickness();
    G4double SetInsetRadius();
    G4Material* GetTargetMaterial();       
    G4Material* GetShieldMaterial(); 
    G4LogicalVolume* GetLogicTarget();
    
    G4double GetDetectorLength();
    G4double GetDetectorThickness();
    G4double GetDetectorRadius();
    G4Material* GetDetectorMaterial();                 
    G4LogicalVolume* GetLogicDetector();      
                       
  private:
  
    G4double           fTargetLength; 
    G4double           fTargetRadius;
    G4double           fTargetPressure;
    G4double           fShieldThickness;
    G4double           fShieldRadius;
    G4double           fTubeThickness;
    G4double           fTubeRadius;
    G4double           fLensRadius;
    G4double           fTubeLength;
    G4double           fFlangeThickness;
    G4double           fBufferThickness;
    G4double           fBufferInnerRadius;
    G4double           fBufferOuterRadius;
    G4double           fBufferLength;
    G4double           fFieldCageThickness;
    G4double           fFieldCageInnerRadius;
    G4double           fFieldCageOuterRadius;
    G4double           fFieldCageLength;
    G4double           fElectrodeThickness;
    G4double           fElectrodeInnerRadius;
    G4double           fElectrodeOuterRadius;
    G4double           fBaffleThickness;
    G4double           fBaffleEndGap;
    G4double           fBaffleInnerRadius;
    G4double           fBaffleOuterRadius;
  

    G4double           fG10Thickness;
    G4double           fInsetRadius;
    G4Material*        fTubeMater;
    G4Material*        fTargetMater;
    G4Material*        fShieldMater;
    G4Material*        fFlangeMater;
    G4Material*        fBufferMater;
    G4Material*        fELPMater;
    G4Material*        fELPPMater;
    G4Material*        fECMater;
    G4Material*        fFieldCageMater;
    G4Material*        fBaffleMater;
    G4Material*        fBaffleFlangeMater;
    G4Material*        fG10Mater;
    G4LogicalVolume*   fLogicTarget;
    G4LogicalVolume*   fLogicDetector;
    G4LogicalVolume*   fLogicTube;
    G4LogicalVolume*   fLogicTubeECN;
    G4LogicalVolume*   fLogicTubeECS;
    G4LogicalVolume*   fLogicShield;
    G4LogicalVolume*   fLogicShieldECN;
    G4LogicalVolume*   fLogicShieldECS;
    G4LogicalVolume*   fLogicFlange;
    G4LogicalVolume*   fLogicFlangeECN;
    G4LogicalVolume*   fLogicFlangeECS;
    G4LogicalVolume*   fLogicBuffer;
    G4LogicalVolume*   fLogicFieldCage;
    G4LogicalVolume*   fLogicBaffle;
    G4LogicalVolume*   fLogicBaffleFlange;
    G4LogicalVolume*   fLogicELP;
    G4LogicalVolume*   fLogicELPP;
    G4LogicalVolume*   fLogicEC;
    G4LogicalVolume*   fLogicG10;
                 
               
    G4double           fWorldLength;
    G4double           fWorldRadius;
    G4Material*        fWorldMater;     
    G4VPhysicalVolume* fPhysiWorld;
                
    DetectorMessenger* fDetectorMessenger;

  private:
    
    void               DefineMaterials();
    G4VPhysicalVolume* ConstructVolumes();     
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


#endif

