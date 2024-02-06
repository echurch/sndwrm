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
#define DetectorConstruction_h 
#include "G4Material.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "DetectorMessenger.hh"
#include "G4PVPlacement.hh"
#include "G4SubtractionSolid.hh"
#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"
#include "G4OpticalSurface.hh"
#include "G4LogicalSkinSurface.hh"
//class G4LogicalVolume;
//class G4Material;
//class DetectorMessenger;
//class G4Box;

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
    void SetShieldThickness (G4double value);
    void SetInsetRadius (G4double value);
    void SetTargetMaterial (G4String);
    
    void SetDetectorLength(G4double value);           
    void SetDetectorThickness(G4double value);  
    void SetDetectorRadius(G4double value);  
    void SetDetectorMaterial(G4String);               
                   
    void PrintParameters();
    
  public:
    
    G4double GetTargetLength();
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
    G4double           fShieldThickness;
    G4double           fWoodThickness;
    G4double           fInsetRadius;
    G4Material*        fTargetMater;
    G4Material*        fShieldMater;
    G4Material*        fWoodMater;
    G4LogicalVolume*   fLogicTarget;
    G4LogicalVolume*   fLogicShield;
    G4LogicalVolume*   fLogicWood;
                 
    G4double           fDetectorLength;
    G4double           fDetectorThickness;
    G4double           fDetectorRadius;
    G4Material*        fDetectorMater;
    G4LogicalVolume*   fLogicDetector;
               
    G4double           fWorldLength;
    G4double           fWorldRadius;
    G4Material*        fWorldMater;     
    
     G4double fthickness;

  G4double fLatY;
  G4double fLatZ;

  G4double fwindow;
  G4double fwindowY;
  G4double fwindowZ;

  G4double      fWorldSizeX;
  G4double      fWorldSizeY;
  G4double      fWorldSizeZ;
  G4double      fCryostat_x;
  G4double      fCryostat_y;
  G4double      fCryostat_z;
  G4double      fColdSkinThickness;
  G4double      fWarmSkinThickness;
  G4double      fFC_x;
  G4double      fFC_y;
  G4double      fFC_z;
  G4double      fFCOut_x;
  G4double      fFCOut_y;
  G4double      fFCOut_z;
  G4double      fAra_x;
  G4double      fAra_y;
  G4double      fAra_z;
  G4double      fAra_offset;
  G4double      fAras_yspacing;

  G4double      fCathode_x;
  G4double      fCathode_z;
  G4double      fAPA_x;
  G4double      fAPA_y;
  G4double      fAPA_z;

  G4double fptp_width;

  G4double fvert_bar_x;
  G4double fvert_bar_y;
  G4double fvert_bar_z;


// Materials

   G4Material*   fDefaultMaterial;
   G4Material*   fSteel;
   G4Material*   fDUNESteel;
   G4Material*   fAluminium;
   G4Material*   fG10;
   G4Material*   fBase;
   G4Material*   facrylic;
   G4Material*   fPTP;
   G4Material*   fMylar;

   G4VPhysicalVolume* fPhysiWorld;
   G4LogicalVolume*   fLogicWorld;
   G4Box*             fSolidWorld;

   G4VPhysicalVolume* fPhysiVol;
   G4LogicalVolume*   fLogicVol;
   G4Box*             fSolidVol;



    DetectorMessenger* fDetectorMessenger;

  private:
    
    void               DefineMaterials();
    G4VPhysicalVolume* ConstructVolumes();        
    G4VPhysicalVolume* ConstructLine();
 
 
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


#endif

