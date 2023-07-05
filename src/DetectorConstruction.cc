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
/// \file DetectorConstruction.cc
/// \brief Implementation of the DetectorConstruction class
//
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "DetectorConstruction.hh"
#include "DetectorMessenger.hh"
#include "MaterialPropertyLoader.hh"

#include "G4Material.hh"
#include "G4NistManager.hh"

#include "G4Tubs.hh"
#include "G4Box.hh"
#include "G4SubtractionSolid.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"

#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"
#include "G4RunManager.hh"

#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "G4UnitsTable.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction()
:G4VUserDetectorConstruction(),
 fTargetMater(0), fLogicTarget(0), 
 fDetectorMater(0), fLogicDetector(0), 
 fWorldMater(0), fPhysiWorld(0),
 fDetectorMessenger(0)
{
  fTargetLength      = 1*cm; 
  fTargetRadius      = 0.5*cm;
  fDetectorLength    = 5*cm; 
  fDetectorThickness = 2*cm;
  fInsetRadius = 0.0*cm;
  fDetectorRadius = 0.0*cm;

  // 13-Feb-2020 Make the acrylic walls go out 1m beyond fiducial.
  fAcrylicRadius    = 3.0*m; fAcrylicRadius += 0.1*m;
  fAcrylicLength    = 40.0*m;  fAcrylicLength += 2.0*m;
  fAcrylicThickness = 5.0 *cm;
  /*
  fAcrylicThickness = 15.0 *cm;
  fAcrylicThickness = 20.0 *cm;
  */

  fTPBThickness = 3.0 *um; // 13-Jan-2021, ala DEAP
  fTPBThickness = 0.0 ; //  14-Jan-2021. Instead of simulating WLS-TPB, turn on reflective Acrylic. 

  fWorldLength = std::max(fTargetLength,fDetectorLength);
  fWorldRadius = fTargetRadius + 1.0*m;
      
  // These 5 read in from g4mac file.
  fSiPMsOnAcrylic = true;
  fSiPMsOnCathode = false;
  fSiPMSize = 24.* cm;
  fSiPMThickness = 1. * cm;
  fSiPMPhotoCathodeCoverage = 0.8;

  fAPEX = true;
  DefineMaterials();
    
  fDetectorMessenger = new DetectorMessenger(this);
  fMPL = new MaterialPropertyLoader();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction()
{ delete fDetectorMessenger;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct()
{

  return ConstructVolumes();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::DefineMaterials()
{
  // build materials
  //
  G4Element* H  = new G4Element("Hydrogen", "H", 1, 1.*g/mole);
  G4Element* N  = new G4Element("Nitrogen", "N", 7, 14.01*g/mole);
  G4Element* O  = new G4Element("Oxygen",   "O", 8, 16.00*g/mole);
  G4Element* C  = new G4Element("Carbon",   "C", 6, 12.00*g/mole);
  G4Element* Si  = new G4Element("Silicon", "Si", 14, 28.0855*g/mole);
  G4Element* Ni  = new G4Element("Nickel",  "Ni", 28, 58.69*g/mole);
  G4Element* Cr  = new G4Element("Chromium", "Cr", 24, 51.996*g/mole);
  G4Element* Fe  = new G4Element("Iron",     "Fe", 26, 55.84*g/mole);
  G4Element* Mn  = new G4Element("Manganese","Mn", 25, 54.94*g/mole);
  G4Element* Cu  = new G4Element("Copper","Cu", 29, 63.55*g/mole);

  G4int ncomponents; G4double fractionmass;      
  G4Material* Air20 = new G4Material("Air", 1.205*mg/cm3, ncomponents=2,
				     kStateGas, 293.*kelvin, 1.*atmosphere);
  Air20->AddElement(N, fractionmass=0.7);
  Air20->AddElement(O, fractionmass=0.3);

  G4int number_of_atoms;
  G4Material* H2O = new G4Material("Water", 1.0*g/cm3, ncomponents=2,
				     kStateLiquid, 293.*kelvin, 1.*atmosphere);
  H2O->AddElement(H, number_of_atoms=2);
  H2O->AddElement(O, number_of_atoms=1);


  //Plastic. From cbsim 
  G4Material *plastic=new G4Material("Poly",0.96*g/cm3, ncomponents=2);
  plastic->AddElement(C,number_of_atoms=2);
  plastic->AddElement(H,number_of_atoms=4);

  //Foam. From https://indico.fnal.gov/event/20144/session/19/contribution/267/material/slides/1.pdf 
  G4Material *foam=new G4Material("Foam",0.09*g/cm3, ncomponents=4);
  foam->AddElement(C,number_of_atoms=54);
  foam->AddElement(O,number_of_atoms=15);
  foam->AddElement(N,number_of_atoms=4);
  foam->AddElement(H,number_of_atoms=60);

  G4Material *wood=new G4Material("Wood",0.5*g/cm3, ncomponents=3);
  wood->AddElement(C,number_of_atoms=50);
  wood->AddElement(O,number_of_atoms=44);
  wood->AddElement(H,number_of_atoms=6);

  G4Material *ss=new G4Material("DUNESteel"/*SS407L*/,7.93*g/cm3, ncomponents=7);
  ss->AddElement(Fe,fractionmass=95.8/100.);
  ss->AddElement(Mn,fractionmass=1.8/100.);
  ss->AddElement(Ni,fractionmass=0.8/100.);
  ss->AddElement(Si,fractionmass=0.6/100.);
  ss->AddElement(Cu,fractionmass=0.5/100.);
  ss->AddElement(Cr,fractionmass=0.3/100.);
  ss->AddElement(C, fractionmass=0.2/100.);
  

  G4Material* g10 =  new G4Material("G10", 1.700*g/cm3, ncomponents=4); // NemaG10
  g10->AddElement(Si, number_of_atoms=1);
  g10->AddElement(O , number_of_atoms=2);
  g10->AddElement(C , number_of_atoms=3);
  g10->AddElement(H , number_of_atoms=3);

  // PMMA C5H8O2 ( Acrylic )
  G4Material* Acrylic = new G4Material("Acrylic", 1.19*g/cm3, ncomponents=3);
  Acrylic->AddElement(C, number_of_atoms=5);
  Acrylic->AddElement(H, number_of_atoms=8);
  Acrylic->AddElement(O, number_of_atoms=2);

  G4Material* TPB = new G4Material("TPB",1.079*g/cm3,2);
  TPB->AddElement (C, 28);
  TPB->AddElement (H, 22);


  //
  fWorldMater = Air20;
  fTargetMater   = new G4Material("Argon", 18, 40.00*g/mole, 1.4 *g/cm3, kStateLiquid,  77.*kelvin, 1.*atmosphere);
  fDetectorMater = new G4Material("Argon", 18, 40.00*g/mole, 1.4 *g/cm3, kStateLiquid,  77.*kelvin, 1.*atmosphere);


  G4Material* SiPM = new G4Material("SiPM", 14, 28.086*g/mole, 2.331 *g/cm3, kStateSolid,  300.*kelvin, 1.*atmosphere);

  // quartz (SiO2, crystalline)
  G4Material* Arapuca = new G4Material("Arapuca", 2.64 *g/cm3, ncomponents=2);
  Arapuca-> AddElement(Si, number_of_atoms=1);
  Arapuca-> AddElement(O,  number_of_atoms=2);

  // Aluminium
  G4Material* Aluminium = new G4Material("Aluminium", 13.,26.98*g/mole, 2.7*g/cm3);

  
  //  fShieldMater = H2O;
  fShieldMater =  foam /*plastic*/ ;
  fWoodMater =  wood ;
  fG10Mater = g10;
  fColdSkinMater = ss;
  fAcrylicMater = Acrylic;
  fTPBMater = TPB;
  fSiPMMater = SiPM;
  fArapucaMater = Arapuca;
  fAluminumMater = Aluminium;
  /*

    Target is a Box nested inside World box. Detector is a Box inset from Target walls, nested Target, as tall (Length) as Target.
    EC, 20-Oct-2019.

   */

  


}


G4VPhysicalVolume* DetectorConstruction::ConstructVolumes()
{



  /*
  G4Tubs*
  sWorld = new G4Tubs("World",                                 //name
                 0.,fWorldRadius, 0.5*fWorldLength, 0.,twopi); //dimensions  
  */
  fWorldRadius = 2*fTargetRadius;
  fWorldLength = 2*fTargetLength;

  G4Box* sWorld = new G4Box("World", fWorldRadius,fWorldRadius,fWorldLength/2.);
               
  G4LogicalVolume*
  lWorld = new G4LogicalVolume(sWorld,                  //shape
                             fWorldMater,               //material
                             "World");                  //name

  fPhysiWorld = new G4PVPlacement(0,                    //no rotation
                            G4ThreeVector(),            //at (0,0,0)
                            lWorld,                     //logical volume
                            "World",                    //name
                            0,                          //mother volume
                            false,                      //no boolean operation
                            0);                         //copy number



  // Shield. Now foam+wood outside Target,if I use a positive fShieldThickness.
  // And launch neutrons from TargetRadius+Foam+WoodRadius
  G4Box* sInShield = new G4Box("OutShield",fTargetRadius, fTargetRadius, fTargetLength/2);
  G4Box* sOutShield = new G4Box("InShield", fTargetRadius+fShieldThickness, fTargetRadius+fShieldThickness, fTargetLength/2.+fShieldThickness);
  G4SubtractionSolid *sShield = new G4SubtractionSolid("Shield",sOutShield, sInShield);  

  fLogicShield = new G4LogicalVolume(sShield,       //shape
                             fShieldMater,            //material
                             "Shield");               //name
                               
           new G4PVPlacement(0,                         //no rotation
			     G4ThreeVector(0.,0.,0.),  // fWorldLength/2.-1*fDetectorLength/2.),             //at (0,0,0)
                           fLogicShield,              //logical volume
                           "Shield",                  //name
			   lWorld,                      //mother  volume
                           false,                       //no boolean operation
                           0);                          //copy number

  fWoodThickness = 2.4*cm;
  G4Box* sInWood = new G4Box("OutWood",fTargetRadius+fShieldThickness, fTargetRadius+fShieldThickness, fTargetLength/2+fShieldThickness);
  G4Box* sOutWood = new G4Box("InWood", fTargetRadius+fWoodThickness+fShieldThickness, fTargetRadius+fWoodThickness+fShieldThickness, fTargetLength/2.+fWoodThickness+fShieldThickness);
  G4SubtractionSolid *sWood = new G4SubtractionSolid("Wood",sOutWood, sInWood);  

  fLogicWood = new G4LogicalVolume(sWood,       //shape
                             fWoodMater,            //material
                             "Wood");               //name
                               
           new G4PVPlacement(0,                         //no rotation
			     G4ThreeVector(0.,0.,0.),  // fWorldLength/2.-1*fDetectorLength/2.),             //at (0,0,0)
                           fLogicWood,              //logical volume
                           "Wood",                  //name
			   lWorld,                      //mother  volume
                           false,                       //no boolean operation
                           0);                          //copy number





  // And decay K40s from box this thick
  fG10Thickness = 3.0*mm;
  /*  G10 outside Shield
  G4Box* sInG10 = new G4Box("InG10",fTargetRadius+fWoodThickness+fShieldThickness, fTargetRadius+fWoodThickness+fShieldThickness, fTargetLength/2+fWoodThickness+fShieldThickness);
  G4Box* sOutG10 = new G4Box("OutG10", fTargetRadius+fG10Thickness+fWoodThickness+fShieldThickness, fTargetRadius+fG10Thickness+fWoodThickness+fShieldThickness, fTargetLength/2.+fG10Thickness+fWoodThickness+fShieldThickness);
  G4SubtractionSolid *sG10 = new G4SubtractionSolid("G10",sOutG10, sInG10);  
  */

  //ColdSkin
  fColdSkinThickness = - 1.2 * cm;
  G4Box* sOutColdSkin = new G4Box("OutColdSkin",fTargetRadius, fTargetRadius, fTargetLength/2);
  G4Box* sInColdSkin = new G4Box("InColdSkin", fTargetRadius+fColdSkinThickness, fTargetRadius+fColdSkinThickness, fTargetLength/2.+fColdSkinThickness);
  G4SubtractionSolid *sColdSkin = new G4SubtractionSolid("ColdSkin",sOutColdSkin, sInColdSkin);  

  fLogicColdSkin = new G4LogicalVolume(sColdSkin,       //shape
                             fColdSkinMater,            //material
                             "ColdSkin");               //name
                               
           new G4PVPlacement(0,                         //no rotation
			     G4ThreeVector(0.,0.,0.),  // fWorldLength/2.-1*fDetectorLength/2.),             //at (0,0,0)
                           fLogicColdSkin,              //logical volume
                           "ColdSkin",                  //name
			   lWorld,                      //mother  volume
                           false,                       //no boolean operation
                           0);                          //copy number

                            
  // Target
  //
  /*
  G4Tubs* 
  sTarget = new G4Tubs("Target",                                   //name
                  0., fTargetRadius, 0.5*fTargetLength, 0.,twopi); //dimensions
  */
	   // With a negative fColdSkinThickness this will put the target inside the skin. EC, 1-Nov-2019.

  G4Box* sTarget = new G4Box("Target", fTargetRadius+fColdSkinThickness, fTargetRadius+fColdSkinThickness, fTargetLength/2.+fColdSkinThickness);



  fLogicTarget = new G4LogicalVolume(sTarget,           //shape
                             fTargetMater,              //material
                             "Target");                 //name
                               
           new G4PVPlacement(0,                         //no rotation
			     G4ThreeVector(0.,0.,0.),             //at (0,0,0)
                           fLogicTarget,                //logical volume
                           "Target",                    //name
                           lWorld,                      //mother  volume
                           false,                       //no boolean operation
                           0);                          //copy number



  // G10 inside Shield, if I use a positive fShieldThickness. Remember fColdSkinThickness is negative.
	   /*
  G4Box* sInG10 = new G4Box("InG10",fTargetRadius-fG10Thickness+fColdSkinThickness, fTargetRadius-fG10Thickness+fColdSkinThickness, fTargetLength/2-fG10Thickness+fColdSkinThickness);
  G4Box* sOutG10 = new G4Box("OutG10", fTargetRadius+fColdSkinThickness, fTargetRadius+fColdSkinThickness, fTargetLength/2.+fColdSkinThickness);
  G4SubtractionSolid *sG10 = new G4SubtractionSolid("G10",sOutG10, sInG10);
	   */
	   // simple G10 rectangle hung from top, EC, 30-Jun-2023,
  G4Box* sG10 = new G4Box("G10",(fTargetRadius-fG10Thickness+fColdSkinThickness)/2.,fG10Thickness/2.,fTargetLength/2-fG10Thickness+fColdSkinThickne\
ss);

  fLogicG10 = new G4LogicalVolume(sG10,       //shape
                             fG10Mater,            //material
                             "G10");               //name
                               
         new G4PVPlacement(0,                         //no rotation
			   //			     G4ThreeVector(0.,0.,0.),  // fWorldLength/2.-1*fDetectorLength/2.),             //at (0,0,0)
			   G4ThreeVector(0.,fTargetRadius-fG10Thickness+fColdSkinThickness,0.),  // fWorldLength/2.-1*fDetectorLength/2.),             //at (0,0,0)
                           fLogicG10,              //logical volume
                           "G10",                  //name
			   fLogicTarget,                      //mother  volume
                           false,                       //no boolean operation
                           0,0);                          //copy number

  
  std::cout << "fInsetRadius is " << fInsetRadius << std::endl;
  std::cout << "fShieldThickness is " << fShieldThickness << std::endl;
  std::cout << "fWoodThickness is " << fWoodThickness << std::endl;
  std::cout << "fColdSkinThickness is " << fColdSkinThickness << std::endl;
  std::cout << "fG10Thickness is " << fG10Thickness << std::endl;

  
  fDetectorRadius = fTargetRadius-fInsetRadius;
  /*
  G4Box* sDetector = new G4Box("Detector",fDetectorRadius, fDetectorRadius, fDetectorLength/2);

  
  fLogicDetector = new G4LogicalVolume(sDetector,       //shape
                             fDetectorMater,            //material
                             "Detector");               //name
                               
           new G4PVPlacement(0,                         //no rotation
			     G4ThreeVector(0.,0.,0.),  // fWorldLength/2.-1*fDetectorLength/2.),             //at (0,0,0)
                           fLogicDetector,              //logical volume
                           "Detector",                  //name
			   lWorld,                      //mother  volume
                           false,                       //no boolean operation
                           0);                          //copy number
  */
	   
  // in fact, we want y to go top to bottom to allow charge drift and E-field to be unperturbed. Indeed, see the 13-Feb changes.
  // that might mean even taller than fTargRad to come in under cryoskin on top and bottom, but this should be good approx insofar as n's go. eps to prevent overlap.

  G4Box* sOutAcrylic = new G4Box("InAcrylic",fAcrylicRadius, fTargetRadius , fAcrylicLength/2);
  // 13-Feb-2020 make the top/bottom 0 thickness. 
  // Note we shave fTPBThickness here off the Acrylic thickness. Will use it for TPB below. EC, 13-Jan-2020
  G4Box* sInAcrylic = new G4Box("OutAcrylic", fAcrylicRadius-fAcrylicThickness+fTPBThickness, fTargetRadius, fAcrylicLength/2.-fAcrylicThickness+fTPBThickness);
  G4SubtractionSolid *sAcrylic = new G4SubtractionSolid("Acrylic",sOutAcrylic, sInAcrylic);  

  fLogicAcrylic = new G4LogicalVolume(sAcrylic,       //shape
                             fAcrylicMater,            //material
                             "Acrylic");               //name
  if (fSiPMsOnAcrylic)
         new G4PVPlacement(0,                         //no rotation
			     G4ThreeVector(0.,0.,0.),  // fWorldLength/2.-1*fDetectorLength/2.),             //at (0,0,0)
                           fLogicAcrylic,              //logical volume
                           "Acrylic",                  //name
			   fLogicTarget, //lWorld,                      //mother  volume
                           false,                       //no boolean operation
                           0,0);                          //copy number

  G4Box* sInTPB = new G4Box("InTPB",fAcrylicRadius-fAcrylicThickness, fTargetRadius , fAcrylicLength/2-fAcrylicThickness);
  // 13-Feb-2020 make the top/bottom 0 thickness.
  G4Box* sOutTPB = new G4Box("OutTPB", fAcrylicRadius-fAcrylicThickness+fTPBThickness, fTargetRadius, fAcrylicLength/2.-fAcrylicThickness+fTPBThickness);
  G4SubtractionSolid *sTPB = new G4SubtractionSolid("TPB",sOutTPB, sInTPB);  

  fLogicTPB = new G4LogicalVolume(sTPB,       //shape
                             fTPBMater,            //material
                             "TPB");               //name
                               
         new G4PVPlacement(0,                         //no rotation
			     G4ThreeVector(0.,0.,0.),  // fWorldLength/2.-1*fDetectorLength/2.),             //at (0,0,0)
                           fLogicTPB,              //logical volume
                           "TPB",                  //name
			   lWorld,                      //mother  volume
                           false,                       //no boolean operation
                           0);                          //copy number


  if (fSiPMsOnAcrylic)
    DetSiPMs("Acrylic",fLogicTarget);
  if (fSiPMsOnCathode)
    DetSiPMs("Cathode",fLogicTarget);

  if (fAPEX)
    DetAPEX("Perimeter",fLogicTarget);

  
  // Now for the purpose of tracking optical photons we do the following to the Argon and TPB to endow them w optical physics properties.
  // Must wait till this late, cuz MLP works by looping over all Logical Volumes which are only just now established.
  // Get the logical volume store and assign material properties. MaterialPropLoader() is borrowed, heavily-edited from LArSoft.   

  fMPL->SetPropertiesFromServices();  // fills local LArprop class with hard-coded data cutnpasted from fcl file.
  fMPL->GetPropertiesFromServices();  // Shoves these into local MaterialTables
  fMPL->UpdateGeometry(G4LogicalVolumeStore::GetInstance()); // Finally, loads properties into G4MaterialProperties



  std::cout <<  "DetConst::ConstructVolumes(): Argon " << " Fast Intensity" << std::endl;
  fTargetMater->GetMaterialPropertiesTable()->GetProperty("FASTCOMPONENT") ->DumpValues();
  std::cout << "DetConst::ConstructVolumes(): Argon " << " Slow Intensity" << std::endl;
  fTargetMater->GetMaterialPropertiesTable()->GetProperty("SLOWCOMPONENT") ->DumpValues();
  std::cout << "DetConst::ConstructVolumes(): Argon " << " All Properties" << std::endl;
  fTargetMater->GetMaterialPropertiesTable()->DumpTable();
  std::cout << "Dumping TPB MaterialPropertiesTable:" << std::endl;  
  fTPBMater->GetMaterialPropertiesTable()->DumpTable();
  /*
  if (fTPBMater->GetMaterialPropertiesTable()->ConstPropertyExists("WLSTIMECONSTANT") )
    {
      std::cout <<  "DetConst::ConstructVolumes(): TPB " << "Re-emission time const [ns]" << 
	fTPBMater->GetMaterialPropertiesTable()->GetConstProperty("WLSTIMECONSTANT") << std::endl;
    }
    else
      std::cout <<  "DetConst::ConstructVolumes(): TPB WLSTIMECONST weirdly doesn't exist" << std::endl;
      
  std::cout <<  "DetConst::ConstructVolumes(): TPB " << " ABS length [m]" << std::endl;
  fTPBMater->GetMaterialPropertiesTable()->GetProperty("WLSABSLENGTH") ->DumpValues();
  std::cout << "DetConst::ConstructVolumes(): TPB " << " EMIS Spect" << std::endl;
  fTPBMater->GetMaterialPropertiesTable()->GetProperty("WLSCOMPONENT") ->DumpValues();
  std::cout <<  "DetConst::ConstructVolumes(): TPB " << "Refr Index [eV]" << std::endl;
  fTPBMater->GetMaterialPropertiesTable()->GetProperty("RINDEX") ->DumpValues();
  */

  std::cout << "Dumping Acrylic MaterialPropertiesTable:" << std::endl;
  fAcrylicMater->GetMaterialPropertiesTable()->DumpTable();
  std::cout << "Dumping SiPM MaterialPropertiesTable:" << std::endl;
  fSiPMMater->GetMaterialPropertiesTable()->DumpTable();
  std::cout << "Dumping Arapuca MaterialPropertiesTable:" << std::endl;
  fArapucaMater->GetMaterialPropertiesTable()->DumpTable();
  


  PrintParameters();
  
  //always return the root volume
  //
  return fPhysiWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::PrintParameters()
{
  G4cout << "\n World : Length = " << G4BestUnit(fWorldLength,"Length")
         << " Radius = " << G4BestUnit(fWorldRadius,"Length")  
         << " Material = " << fWorldMater->GetName();
  G4cout << "\n Target : Length = " << G4BestUnit(fTargetLength,"Length")
         << " Radius = " << G4BestUnit(fTargetRadius,"Length")  
         << " Material = " << fTargetMater->GetName();
  G4cout << "\n Detector : Length = " << G4BestUnit(fDetectorLength,"Length")
         << " Radius = " << G4BestUnit(fDetectorRadius,"Length")  
    //         << " Thickness = " << G4BestUnit(fDetectorThickness,"Length")  
         << " Material = " << fDetectorMater->GetName() << G4endl;          


  G4cout << "\n Shield : Thickness = " << G4BestUnit(fShieldThickness,"Length")
         << " Material = " << fShieldMater->GetName() << G4endl;          
  G4cout << "\n Wood : Thickness = " << G4BestUnit(fWoodThickness,"Length")
         << " Material = " << fWoodMater->GetName() << G4endl;          
  G4cout << "\n G10 : Thickness = " << G4BestUnit(fG10Thickness,"Length")
         << " Material = " << fG10Mater->GetName() << G4endl;          
  G4cout << "\n ColdSkin : Thickness = " << G4BestUnit(fColdSkinThickness,"Length")
         << " Material = " << fColdSkinMater->GetName() << G4endl;          
  G4cout << "\n Acrylic : Thickness = " << G4BestUnit(fAcrylicThickness,"Length")
         << " x -size  = " << G4BestUnit(fAcrylicRadius,"Length")  
         << " z Length = " << G4BestUnit(fAcrylicLength,"Length")  
         << " Material = " << fAcrylicMater->GetName() << G4endl;          
  G4cout << "\n SiPM : Thickness = " << G4BestUnit(fSiPMThickness,"Length")
         << " Size = " << G4BestUnit(fSiPMSize,"Length")  
         << " Material = " << fSiPMMater->GetName() << G4endl;          
  G4cout << "\n TPB : Thickness = " << G4BestUnit(fTPBThickness,"Length")
         << " Material = " << fTPBMater->GetName() << G4endl;          

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetTargetMaterial(G4String materialChoice)
{
  // search the material by its name
  G4Material* pttoMaterial =
     G4NistManager::Instance()->FindOrBuildMaterial(materialChoice);   
  
  if (pttoMaterial) { 
    fTargetMater = pttoMaterial;
    if(fLogicTarget) { fLogicTarget->SetMaterial(fTargetMater); }
    G4RunManager::GetRunManager()->PhysicsHasBeenModified();
  } else {
    G4cout << "\n--> warning from DetectorConstruction::SetTargetMaterial : "
           << materialChoice << " not found" << G4endl;
  }              
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetDetectorMaterial(G4String materialChoice)
{
  // search the material by its name
  G4Material* pttoMaterial =
     G4NistManager::Instance()->FindOrBuildMaterial(materialChoice);   
  
  if (pttoMaterial) { 
    fDetectorMater = pttoMaterial;
    if(fLogicDetector) { fLogicDetector->SetMaterial(fDetectorMater); }
    G4RunManager::GetRunManager()->PhysicsHasBeenModified();
  } else {
    G4cout << "\n--> warning from DetectorConstruction::SetDetectorMaterial : "
           << materialChoice << " not found" << G4endl;
  }              
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetTargetRadius(G4double value)
{
  fTargetRadius = value;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetAcrylicRadius(G4double value)
{
  fAcrylicRadius = value;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetAcrylicLength(G4double value)
{
  fAcrylicLength = value;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetInsetRadius(G4double value)
{
  fInsetRadius = value;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetTargetLength(G4double value)
{
  fTargetLength = value;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}

void DetectorConstruction::SetShieldThickness(G4double value)
{
  fShieldThickness = value;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetDetectorThickness(G4double value)
{
  fDetectorThickness = value;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetDetectorLength(G4double value)
{
  fDetectorLength = value;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}

void DetectorConstruction::SetDetectorRadius(G4double value)
{
  fDetectorRadius = value;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double DetectorConstruction::GetTargetRadius()
{

  return fTargetRadius;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double DetectorConstruction::GetAcrylicRadius()
{
  // 5cm is for the G10 and/or coldcryoskin
  assert((fTargetRadius-5*cm)>fAcrylicRadius && "fAcrylicLength shorter than fTargetRadius");
  return fAcrylicRadius;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double DetectorConstruction::GetTargetLength()
{
  return fTargetLength;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double DetectorConstruction::GetAcrylicLength()
{
  // 5cm is for the G10 and/or coldcryoskin
  assert((fTargetLength-5*cm)>fAcrylicLength && "fAcrylicLength shorter than fTargetLength");
  return fAcrylicLength;
}

G4double DetectorConstruction::GetShieldThickness()
{
  return fShieldThickness;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4Material* DetectorConstruction::GetTargetMaterial()
{
  return fTargetMater;
}

G4Material* DetectorConstruction::GetShieldMaterial()
{
  return fShieldMater;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4LogicalVolume* DetectorConstruction::GetLogicTarget()
{
  return fLogicTarget;
}

G4LogicalVolume* DetectorConstruction::GetLogicSiPM()
{
  return fLogicSiPM;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double DetectorConstruction::GetDetectorLength()
{
  return fDetectorLength;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double DetectorConstruction::GetDetectorThickness()
{
  return fDetectorThickness;
}

G4double DetectorConstruction::GetDetectorRadius()
{
  return fDetectorRadius;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4Material* DetectorConstruction::GetDetectorMaterial()
{
  return fDetectorMater;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4LogicalVolume* DetectorConstruction::GetLogicDetector()
{
  return fLogicDetector;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstruction::SetSiPMsOnAcrylic(G4bool value)
{
  fSiPMsOnAcrylic = value;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstruction::SetSiPMsOnCathode(G4bool value)
{
  fSiPMsOnCathode = value;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstruction::SetSiPMThickness(G4double value)
{
  fSiPMThickness = value;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstruction::SetSiPMSize(G4double value)
{
  fSiPMSize = value;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}
void DetectorConstruction::SetSiPMPhotoCathodeCoverage(G4double value)
{
  fSiPMPhotoCathodeCoverage = value;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstruction::SetArapucasInCage(G4bool value)
{
  fAPEX = value;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}



void DetectorConstruction::DetSiPMs(G4String component, G4LogicalVolume* logiMother)
{
  

  const double halfSpacing = (sqrt(fSiPMSize*fSiPMSize/fSiPMPhotoCathodeCoverage)- fSiPMSize)/2. ;

  std::cout << "fTargetRadius, fAcrylicLength, fAcrylicRadius, halfSpacing, SiPMSize+2*halfSpacing" << fTargetRadius << ", " << fAcrylicLength << ", " << fAcrylicRadius << ", "<< halfSpacing <<  ", " << fSiPMSize+2*halfSpacing<< std::endl;

  double area = 2*fTargetRadius*fAcrylicLength;
  double areaL = 4*fTargetRadius*fAcrylicRadius;
  int numLateral = int((2*fAcrylicRadius/(fSiPMSize+2*halfSpacing) ));
  int numLength = int((fAcrylicLength/(fSiPMSize+2*halfSpacing)) );
  int numHeight = int((2*fTargetRadius/(fSiPMSize+2*halfSpacing) ));

  G4Box* sSiPM = new G4Box("SiPM", fSiPMThickness/2, fSiPMSize/2, fSiPMSize/2);
  fLogicSiPM = new G4LogicalVolume(sSiPM,       //shape  
				   fSiPMMater,    //material  
				   "SiPM");      //name                                                                                                                              
  double offsetLat = 0;
  if (!(numLateral%2))
    offsetLat = 0.5;
  double offsetHeight = 0;
  if (!(numHeight%2))
    offsetHeight = 0.5;
  double offsetLength = 0;
  if (!(numLength%2)) 
    offsetLength = 0.5;
  int nSiPM(0);
  int nSiPM2(numHeight*numLength);
  int nSiPM3(2*numHeight*numLength);
  int nSiPM4(2*numHeight*numLength+numHeight*numLength);
  int nSiPMC(0);

  G4RotationMatrix* rotationMatrix = new G4RotationMatrix();

  if (component=="Acrylic") 
    {
      std::cout << "DetectorConstruction::DetSiPMs: Placing " << numHeight*numLength << " SiPMs just inside x Acrylic each wall, each of area " << area/1E6 << " [m^2] of VD detector." << std::endl;


      for (int ii=-numHeight/2; ii<=numHeight/2; ii++)
      {
	for (int jj=-numLength/2; jj<=numLength/2; jj++)
	  {

	    if (!(numHeight%2) && ii==int(numHeight/2)) // if numLateral is even we need to stagger placement by 0.5 integer and skip the last one.
	      continue;
	    if (!(numLength%2) && jj==int(numLength/2)) // if numLength is even we need to stagger placement by 0.5 integer and skip the last one.
	      continue;

 	    // neg x Acrylic wall
	    new G4PVPlacement(0, 
			      G4ThreeVector(-fAcrylicRadius+fAcrylicThickness+fTPBThickness+fSiPMThickness/2,(ii+offsetHeight)*(fSiPMSize+2*halfSpacing),(jj+offsetLength)*(fSiPMSize+2*halfSpacing)),   //cm
			      fLogicSiPM,              //logical volume                         
			      "SiPM",                  //name                         
			      logiMother,                      //mother  volume
			      false,                       //no boolean operation
			      nSiPM++,0);                          //copy number   
	    // pos x Acrylic wall
	    new G4PVPlacement(0, 
			      G4ThreeVector(+fAcrylicRadius-fAcrylicThickness-fTPBThickness-fSiPMThickness/2,(ii+offsetHeight)*(fSiPMSize+2*halfSpacing),(jj+offsetLength)*(fSiPMSize+2*halfSpacing)),  
			      fLogicSiPM,              //logical volume                         
			      "SiPM",                  //name                             
			      logiMother,                      //mother  volume
			      false,                       //no boolean operation
			      nSiPM2++,0);                          //copy number

	  }
      }

      std::cout << "DetectorConstruction::DetSiPMs: Placing " << numLateral*numHeight << " SiPMs just inside each z Acrylic wall, each of area " << areaL/1E6 << " [m^2] of VD detector." << std::endl;
      *rotationMatrix = G4RotationMatrix();
      rotationMatrix->rotateY(90.*deg);
      for (int ii=-numHeight/2; ii<=numHeight/2; ii++)
      {
	for (int jj=-numLateral/2; jj<=numLateral/2; jj++)
	  {

	    if (!(numHeight%2) && ii==int(numHeight/2)) // if numLateral is even we need to stagger placement by 0.5 integer and skip the last one.
	      continue;
	    if (!(numLateral%2) && jj==int(numLateral/2)) // if numLength is even we need to stagger placement by 0.5 integer and skip the last one.
	      continue;

 	    // neg z Acrylic wall
	    new G4PVPlacement(rotationMatrix,                         
			      G4ThreeVector((jj+offsetLat)*(fSiPMSize+2*halfSpacing),(ii+offsetHeight)*(fSiPMSize+2*halfSpacing),(-fAcrylicLength/2+fAcrylicThickness+fTPBThickness+fSiPMThickness/2)),  
			      fLogicSiPM,              //logical volume                         
			      "SiPM",                  //name                                                                                                                         
			      logiMother,                      //mother  volume
			      false,                       //no boolean operation
			      nSiPM3++,0);                          //copy number   

	    // pos z Acrylic wall	    
	    new G4PVPlacement(rotationMatrix,                         
			      G4ThreeVector((jj+offsetLat)*(fSiPMSize+2*halfSpacing),(ii+offsetHeight)*(fSiPMSize+2*halfSpacing),(+fAcrylicLength/2-fAcrylicThickness-fTPBThickness-fSiPMThickness/2)),  
			      fLogicSiPM,              //logical volume                         
			      "SiPM",                  //name                                                                                                                         
			      logiMother,                      //mother  volume
			      false,                       //no boolean operation
			      nSiPM4++,0);                          //copy number


	  }
      }



      *rotationMatrix = G4RotationMatrix();
      rotationMatrix->rotateY(90.*deg);

    }

  else if (component=="Cathode") 
    {
      std::cout << "DetectorConstruction::DetSiPMs: Placing " << numLateral*numLength << " SiPMs in central VD cathode." << std::endl;
      *rotationMatrix = G4RotationMatrix();
      rotationMatrix->rotateZ(90.*deg);

      if (nSiPM4 > (2*numHeight*numLength+numHeight*numLength)) // if this is true its cuz we've already placed a bunch of SiPMs and should start our Cathode count at this value.
	nSiPMC = nSiPM4;

      for (int ii=-numLateral/2; ii<=numLateral/2; ii++)
      {
	for (int jj=-numLength/2; jj<=numLength/2; jj++)
	  {

	    if (!(numLateral%2) && ii==int(numLateral/2)) // if numLateral is even we need to stagger placement by 0.5 integer and skip the last one.
	      continue;
	    if (!(numLength%2) && jj==int(numLength/2)) // if numLength is even we need to stagger placement by 0.5 integer and skip the last one.
	      continue;

	    new G4PVPlacement(rotationMatrix, 
			      G4ThreeVector((ii+offsetLat)*(fSiPMSize+2*halfSpacing),0.0,(jj+offsetLength)*(fSiPMSize+2*halfSpacing)),   //cm
			      fLogicSiPM,              //logical volume                         
			      "SiPM",                  //name                           
			      logiMother,                      //mother  volume
			      false,                       //no boolean operation
			      nSiPMC++,0);                          //copy number   

	  }
      }
      
    }

  else
    {
      std::cout << "DetSiPMs: No entiendo where to poner los SiPMs."  << std::endl;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


void DetectorConstruction::DetAPEX(G4String component, G4LogicalVolume* logiMother)
{
  
  // eventually, need to be grabbed by detectorMessenger
  const double fArapucaWallgap = 40.0*cm ; 
  const double fArapucaHeight =  8.0*cm ; 
  const double fArapucaThickness = 0.6*cm;
  const double fArapucaSpacing = 4.0*cm ; 
  const double fArapucaLength  = 200.0*cm;
  const double fArapucaLengthShort  = 40.0*cm;

  const double fHbarFrontHeight = 6*cm;
  const double fHbarBackHeight = fArapucaHeight+fArapucaSpacing; // This is too high, but it will enforce light-tightness for light behind Arapucas.
  const double fHbarThickness  = 0.5*cm;
  //  std::cout << "fTargetRadius, fAcrylicLength, fAcrylicRadius, halfSpacing, SiPMSize+2*halfSpacing" << fTargetRadius << ", " << fAcrylicLength << ", " << fAcrylicRadius << ", "<< halfSpacing <<  ", " << fSiPMSize+2*halfSpacing<< std::endl;

  int numHeight = int ((2*(fTargetRadius-fArapucaWallgap)/(fArapucaHeight+fArapucaSpacing) ));
  int numLength = int ((2*(fTargetLength/2.-fArapucaWallgap)/(fArapucaLength+fArapucaSpacing) ));
  int numLengthShortwall = int ((2*(fTargetRadius-fArapucaWallgap)/(fArapucaLengthShort+fArapucaSpacing) )); // allow the shortwall to use 40cm length bars

  // force an even number of arapucas about center, for symmetry. -- EC
  if (numHeight%2) // odd
    {numHeight--;}
  if (numLength%2) // odd
    {numLength--;}
  if (numLengthShortwall%2) // odd
    {numLengthShortwall--;}


  G4Box* sArapuca = new G4Box("Arapuca", fArapucaThickness/2, fArapucaHeight/2, fArapucaLength/2);
  fLogicArapuca = new G4LogicalVolume(sArapuca,       //shape  
				   fArapucaMater,    //material  
				   "Arapuca");      //name                                                                                                                             
  G4Box* sArapuca2 = new G4Box("Arapuca", fArapucaThickness/2, fArapucaHeight/2, fArapucaLengthShort/2);
  fLogicArapuca2 = new G4LogicalVolume(sArapuca2,       //shape  
				   fArapucaMater,    //material  
				   "Arapuca");      //name

  G4Box* sHBack = new G4Box("Aluminium", fHbarThickness/2, fHbarBackHeight/2, fArapucaLength/2);
  G4Box* sHIbar = new G4Box("Aluminium", fArapucaThickness*1.5/2, fHbarThickness/2, fArapucaLength/2);
  G4Box* sHFront = new G4Box("Aluminium", fHbarThickness/2, fHbarFrontHeight/2, fArapucaLength/2);
  G4Box* sHBackEnd = new G4Box("Aluminium", fHbarThickness/2, fHbarBackHeight/2, fArapucaLengthShort/2);
  G4Box* sHIbarEnd = new G4Box("Aluminium", fArapucaThickness*1.5/2, fHbarThickness/2, fArapucaLengthShort/2);
  G4Box* sHFrontEnd = new G4Box("Aluminium", fHbarThickness/2, fHbarFrontHeight/2, fArapucaLengthShort/2);
  fLogicHBack = new G4LogicalVolume(sHBack,       //shape  
				   fAluminumMater,    //material  
				   "Hbar");      //name                                                                                                                        
  fLogicHIbar = new G4LogicalVolume(sHIbar,       //shape  
				   fAluminumMater,    //material  
				   "Hbar");      //name                                                                                                                        
  fLogicHFront = new G4LogicalVolume(sHFront,       //shape  
				   fAluminumMater,    //material  
				   "Hbar");      //name                                                                                                                        
  fLogicHBackEnd = new G4LogicalVolume(sHBackEnd,       //shape  
				   fAluminumMater,    //material  
				   "Hbar");      //name                                                                                                                        
  fLogicHIbarEnd = new G4LogicalVolume(sHIbarEnd,       //shape  
				   fAluminumMater,    //material  
				   "Hbar");      //name                                                                                                                        
  fLogicHFrontEnd = new G4LogicalVolume(sHFrontEnd,       //shape  
				   fAluminumMater,    //material  
				   "Hbar");      //name                                                                                                                        
  
  
  int nArapuca(0);
  int nArapuca2((numHeight-1)*(numLength-1));
  int nArapuca3(2*((numHeight-1)*(numLength-1)));
  int nArapuca4(2*(numHeight-1)*(numLength-1)+(numHeight-1)*(numLengthShortwall-1));
  int nArapucaC(0);
  int nhbar(0);
  
  G4RotationMatrix* rotationMatrix = new G4RotationMatrix();

  
  if (component=="Perimeter") 
    {
      std::cout << "DetectorConstruction::DetAPEX: Placing " << numHeight*numLength << " Arapucas and H-bars just inside each x cryostat wall of VD detector." << std::endl;


      for (int ii=-numHeight/2; ii<=numHeight/2; ii++)	
      {
	if (ii==0) continue;
	for (int jj=-numLength/2; jj<=numLength/2; jj++)	
	  {
	    if (jj == 0) continue;
	    if (!(numHeight%2) && ii==int(numHeight/2)) // if numLateral is even we need to stagger placement by 0.5 integer and skip the last one.
	      continue;

	    // neg x Target (inside cryo) wall
	    int signi = (ii<0) ? -1 : (ii>0) ;
	    int signj = (jj<0) ? -1 : (jj>0) ;
	    // neg x Target (inside cryo) wall
	    new G4PVPlacement(0, 
			      G4ThreeVector(-fTargetRadius+fArapucaWallgap,ii*fArapucaSpacing+signi*(fabs(ii)-0.5)*fArapucaHeight,jj*fArapucaSpacing+signj*(fabs(jj)-0.5)*fArapucaLength),   //cm
			      fLogicArapuca,              //logical volume                         
			      "Arapuca",                  //name                         
			      logiMother,                      //mother  volume
			      false,                       //no boolean operation
			      nArapuca++,0);                          //copy number   
	    // HIbar, HFront, HBack
	    new G4PVPlacement(0, 
			      G4ThreeVector(-fTargetRadius+fArapucaWallgap,fArapucaHeight/2.+fArapucaSpacing/2. +ii*fArapucaSpacing+signi*(fabs(ii)-0.5)*fArapucaHeight,jj*fArapucaSpacing+signj*(fabs(jj)-0.5)*fArapucaLength),   //cm
			      fLogicHIbar,              //logical volume                         
			      "Hbar",                  //name                         
			      logiMother,                      //mother  volume
			      false,                       //no boolean operation
			      nhbar++,0);                          //copy number
	    new G4PVPlacement(0, 
			      G4ThreeVector(-1.5*fArapucaThickness/2.-fHbarThickness/2. -fTargetRadius+fArapucaWallgap,fArapucaHeight/2.+fArapucaSpacing/2. +ii*fArapucaSpacing+signi*(fabs(ii)-0.5)*fArapucaHeight,jj*fArapucaSpacing+signj*(fabs(jj)-0.5)*fArapucaLength),   //cm
			      fLogicHBack,              //logical volume                         
			      "Hbar",                  //name                         
			      logiMother,                      //mother  volume
			      false,                       //no boolean operation
			      nhbar++,0);                          //copy number
	    new G4PVPlacement(0, 
			      G4ThreeVector(+1.5*fArapucaThickness/2.+fHbarThickness/2. -fTargetRadius+fArapucaWallgap,fArapucaHeight/2.+fArapucaSpacing/2. +ii*fArapucaSpacing+signi*(fabs(ii)-0.5)*fArapucaHeight,jj*fArapucaSpacing+signj*(fabs(jj)-0.5)*fArapucaLength),   //cm
			      fLogicHFront,              //logical volume                         
			      "Hbar",                  //name                         
			      logiMother,                      //mother  volume
			      false,                       //no boolean operation
			      nhbar++,0);                          //copy number



	    
            // pos x Target (inside cryo) wall
	    new G4PVPlacement(0, 
			      G4ThreeVector(fTargetRadius-fArapucaWallgap,ii*fArapucaSpacing+signi*(fabs(ii)-0.5)*fArapucaHeight,jj*fArapucaSpacing+signj*(fabs(jj)-0.5)*fArapucaLength),   //cm
			      fLogicArapuca,              //logical volume                         
			      "Arapuca",                  //name                             
			      logiMother,                      //mother  volume
			      false,                       //no boolean operation
			      nArapuca2++,0);                          //copy number

	    // HIbar, HFront, HBack
	    new G4PVPlacement(0, 
			      G4ThreeVector(fTargetRadius-fArapucaWallgap,fArapucaHeight/2.+fArapucaSpacing/2. +ii*fArapucaSpacing+signi*(fabs(ii)-0.5)*fArapucaHeight,jj*fArapucaSpacing+signj*(fabs(jj)-0.5)*fArapucaLength),   //cm
			      fLogicHIbar,              //logical volume                         
			      "Hbar",                  //name                         
			      logiMother,                      //mother  volume
			      false,                       //no boolean operation
			      nhbar++,0);                          //copy number
	    new G4PVPlacement(0, 
			      G4ThreeVector(+1.5*fArapucaThickness/2.+fHbarThickness/2. +fTargetRadius-fArapucaWallgap,fArapucaHeight/2.+fArapucaSpacing/2. +ii*fArapucaSpacing+signi*(fabs(ii)-0.5)*fArapucaHeight,jj*fArapucaSpacing+signj*(fabs(jj)-0.5)*fArapucaLength),   //cm
			      fLogicHBack,              //logical volume                         
			      "Hbar",                  //name                         
			      logiMother,                      //mother  volume
			      false,                       //no boolean operation
			      nhbar++,0);                          //copy number
	    new G4PVPlacement(0, 
			      G4ThreeVector(-1.5*fArapucaThickness/2.-fHbarThickness/2. +fTargetRadius-fArapucaWallgap,fArapucaHeight/2.+fArapucaSpacing/2. +ii*fArapucaSpacing+signi*(fabs(ii)-0.5)*fArapucaHeight,jj*fArapucaSpacing+signj*(fabs(jj)-0.5)*fArapucaLength),   //cm
			      fLogicHFront,              //logical volume                         
			      "Hbar",                  //name                         
			      logiMother,                      //mother  volume
			      false,                       //no boolean operation
			      nhbar++,0);                          //copy number
	    
	  }
      }

      std::cout << "DetectorConstruction::DetAPEX: Placing " << numLengthShortwall*numHeight << " Arapucas  -- with 1/5th the 200cm length -- just inside each z cryostat wall of VD detector." << std::endl;
      *rotationMatrix = G4RotationMatrix();
      rotationMatrix->rotateY(90.*deg);
      for (int ii=-numHeight/2; ii<=numHeight/2; ii++)
      {
	if (ii == 0) continue;
	for (int jj=-numLengthShortwall/2; jj<=numLengthShortwall/2; jj++)
	  {
	    if (jj == 0) continue;
 	    // neg z Target (inside cryo) wall
	    int signi = (ii<0) ? -1 : (ii>0) ;
	    int signj = (jj<0) ? -1 : (jj>0) ;
	    // neg z Target (inside cryo) wall	    
	    new G4PVPlacement(rotationMatrix,                         
			      G4ThreeVector(jj*fArapucaSpacing+signj*(fabs(jj)-0.5)*fArapucaLengthShort,ii*fArapucaSpacing+signi*(fabs(ii)-0.5)*fArapucaHeight,-fTargetLength/2+fArapucaWallgap),
			      fLogicArapuca2,              //logical volume                         
			      "Arapuca",                  //name             
			      logiMother,                      //mother  volume
			      false,                       //no boolean operation
			      nArapuca3++,0);                          //copy number   

	    // HIbar, HFront, HBack
	    new G4PVPlacement(rotationMatrix, 
			      G4ThreeVector(jj*fArapucaSpacing+signj*(fabs(jj)-0.5)*fArapucaLengthShort,fArapucaHeight/2.+fArapucaSpacing/2. +ii*fArapucaSpacing+signi*(fabs(ii)-0.5)*fArapucaHeight,-fTargetLength/2+fArapucaWallgap), //cm
			      fLogicHIbar,              //logical volume                         
			      "Hbar",                  //name                         
			      logiMother,                      //mother  volume
			      false,                       //no boolean operation
			      nhbar++,0);                          //copy number
	    new G4PVPlacement(rotationMatrix, 
			      G4ThreeVector(jj*fArapucaSpacing+signj*(fabs(jj)-0.5)*fArapucaLengthShort,fArapucaHeight/2.+fArapucaSpacing/2. +ii*fArapucaSpacing+signi*(fabs(ii)-0.5)*fArapucaHeight,-1.5*fArapucaThickness/2.-fHbarThickness/2.-fTargetLength/2+fArapucaWallgap), //cm
			      fLogicHBack,              //logical volume                         
			      "Hbar",                  //name                         
			      logiMother,                      //mother  volume
			      false,                       //no boolean operation
			      nhbar++,0);                          //copy number
	    new G4PVPlacement(rotationMatrix, 
			      G4ThreeVector(jj*fArapucaSpacing+signj*(fabs(jj)-0.5)*fArapucaLengthShort,fArapucaHeight/2.+fArapucaSpacing/2. +ii*fArapucaSpacing+signi*(fabs(ii)-0.5)*fArapucaHeight,+1.5*fArapucaThickness/2.+fHbarThickness/2.-fTargetLength/2+fArapucaWallgap), //cm
			      fLogicHFront,              //logical volume                         
			      "Hbar",                  //name                         
			      logiMother,                      //mother  volume
			      false,                       //no boolean operation
			      nhbar++,0);                          //copy number



	    // pos z Target (inside cryo) wall	    
	    new G4PVPlacement(rotationMatrix,                         
			      G4ThreeVector(jj*fArapucaSpacing+signj*(fabs(jj)-0.5)*fArapucaLengthShort,ii*fArapucaSpacing+signi*(fabs(ii)-0.5)*fArapucaHeight,fTargetLength/2-fArapucaWallgap),  
			      fLogicArapuca2,              //logical volume                         
			      "Arapuca",                  //name                     
			      logiMother,                      //mother  volume
			      false,                       //no boolean operation
			      nArapuca4++,0);                          //copy number

	    // HIbar, HFront, HBack
	    new G4PVPlacement(rotationMatrix, 
			      G4ThreeVector(jj*fArapucaSpacing+signj*(fabs(jj)-0.5)*fArapucaLengthShort,fArapucaHeight/2.+fArapucaSpacing/2. +ii*fArapucaSpacing+signi*(fabs(ii)-0.5)*fArapucaHeight,fTargetLength/2-fArapucaWallgap), //cm
			      fLogicHIbar,              //logical volume                         
			      "Hbar",                  //name                         
			      logiMother,                      //mother  volume
			      false,                       //no boolean operation
			      nhbar++,0);                          //copy number
	    new G4PVPlacement(rotationMatrix, 
			      G4ThreeVector(jj*fArapucaSpacing+signj*(fabs(jj)-0.5)*fArapucaLengthShort,fArapucaHeight/2.+fArapucaSpacing/2. +ii*fArapucaSpacing+signi*(fabs(ii)-0.5)*fArapucaHeight,+1.5*fArapucaThickness/2.+fHbarThickness/2.+fTargetLength/2-fArapucaWallgap), //cm
			      fLogicHBack,              //logical volume                         
			      "Hbar",                  //name                         
			      logiMother,                      //mother  volume
			      false,                       //no boolean operation
			      nhbar++,0);                          //copy number
	    new G4PVPlacement(rotationMatrix, 
			      G4ThreeVector(jj*fArapucaSpacing+signj*(fabs(jj)-0.5)*fArapucaLengthShort,fArapucaHeight/2.+fArapucaSpacing/2. +ii*fArapucaSpacing+signi*(fabs(ii)-0.5)*fArapucaHeight,-1.5*fArapucaThickness/2.-fHbarThickness/2.+fTargetLength/2-fArapucaWallgap), //cm
			      fLogicHFront,              //logical volume                         
			      "Hbar",                  //name                         
			      logiMother,                      //mother  volume
			      false,                       //no boolean operation
			      nhbar++,0);                          //copy number

	  }
      }



      *rotationMatrix = G4RotationMatrix();
      rotationMatrix->rotateY(90.*deg);

    }

  else if (component=="Cathode") 
    {
      std::cout << "DetectorConstruction::DetAPEX: Placing  Arapucas in central VD cathode. Code Not built: QUITTING." << std::endl;
      exit(0);
      *rotationMatrix = G4RotationMatrix();
      rotationMatrix->rotateZ(90.*deg);

      /*
      if (nSiPM4 > (2*numHeight*numLength+numHeight*numLength)) // if this is true its cuz we've already placed a bunch of SiPMs and should start our Cathode count at this value.
	nSiPMC = nSiPM4;

      for (int ii=-numLateral/2; ii<=numLateral/2; ii++)
      {
	for (int jj=-numLength/2; jj<=numLength/2; jj++)
	  {

	    if (!(numLateral%2) && ii==int(numLateral/2)) // if numLateral is even we need to stagger placement by 0.5 integer and skip the last one.
	      continue;
	    if (!(numLength%2) && jj==int(numLength/2)) // if numLength is even we need to stagger placement by 0.5 integer and skip the last one.
	      continue;

	    new G4PVPlacement(rotationMatrix, 
			      G4ThreeVector((ii+offsetLat)*(fSiPMSize+2*halfSpacing),0.0,(jj+offsetLength)*(fSiPMSize+2*halfSpacing)),   //cm
			      fLogicSiPM,              //logical volume                         
			      "SiPM",                  //name                           
			      logiMother,                      //mother  volume
			      false,                       //no boolean operation
			      nSiPMC++,0);                          //copy number   

	  }
      }
      */
    }

  else
    {
      std::cout << "DetAPEX: No entiendo where to poner los bars/arapucas."  << std::endl;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........OAoooOO0OOooo......
