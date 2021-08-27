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
#include "G4UserLimits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction()
:G4VUserDetectorConstruction(),
 fTargetMater(0), fLogicTarget(0), 
 fWorldMater(0), fPhysiWorld(0),
 fDetectorMessenger(0)
{
  fInsetRadius = 0.0*cm;

  fWorldLength = fTubeLength + 1.0*m;
  fWorldRadius = fTubeRadius + 1.0*m;

  fTubeThickness = 4.0*mm;

  fShieldThickness = (455.5-431.8)*mm;


  fFieldCageInnerRadius = 86.1/2.*mm;
  fFieldCageOuterRadius = fFieldCageInnerRadius+5.0*fTubeThickness;
  fFieldCageLength = (290.2-129.9)*mm;

  fBufferInnerRadius = fFieldCageInnerRadius;
  fBufferOuterRadius = fFieldCageOuterRadius;
  fBufferLength = (431.8-366.0 -1 )*mm;

  fLensRadius = 20.0*mm;
  fElectrodeThickness = 3*fTubeThickness;
  fElectrodeInnerRadius = 85.0/2. *mm;
  fElectrodeOuterRadius = fElectrodeInnerRadius + 3.*fTubeThickness;
  
  fBaffleInnerRadius = fFieldCageOuterRadius;
  fBaffleOuterRadius = fBaffleInnerRadius + fTubeThickness/2.;
  fBaffleEndGap = 5.0*mm;

  DefineMaterials();
    
  fDetectorMessenger = new DetectorMessenger(this);

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
  G4Element* Cu  = new G4Element("Copper", "Cu", 29, 63.5*g/mole);
  G4Element* Cr = new G4Element("Chromium" ,"Cr" , 24., 51.99*g/mole);
  G4Element* Fe = new G4Element("Iron"     ,"Fe" , 26., 55.84*g/mole);
  G4Element* Ni = new G4Element("Nickel"   ,"Ni" , 28., 58.69*g/mole);

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

  G4Material* g10 =  new G4Material("NemaG10", 1.700*g/cm3, ncomponents=4); // Nema Arkani-Hamed G10?
  g10->AddElement(Si, number_of_atoms=1);
  g10->AddElement(O , number_of_atoms=2);
  g10->AddElement(C , number_of_atoms=3);
  g10->AddElement(H , number_of_atoms=3);

  G4Material* ssteel = 
    new G4Material("Stainless-Steel", 8*g/cm3, ncomponents=3);
  ssteel->AddElement(Fe, 74*perCent);
  ssteel->AddElement(Cr, 18*perCent);
  ssteel->AddElement(Ni,  8*perCent);

  G4Material* copper = 
    new G4Material("Copper", 8.96*g/cm3, ncomponents=1);
  copper->AddElement(Cu, 100*perCent);


  //
  fWorldMater = Air20;
  //  fGapMater = Air20;


  fTubeMater = ssteel;
  fShieldMater = ssteel;

  fBufferMater = ssteel;
  fELPMater = ssteel;
  fELPPMater = ssteel;
  fECMater = ssteel;
  fFieldCageMater = ssteel;
  fBaffleMater = ssteel;
  fBaffleFlangeMater = ssteel;

  fG10Mater = g10;

  /*

    Target is a Box nested inside World box. Detector is a Box inset from Target walls, nested Target, as tall (Length) as Target.
    EC, 20-Oct-2019.

   */

}


G4VPhysicalVolume* DetectorConstruction::ConstructVolumes()
{


  // This G4Material requires having read up the TargetPressure from config file, which hasn't been done by  DefineMaterials() above.
  G4double gasFactor;
  gasFactor = 5.894; // 1 bar: gm/L 
  gasFactor /= 1E3; // gm/cm^3
  gasFactor *= fTargetPressure / atmosphere ; // presuming wrongly ideal gas law

  fTargetMater   = new G4Material("Xenon", 54, 131.3*g/mole, gasFactor *g/cm3, kStateGas,  77.*kelvin, fTargetPressure/atmosphere);

  /*
  G4Tubs*
  sWorld = new G4Tubs("World",                                 //name
                 0.,fWorldRadius, 0.5*fWorldLength, 0.,twopi); //dimensions  
  */
  fWorldRadius = 5*fTubeRadius;
  fWorldLength = 5*fTubeLength;

  G4Tubs* sWorld = new G4Tubs("World", 0.0, fWorldRadius, fWorldLength/2., 0., twopi);
               
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
                            0, true);                         //copy number


  // Tube
  //
  

  G4Tubs*  sTube = new G4Tubs("Tube",                                   //name
                  fTubeRadius, fTubeRadius+fTubeThickness, 0.5*fTubeLength, 0.,twopi); //dimensions
    fLogicTube = new G4LogicalVolume(sTube,           //shape
                             fTubeMater,              //material
                             "Tube");                 //name
                               
           new G4PVPlacement(0,                         //no rotation
			     G4ThreeVector(0.,0.,0.),             //at (0,0,0)
                           fLogicTube,                //logical volume
                           "Tube",                    //name
                           lWorld,                      //mother  volume
                           false,                       //no boolean operation
                           0, true);                          //copy number

  // Tube Endcaps
  G4Tubs*  sTubeECN = new G4Tubs("TubeECN",                                   //name
                  fTubeRadius, fTubeRadius+4.*fTubeThickness, fTubeThickness, 0.,twopi); //dimensions
    fLogicTubeECN = new G4LogicalVolume(sTubeECN,           //shape
                             fTubeMater,              //material
                             "TubeECN");                 //name
                               
           new G4PVPlacement(0,                         //no rotation
			     G4ThreeVector(0.,0.,fTubeLength/2.-fTubeThickness),             //at (0,0,0)
                           fLogicTubeECN,                //logical volume
                           "TubeECN",                    //name
                           lWorld,                      //mother  volume
                           false,                       //no boolean operation
                           0, true);                          //copy number

  G4Tubs*  sTubeECS = new G4Tubs("TubeECS",                                   //name
                  fTubeRadius, fTubeRadius+2.*fTubeThickness, fTubeThickness, 0.,twopi); //dimensions
    fLogicTubeECS = new G4LogicalVolume(sTubeECS,           //shape
                             fTubeMater,              //material
                             "TubeECS");                 //name
                               
           new G4PVPlacement(0,                         //no rotation
			     G4ThreeVector(0.,0.,-fTubeLength/2.+fBaffleEndGap+fTubeThickness),             //at (0,0,0)
                           fLogicTubeECS,                //logical volume
                           "TubeECS",                    //name
                           lWorld,                      //mother  volume
                           false,                       //no boolean operation
                           0, true);                          //copy number


  G4Tubs*  sFlangeECN = new G4Tubs("FlangeECN",                                   //name
				   fTubeRadius+2.*fTubeThickness,fTubeRadius+2.*fTubeThickness+2.*fTubeThickness, 2.*fTubeThickness, 0.,twopi); //dimensions
    fLogicFlangeECN = new G4LogicalVolume(sFlangeECN,           //shape
                             fTubeMater,              //material
                             "FlangeECN");                 //name
                               
           new G4PVPlacement(0,                         //no rotation
			     G4ThreeVector(0.,0.,fTubeLength/2.),
                           fLogicFlangeECN,                //logical volume
                           "FlangeECN",                    //name
                           lWorld,                      //mother  volume
                           false,                       //no boolean operation
                           0, true);                          //copy number

  G4Tubs*  sFlangeECS = new G4Tubs("FlangeECS",                                   //name
				   fTubeRadius+2.*fTubeThickness,fTubeRadius+2.*fTubeThickness+2.*fTubeThickness, 2.*fTubeThickness, 0.,twopi); //dimensions
    fLogicFlangeECS = new G4LogicalVolume(sFlangeECS,           //shape
                             fTubeMater,              //material
                             "FlangeECS");                 //name
                               
           new G4PVPlacement(0,                         //no rotation
			     G4ThreeVector(0.,0.,-fTubeLength/2.),             //at (0,0,0)
                           fLogicFlangeECS,                //logical volume
                           "FlangeECS",                    //name
                           lWorld,                      //mother  volume
                           false,                       //no boolean operation
                           0, true);                          //copy number


  // shield Endcaps
  G4Tubs*  sShieldECN = new G4Tubs("ShieldECN",                                   //name
                  0.0, fTubeRadius+4.*fTubeThickness, fShieldThickness/2., 0.,twopi); //dimensions
    fLogicShieldECN = new G4LogicalVolume(sShieldECN,           //shape
                             fTubeMater,              //material
                             "ShieldECN");                 //name
                               
           new G4PVPlacement(0,                         //no rotation
			     G4ThreeVector(0.,0.,fTubeLength/2.+4.*fTubeThickness+fShieldThickness/2.),             //at (0,0,0)
                           fLogicShieldECN,                //logical volume
                           "ShieldECN",                    //name
                           lWorld,                      //mother  volume
                           false,                       //no boolean operation
                           0, true);                          //copy number

  G4Tubs*  sShieldECS = new G4Tubs("ShieldECS",                                   //name
		  fLensRadius,fTubeRadius+4.*fTubeThickness, fShieldThickness/2., 0.,twopi); //dimensions
    fLogicShieldECS = new G4LogicalVolume(sShieldECS,           //shape
                             fTubeMater,              //material
                             "ShieldECS");                 //name

           new G4PVPlacement(0,                         //no rotation
			     G4ThreeVector(0.,0.,-fTubeLength/2.-4.*fTubeThickness-fShieldThickness/2.),             //at (0,0,0)
                           fLogicShieldECS,                //logical volume
                           "ShieldECS",                    //name
                           lWorld,                      //mother  volume
                           false,                       //no boolean operation
                           0, true);                          //copy number


  G4Tubs* sTarget = new G4Tubs("Target", 0.0, fTubeRadius, fTubeLength/2.+4.*fTubeThickness,0.0,twopi);
  fLogicTarget = new G4LogicalVolume(sTarget,           //shape
                             fTargetMater,              //material
                             "Target");                 //name
                               
           new G4PVPlacement(0,                         //no rotation
			     G4ThreeVector(0.,0.,0.),             //at (0,0,0)
                           fLogicTarget,                //logical volume
                           "Target",                    //name
                           lWorld,                      //mother  volume
                           false,                       //no boolean operation
                           0, true);                          //copy number



  // Above is all the outer diameter stuff. Now we move to the parts inside the GXe

  fTargetLength = fTubeLength + (2.+1.)*fTubeThickness;
  //  fTargetCenter = G4ThreeVector(0.0,0.0,fTargetLength/2.);

  G4Tubs* sBuffer = new G4Tubs("Buffer", fBufferInnerRadius, fBufferOuterRadius, fBufferLength/2.,0.0,twopi);
  fLogicBuffer = new G4LogicalVolume(sBuffer,           //shape
				     fBufferMater,              //material
				     "Buffer");                 //name
                               
          new G4PVPlacement(0,                         //no rotation
			     G4ThreeVector(0.,0.,366.0+fBufferLength/2.-fTubeLength/2.),             //at (0,0,0)
                           fLogicBuffer,                //logical volume
                           "Buffer",                    //name
                           fLogicTarget,                      //mother  volume
                           false,                       //no boolean operation
                           0, true);                          //copy number


  G4Tubs* sFieldCage = new G4Tubs("FieldCage", fFieldCageInnerRadius, fFieldCageOuterRadius, fFieldCageLength/2.,0.0,twopi);
  fLogicFieldCage = new G4LogicalVolume(sFieldCage,           //shape
				     fFieldCageMater,              //material
				     "FieldCage");  

  G4double zMidELP = 326.1*mm;
  G4double zMidELPP = 331.8*mm;
  G4double zMidFC = 290.2-fTubeLength/2.-fFieldCageLength/2.; // z coord of center of Xenon volume: could be a negative number
          new G4PVPlacement(0,                         //no rotation
			     G4ThreeVector(0.,0.,zMidFC),             //at (0,0,0)
                           fLogicFieldCage,                //logical volume
                           "FieldCage",                    //name
                           fLogicTarget,                      //mother  volume
                           false,                       //no boolean operation
                           0, true);                          //copy number

  G4Tubs* sELP = new G4Tubs("ELP", fElectrodeInnerRadius, fElectrodeOuterRadius, fElectrodeThickness,0.0,twopi);
         fLogicELP = new G4LogicalVolume(sELP,           //shape
				     fELPMater,              //material
				     "ELP");  

          new G4PVPlacement(0,                         //no rotation
			     G4ThreeVector(0.,0.,zMidELP-fTubeLength/2.-fElectrodeThickness/2.),             //at (0,0,0)
                           fLogicELP,                //logical volume
                           "ELP",                    //name
                           fLogicTarget,                      //mother  volume
                           false,                       //no boolean operation
                           0, true);                          //copy number


  G4Tubs* sELPP = new G4Tubs("ELPP", fElectrodeInnerRadius, fElectrodeOuterRadius, fElectrodeThickness,0.0,twopi);
         fLogicELPP = new G4LogicalVolume(sELPP,           //shape
				     fELPPMater,              //material
				     "ELPP");  

          new G4PVPlacement(0,                         //no rotation
			     G4ThreeVector(0.,0.,zMidELPP-fTubeLength/2. + fElectrodeThickness/2.),             //at (0,0,0)
                           fLogicELPP,                //logical volume
                           "ELPP",                    //name
                           fLogicTarget,                      //mother  volume
                           false,                       //no boolean operation
                           0, true);                          //copy number


  G4Tubs* sEC = new G4Tubs("EC", fElectrodeInnerRadius, fElectrodeOuterRadius, fElectrodeThickness,0.0,twopi);
  G4double zMidEC = 129.9*mm;
         fLogicEC = new G4LogicalVolume(sEC,           //shape
				     fECMater,              //material
				     "EC");  

          new G4PVPlacement(0,                         //no rotation
			     G4ThreeVector(0.,0.,zMidEC-fTubeLength/2. - fElectrodeThickness/2.),    
                           fLogicEC,                //logical volume
                           "EC",                    //name
                           fLogicTarget,                      //mother  volume
                           false,                       //no boolean operation
                           0, true);                          //copy number


  G4double zMidBaffle = (109.0*mm-fBaffleEndGap)/2.+fBaffleEndGap;
  G4Tubs* sBaffle = new G4Tubs("Baffle", fBaffleInnerRadius, fBaffleOuterRadius, zMidBaffle,0.0,twopi);
         fLogicBaffle = new G4LogicalVolume(sBaffle,           //shape
				     fBaffleMater,              //material
				     "Baffle");  

          new G4PVPlacement(0,                         //no rotation
			     G4ThreeVector(0.,0.,-fTubeLength/2.+zMidBaffle+fBaffleEndGap),   
                           fLogicBaffle,                //logical volume
                           "Baffle",                    //name
                           fLogicTarget,                      //mother  volume
                           false,                       //no boolean operation
                           0, true);                          //copy number


  G4Tubs* sBaffleFlange = new G4Tubs("BaffleFlange", fBaffleOuterRadius, fBaffleOuterRadius+4.0*fTubeThickness, fTubeThickness/2.,0.0,twopi);
         fLogicBaffleFlange = new G4LogicalVolume(sBaffleFlange,           //shape
				     fBaffleMater,              //material
				     "BaffleFlange");  

          new G4PVPlacement(0,                         //no rotation
			     G4ThreeVector(0.,0.,fBaffleEndGap+fTubeThickness/2./2.-fTubeLength/2.),
                           fLogicBaffleFlange,                //logical volume
                           "BaffleFlange",                    //name
                           fLogicTarget,                      //mother  volume
                           false,                       //no boolean operation
                           0, true);                          //copy number








	   G4UserLimits *lim = new G4UserLimits();
	   G4double maxStep = 1.0 * mm;
	   lim->SetMaxAllowedStep(maxStep);
	   fLogicTarget->SetUserLimits(lim);




  
  std::cout << "fTubeRadius is " << fTubeRadius << std::endl;
  std::cout << "fTubeThickness is " << fTubeThickness << std::endl;
  std::cout << "fTubeLength is " << fTubeLength << std::endl;
  std::cout << "fShieldThickness is " << fShieldThickness << std::endl;




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
  G4cout << "\n Tube : Thickness = " << G4BestUnit(fTubeThickness,"Length")
         << " Material = " << fTubeMater->GetName() << G4endl;          
  G4cout << "\n Shield : Thickness = " << G4BestUnit(fShieldThickness,"Length");
    //         << " Material = " << fShieldMater->GetName() << G4endl;          
  
  G4cout << "\n Target : Length = " << G4BestUnit(fTargetLength,"Length")
         << " Radius = " << G4BestUnit(fTargetRadius,"Length")  
	 << " Material = " << fTargetMater->GetName()
         << " Pressure = " << G4BestUnit(fTargetPressure,"Pressure")  ;
  
  
  G4cout << "\n Flange : Length = " << G4BestUnit(fFlangeThickness,"Length")
         << " Radius = " << G4BestUnit(fFlangeThickness,"Length")  ;
    //         << " Material = " << fFlangeMater->GetName();
  
  //  G4cout << "\n G10 : Thickness = " << G4BestUnit(fG10Thickness,"Length")
  //     << " Material = " << fG10Mater->GetName() << G4endl;          

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
    G4RunManager::GetRunManager()->PhysicsHasBeenModified();
  } else {
    G4cout << "\n--> warning from DetectorConstruction::SetDetectorMaterial : "
           << materialChoice << " not found" << G4endl;
  }              
}




//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetInsetRadius(G4double value)
{
  fInsetRadius = value;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}


void DetectorConstruction::SetTargetPressure(G4double value)
{
  fTargetPressure = value;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}

void DetectorConstruction::SetFlangeThickness(G4double value)
{
  fFlangeThickness = value;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}

void DetectorConstruction::SetFlangeMaterial(G4String value)
{

  G4Material* pttoMaterial =
     G4NistManager::Instance()->FindOrBuildMaterial(value);   

  G4RunManager::GetRunManager()->ReinitializeGeometry();
  fFlangeMater = pttoMaterial;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


void DetectorConstruction::SetTubeLength(G4double value)
{
  fTubeLength = value;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}
void DetectorConstruction::SetTubeRadius(G4double value)
{
  fTubeRadius = value;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}
void DetectorConstruction::SetTubeThickness(G4double value)
{
  fTubeThickness = value;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}

void DetectorConstruction::SetShieldThickness(G4double value)
{
  fShieldThickness = value;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}
void DetectorConstruction::SetShieldRadius(G4double value)
{
  fShieldRadius = value;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double DetectorConstruction::GetTargetLength()
{
  return fTargetLength;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double DetectorConstruction::GetTargetRadius()
{
  return fTargetRadius;
}
G4double DetectorConstruction::GetShieldThickness()
{
  return fShieldThickness;
}

G4double DetectorConstruction::GetTargetPressure()
{
  return fTargetPressure;
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

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4LogicalVolume* DetectorConstruction::GetLogicDetector()
{
  return fLogicDetector;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......



