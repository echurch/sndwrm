// Module 3 light simulation
// Authors: L. Paulucci & F. Marinho & E. Church
// Date:  2024
//
// Added modifications should be reported to the original authors for updating authorship

#include "DetectorConstruction.hh"
#include "MaterialPropertyLoader.hh"
#include "G4LogicalVolumeStore.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4Tubs.hh"
#include "G4EllipticalTube.hh"
#include "G4Orb.hh"
#include "G4Sphere.hh"
#include "G4NistManager.hh"
#include "G4MultiUnion.hh"
#include "G4Color.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"

#include <string>

#include "globals.hh"
#include <CLHEP/Geometry/Transform3D.h>


DetectorConstruction::DetectorConstruction()
  
  :fDefaultMaterial(NULL),
   fPhysiWorld(NULL),fLogicWorld(NULL),fSolidWorld(NULL),
   fPhysiVol(NULL),fLogicVol(NULL),fSolidVol(NULL)
  
{
  fWorldSizeX=50.0; //in meters
  fWorldSizeY=50.0; //in meters
  fWorldSizeZ=75.0;
  
  fCryostat_x = 14.8; //in meters
  //fCryostat_y = 6.5; //in meters
  fCryostat_y = 13.0; //in meters
  fCryostat_z = 60.0; //in meters
  fFC_x = 13.5; //in meters
  fFC_y = fCryostat_y; //in meters
  fFC_z = fCryostat_z; //in meters
  
  fCathode_x = 13.5; //Cathode size in m
  fCathode_z = 60.0; //Cathode size in m

  fLatY = 6.5; //APA internal size in cm
  fLatZ = 60.0;

  fthickness=0.10; //m

  fFCOut_x=0.0345;//m
  fFCOut_y=0.04; //m
  fFCOut_z= 3.0/2; //m

  fAra_x = 0.007; //Arapuca window size in m  
  fAra_y = 0.50; //Arapuca window size in m
  fAra_z = 0.50; //Arapuca window size in m

  fAra_offset = 0.07; //sensor offset wrt FC profiles in m
  fAras_yspacing = 0.005; //vertical spacing between XARAPUCAS in m
  fptp_width = 2.0*um;
  //fptp_width = 10*cm;

  fvert_bar_x = 0.075; //in m
  fvert_bar_y = fFC_y/2-0.1;//making it slightly smaller to avoid overlap on cathode
  fvert_bar_z = 0.075;

  // IBeams
  fIFlangeWidth = 0.402; // all m here.
  fIFlangeThick = 0.040;
  fIFlangeWaist = 0.022;
  fIFlangeHeight = 1108./1000.;
  fITopLength = 17832./1000. + 2* (fIFlangeHeight/2.0);
  fISideLength = 16732./1000. - 2* (fIFlangeHeight/2.0) ; // need a little space with these side beams

  fIPortHoleRad = 0.800/2;
  fISidePortLoc = 5907./1000. - fIFlangeHeight/2. ;
  fIPortSpacing = 4000.0/1000.0 ;
  fIBotPortLoc = 5000.0/1000.0;

  fht = 16732/1000./2 + 0.045; // m, for beam and belt placement
  fst = 17832/1000./2 + 0.030;
  fzpl = 64732./1000.;
  fSpacing = 64732./1000./41; // m
  
  fMPL = new MaterialPropertyLoader();  
}

  
  

DetectorConstruction::~DetectorConstruction()
{
  delete fDefaultMaterial;

}

G4VPhysicalVolume* DetectorConstruction::Construct()
  
{

  // Now for the purpose of tracking optical photons we do the following to the Argon and TPB to endow them w optical physics properties.
  // Must wait till this late, cuz MLP works by looping over all Logical Volumes which are only just now established.
  // Get the logical volume store and assign material properties. MaterialPropLoader() is borrowed, heavily-edited from LArSoft.   

  DefineMaterials();
  return ConstructLine();
}


// https://indico.cern.ch/event/698002/contributions/2868259/attachments/1591642/2519098/AC_introduction_to_the_cryostat_design_warm_vessel.pdf
void DetectorConstruction::IBeams()
{

  G4double Offset(0.001*m);
  G4Box* IBeamTopFlange = new G4Box("IBeamTopFlange",(fIFlangeWidth/2.0)*m, (fIFlangeThick/2.0)*m,(fITopLength/2.0)*m); 
  G4Box* IBeamTopMid = new G4Box("IBeamTopMid",(fIFlangeWaist/2.0)*m, (fIFlangeHeight/2.)*m,(fITopLength/2.)*m);
  G4Box* IBeamSideFlange = new G4Box("IBeamTopFlange",(fIFlangeWidth/2.0)*m, (fIFlangeThick/2.0)*m,(fISideLength/2.0)*m); 
  G4Box* IBeamSideMidtmp0 = new G4Box("IBeamSideMid",(fIFlangeWaist/2.0)*m, (fIFlangeHeight/2.)*m,(fISideLength/2.)*m); 
  G4Tubs* IBeamPort = new G4Tubs("IBeamPortHole",0.,fIPortHoleRad*m,(fIFlangeThick/2.0)*m,0.0,2.0*CLHEP::pi);
  G4RotationMatrix* fc = new G4RotationMatrix();
  G4ThreeVector* axisfc = new G4ThreeVector(0.0,1.0,0.0);
  fc->rotate(CLHEP::pi/2.,axisfc);
  G4RotationMatrix* fc2 = new G4RotationMatrix();
  G4ThreeVector* axisfc2 = new G4ThreeVector(1.0,0.0,0.0);
  G4ThreeVector* axisfc3 = new G4ThreeVector(0.0,0.0,1.0);
  fc2->rotate(-CLHEP::pi/2,axisfc2);
  fc2->rotate(CLHEP::pi/2,axisfc3);
  G4RotationMatrix* fc3 = new G4RotationMatrix();
  fc3->rotate(-CLHEP::pi/2,axisfc2);
  
  G4SubtractionSolid* IBeamBotMidtmp = new G4SubtractionSolid("IBeamBottomtmp", IBeamTopMid, IBeamPort,fc, G4ThreeVector(0.0,0.0,fIPortSpacing/2.*m));
  G4SubtractionSolid* IBeamBotMid = new G4SubtractionSolid("IBeamBottom", IBeamBotMidtmp, IBeamPort,fc, G4ThreeVector(0.0,0.0,-fIPortSpacing/2.*m));
  G4SubtractionSolid* IBeamSideMidtmp1 = new G4SubtractionSolid("IBeamSidetmp", IBeamSideMidtmp0, IBeamPort,fc, G4ThreeVector(0.0,0.0,(fISideLength/2.+fIFlangeHeight/2.-5907./1000.)*m));
  G4SubtractionSolid* IBeamSideMidtmp2 = new G4SubtractionSolid("IBeamSidetmp2", IBeamSideMidtmp1, IBeamPort,fc, G4ThreeVector(0.0,0.0,(fISideLength/2.+fIFlangeHeight/2.-5907./1000.-fIPortSpacing)*m));
  G4SubtractionSolid* IBeamSideMid = new G4SubtractionSolid("IBeamSide", IBeamSideMidtmp2, IBeamPort,fc, G4ThreeVector(0.0,0.0,(fISideLength/2.+fIFlangeHeight/2.-5907./1000.-2.0*fIPortSpacing)*m));

  
  HepGeom::Transform3D tnull, tr1, tr2;
  tnull = HepGeom::TranslateY3D(0.0);
  tr1 = HepGeom::TranslateY3D(( fIFlangeHeight/2.+fIFlangeThick/2.0)*m);
  tr2 = HepGeom::TranslateY3D((-fIFlangeHeight/2.-fIFlangeThick/2.0)*m);

  G4MultiUnion* fBeamTopVol = new G4MultiUnion("IBeamTopVol");
  G4MultiUnion* fBeamBotVol = new G4MultiUnion("IBeamBotVol");
  G4MultiUnion* fBeamSideVol = new G4MultiUnion("IBeamSideVol");

  fBeamTopVol->AddNode(IBeamTopMid, tnull);
  fBeamTopVol->AddNode(IBeamTopFlange,tr1);
  fBeamTopVol->AddNode(IBeamTopFlange,tr2);
  fBeamTopVol->Voxelize();
  G4LogicalVolume* fIBeamTopLog = new G4LogicalVolume(fBeamTopVol, fDUNESteel, "IBeamTopLog");
  fBeamBotVol->AddNode(IBeamBotMid, tnull);
  fBeamBotVol->AddNode(IBeamTopFlange,tr1);
  fBeamBotVol->AddNode(IBeamTopFlange,tr2);
  fBeamBotVol->Voxelize();
  G4LogicalVolume* fIBeamBotLog = new G4LogicalVolume(fBeamBotVol, fDUNESteel, "IBeamBotLog");
  fBeamSideVol->AddNode(IBeamSideMid, tnull);
  fBeamSideVol->AddNode(IBeamSideFlange,tr1);
  fBeamSideVol->AddNode(IBeamSideFlange,tr2);
  fBeamSideVol->Voxelize();
  G4LogicalVolume* fIBeamSideLog = new G4LogicalVolume(fBeamSideVol, fDUNESteel, "IBeamSideLog");

  G4VisAttributes* simpleBoxAtt= new G4VisAttributes(G4Colour::Green());
  simpleBoxAtt->SetDaughtersInvisible(true);
  simpleBoxAtt->SetForceSolid(true);
  simpleBoxAtt->SetForceAuxEdgeVisible(true);
  fIBeamTopLog->SetVisAttributes(simpleBoxAtt);
  fIBeamBotLog->SetVisAttributes(simpleBoxAtt);
  fIBeamSideLog->SetVisAttributes(simpleBoxAtt);

  const double ht =  fht;  //m  
  const double st =  fst;  //m

  const double zbsp = fSpacing; //m
    
  //  Top, Bottom, sides
  int cpIT(0), cpIB(0), cpIL(0), cpIR(0);
  double zpl(0.0);
  double mIL (0);
  
  for (size_t ii=0;ii<=19;ii++) {
    
     new G4PVPlacement(fc,G4ThreeVector(0,(ht)*m,(zpl)*m),"IBeamTop",
							 fIBeamTopLog,      //its logical volume   
							 fPhysOuterAir,           //its mother  volume
							 false,                 //no boolean operation
							 cpIT++, // copyNo
							 true); //check for overlaps
     new G4PVPlacement(fc,G4ThreeVector(0,(-ht)*m,(zpl)*m),"IBeamBot",
							 fIBeamBotLog,      //its logical volume   
							 fPhysOuterAir,           //its mother  volume
							 false,                 //no boolean operation
							 cpIB++, // copyNo
							 true); //check for overlaps
     new G4PVPlacement(fc2,G4ThreeVector((-st)*m,0,(zpl)*m),"IBeamLeft",
							 fIBeamSideLog,      //its logical volume   
							 fPhysOuterAir,           //its mother  volume
							 false,                 //no boolean operation
							 cpIL++, // copyNo
							 true); //check for overlaps
     mIL += fIBeamSideLog->GetMass()/kg;
     new G4PVPlacement(fc2,G4ThreeVector((+st)*m,0,(zpl)*m),"IBeamRight",
							 fIBeamSideLog,      //its logical volume   
							 fPhysOuterAir,           //its mother  volume
							 false,                 //no boolean operation
							 cpIR++, // copyNo
							 true); //check for overlaps

    
    if (ii==0) {  zpl+=zbsp; continue; }
    new G4PVPlacement(fc,G4ThreeVector(0,(ht)*m,(-zpl)*m),"IBeamTop",
							 fIBeamTopLog,      //its logical volume   
							 fPhysOuterAir,           //its mother  volume
							 false,                 //no boolean operation
							 cpIT++, // copyNo
							 true); //check for overlaps
    new G4PVPlacement(fc,G4ThreeVector(0,(-ht)*m,(-zpl)*m),"IBeamBot",
							 fIBeamBotLog,      //its logical volume   
							 fPhysOuterAir,           //its mother  volume
							 false,                 //no boolean operation
							 cpIB++, // copyNo
							 true); //check for overlaps
    new G4PVPlacement(fc2,G4ThreeVector((-st)*m,0,(-zpl)*m),"IBeamLeft",
							 fIBeamSideLog,      //its logical volume   
							 fPhysOuterAir,           //its mother  volume
							 false,                 //no boolean operation
							 cpIL++, // copyNo
							 true); //check for overlaps
    new G4PVPlacement(fc2,G4ThreeVector((+st)*m,0,(-zpl)*m),"IBeamRight",
							 fIBeamSideLog,      //its logical volume   
							 fPhysOuterAir,           //its mother  volume
							 false,                 //no boolean operation
							 cpIR++, // copyNo
							 true); //check for overlaps

    zpl+=zbsp;
  }



  // Front face, back face
  int cpIF(0), cpIBk(0);
  double xpl(0.0);
  zpl = fzpl/2.; //m
  for (size_t ii=0;ii<=4;ii++) {
    // use zpl where it's finished at large +-ive value.
     new G4PVPlacement(fc3,G4ThreeVector((+xpl)*m,0,(zpl)*m),"IBeamFront",
							 fIBeamSideLog,      //its logical volume   
							 fPhysOuterAir,           //its mother  volume
							 false,                 //no boolean operation
							 cpIF++, // copyNo
							 true); //check for overlaps
     new G4PVPlacement(fc3,G4ThreeVector((+xpl)*m,0,(-zpl)*m),"IBeamBack",
							 fIBeamSideLog,      //its logical volume   
							 fPhysOuterAir,           //its mother  volume
							 false,                 //no boolean operation
							 cpIBk++, // copyNo
							 true); //check for overlaps

    if (ii==0) {  xpl+=zbsp; continue; }
     new G4PVPlacement(fc3,G4ThreeVector((-xpl)*m,0,(+zpl)*m),"IBeamFront",
							 fIBeamSideLog,      //its logical volume   
							 fPhysOuterAir,           //its mother  volume
							 false,                 //no boolean operation
							 cpIF++, // copyNo
							 true); //check for overlaps
     new G4PVPlacement(fc3,G4ThreeVector((-xpl)*m,0,(-zpl)*m),"IBeamBack",
							 fIBeamSideLog,      //its logical volume   
							 fPhysOuterAir,           //its mother  volume
							 false,                 //no boolean operation
							 cpIBk++, // copyNo
							 true); //check for overlaps
     xpl+=zbsp;
  }

  std::cout << "DetectorConstruction::IBeams(): mass of Left IBeams [kg] " << mIL << std::endl;
  
}

void DetectorConstruction::Belts()
{
  // Just two belts, one with a hole, one without.

  const double ht =  fht ; //m
  const double st =  fst ; //m

  G4Box* BeltFlange = new G4Box("BeltFlange", ((0.200)/2.0)*m, (fIFlangeWaist/2.0)*m, (fSpacing/2.-fIFlangeWaist/2.-0.001)*m ); 
  G4Box* BeltMid = new G4Box("BeltMid",(fIFlangeWaist/2.0)*m, (fIFlangeHeight/4.)*m, (fSpacing/2.-fIFlangeWaist/2.-0.001)*m);
  G4Tubs* IBeamPort = new G4Tubs("BeltPortHole",0.,0.25*m,(fIFlangeThick/2.0)*m,0.0,2.0*CLHEP::pi); // 0.6m??
  G4RotationMatrix* fc = new G4RotationMatrix();
  G4RotationMatrix* fc3 = new G4RotationMatrix();
  G4ThreeVector* axisfc = new G4ThreeVector(0.0,0.0,1.0);
  G4RotationMatrix* fc2 = new G4RotationMatrix();
  G4ThreeVector* axisfc2 = new G4ThreeVector(0.0,1.0,0.0);
  G4ThreeVector* axisfc3 = new G4ThreeVector(1.0,0.0,0.0);
  fc->rotate(CLHEP::pi/2.,axisfc);
  fc2->rotate(CLHEP::pi/2.,axisfc2);
  fc3->rotate(CLHEP::pi/2.,axisfc);
  fc3->rotate(CLHEP::pi/2.,axisfc3);
  G4SubtractionSolid* BeltHole = new G4SubtractionSolid("BeltHole", BeltMid, IBeamPort, fc2, G4ThreeVector(0.0,0.0,0.0));

  HepGeom::Transform3D tnull, tr1, tr2;
  tnull = HepGeom::TranslateY3D(0.0);
  tr1 = HepGeom::TranslateY3D( (fIFlangeHeight/4.0)*m);
  tr2 = HepGeom::TranslateY3D(-(fIFlangeHeight/4.0)*m);

  G4MultiUnion* BeltHoleUni = new G4MultiUnion("BeltHoleUni");
  BeltHoleUni->AddNode(BeltHole, tnull);
  BeltHoleUni->AddNode(BeltFlange,tr1);
  BeltHoleUni->AddNode(BeltFlange,tr2);
  BeltHoleUni->Voxelize();
  G4LogicalVolume* BeltHoleUniLog = new G4LogicalVolume(BeltHoleUni, fDUNESteel, "BeltHoleUniLog");
  G4MultiUnion* BeltUni = new G4MultiUnion("BeltUni");
  BeltUni->AddNode(BeltMid, tnull);
  BeltUni->AddNode(BeltFlange,tr1);
  BeltUni->AddNode(BeltFlange,tr2);
  BeltUni->Voxelize();
  G4LogicalVolume* BeltUniLog = new G4LogicalVolume(BeltUni, fDUNESteel, "BeltUniLog");
  G4VisAttributes* simpleBoxAtt= new G4VisAttributes(G4Colour::Cyan());
  G4VisAttributes* simpleBoxHoleAtt= new G4VisAttributes(G4Colour::Brown());
  simpleBoxAtt->SetDaughtersInvisible(true);
  simpleBoxAtt->SetForceSolid(true);
  simpleBoxAtt->SetForceAuxEdgeVisible(true);
  simpleBoxHoleAtt->SetDaughtersInvisible(true);
  simpleBoxHoleAtt->SetForceSolid(true);
  simpleBoxHoleAtt->SetForceAuxEdgeVisible(true);
  BeltUniLog->SetVisAttributes(simpleBoxAtt);
  BeltHoleUniLog->SetVisAttributes(simpleBoxHoleAtt);

  
  const double zbsp = fSpacing; //m
  double zpl(zbsp/2.);
  double xpl(0.);
  double eps(0.215);
  int cpIT(0), cpIB(0), cpIL(0), cpIR(0),cpBlt(0);  
  // Top, bottom

 
  for (size_t ii=0;ii<=19 /*20*/;ii++) {
    // loop on x for top and bottom
    for (int jj=-5;jj<5;jj++) {

      new G4PVPlacement(0,G4ThreeVector((jj+0.5)*zbsp*m,(-ht+eps)*m,(zpl)*m),"BeltBot",
							 BeltHoleUniLog,      //its logical volume   
							 fPhysOuterAir,           //its mother  volume
							 false,                 //no boolean operation
							 cpIB++, // copyNo
							 true); //check for overlaps
      new G4PVPlacement(0,G4ThreeVector((jj+0.5)*zbsp*m,(-ht+eps)*m,(-zpl)*m),"BeltBot",
							 BeltHoleUniLog,      //its logical volume   
							 fPhysOuterAir,           //its mother  volume
							 false,                 //no boolean operation
							 cpIB++, // copyNo
							 true); //check for overlaps
      // belts dodge the flanges on top
      if (std::abs(jj)==1 || std::abs(jj)==2 || std::abs(jj)==4) { 
	new G4PVPlacement(0,G4ThreeVector((jj+0.5)*zbsp*m,(ht-eps)*m,(zpl)*m),"BeltTop",
							 BeltUniLog,      //its logical volume   
							 fPhysOuterAir,           //its mother  volume
							 false,                 //no boolean operation
							 cpIT++, // copyNo
							 true); //check for overlaps
	new G4PVPlacement(0,G4ThreeVector((jj+0.5)*zbsp*m,(ht-eps)*m,(-zpl)*m),"BeltTop",
							 BeltUniLog,      //its logical volume   
							 fPhysOuterAir,           //its mother  volume
							 false,                 //no boolean operation
							 cpIT++, // copyNo
							 true); //check for overlaps
      }

    }
    // left and right sides. bot 3 y-levels just below the port hole heights, top one at the top
    for (size_t jj=0;jj<4;jj++) {
      double y;
      G4LogicalVolume* belt = BeltHoleUniLog;
     // loop on y bot to top for these  two side walls -- 2 of 3 is nohole
      if (jj==3) { // top 
	y = (ht-eps)*m ;
	belt = BeltUniLog;
	if (ii % 3)  belt = BeltHoleUniLog;
      }
      if (jj==2) { // bot hole
	y = -(fISideLength/2.+fIFlangeHeight/2.-5907./1000.+0.0*fIPortSpacing + 9*fIPortHoleRad)/2.*m;
      }
      if (jj==1) { // up 1 hole
	y = -(fISideLength/2.+fIFlangeHeight/2.-5907./1000.-2.0*fIPortSpacing + 9*fIPortHoleRad)/2.*m;
	belt = BeltUniLog;
	if ((ii+1) % 3)  belt = BeltHoleUniLog;
      }
      if (jj==0) { // top hole
	y = -(fISideLength/2.+fIFlangeHeight/2.-5907./1000.-4.0*fIPortSpacing + 9*fIPortHoleRad)/2.*m;
	belt = BeltUniLog;
	if ((ii+2) % 3)  belt = BeltHoleUniLog;
      }

      if (ii==20) { // nohole on extreme ends for all four levels
	belt = BeltUniLog;
      }

      new G4PVPlacement(fc,G4ThreeVector(-st*m,y,(-zpl)*m),"BeltLeft",
							 belt,      //its logical volume   
							 fPhysOuterAir,           //its mother  volume
							 false,                 //no boolean operation
							 cpBlt++, // copyNo
							 true); //check for overlaps
      new G4PVPlacement(fc,G4ThreeVector( st*m,y,(-zpl)*m),"BeltRight",
							 belt,      //its logical volume   
							 fPhysOuterAir,           //its mother  volume
							 false,                 //no boolean operation
							 cpBlt++, // copyNo
							 true); //check for overlaps
      new G4PVPlacement(fc,G4ThreeVector(-st*m,y,(+zpl)*m),"BeltLeft",
							 belt,      //its logical volume   
							 fPhysOuterAir,           //its mother  volume
							 false,                 //no boolean operation
							 cpBlt++, // copyNo
							 true); //check for overlaps
      new G4PVPlacement(fc,G4ThreeVector( st*m,y,(+zpl)*m),"BeltRight",
							 belt,      //its logical volume   
							 fPhysOuterAir,           //its mother  volume
							 false,                 //no boolean operation
							 cpBlt++, // copyNo
							 true); //check for overlaps
      
    }
    zpl+=zbsp;
  }
  
  
  
  // Front face, back face
  int cpBF(0), cpBBk(0);
  xpl = zbsp/2.; //m
  zpl = fzpl/2.+0.100; //m
  for (size_t ii=0;ii<5;ii++) {
  for (size_t jj=0;jj<4;jj++) {
      double y;
      G4LogicalVolume* belt = BeltHoleUniLog;
     // loop on y bot to top for these  two side walls -- 2 of 3 is nohole
      if (jj==3) { // top 
	y = ht*m;
	belt = BeltUniLog;
	if (ii % 3)  belt = BeltHoleUniLog;
      }
      if (jj==2) { // bot hole
	y = -(fISideLength/2.+fIFlangeHeight/2.-5907./1000.+0.0*fIPortSpacing + 9*fIPortHoleRad)/2.*m;
      }
      if (jj==1) { // up 1 hole
	y = -(fISideLength/2.+fIFlangeHeight/2.-5907./1000.-2.0*fIPortSpacing + 9*fIPortHoleRad)/2.*m;
	belt = BeltUniLog;
	if ((ii+1) % 3)  belt = BeltHoleUniLog;
      }
      if (jj==0) { // top hole
	y = -(fISideLength/2.+fIFlangeHeight/2.-5907./1000.-4.0*fIPortSpacing + 9*fIPortHoleRad)/2.*m;
	belt = BeltUniLog;
	if ((ii+2) % 3)  belt = BeltHoleUniLog;
      }

      if (ii==20) { // nohole on extreme ends for all four levels
	belt = BeltUniLog;
      }

      new G4PVPlacement(fc3,G4ThreeVector(-xpl*m,y,(-zpl)*m),"BeltBack",
							 belt,      //its logical volume   
							 fPhysOuterAir,           //its mother  volume
							 false,                 //no boolean operation
							 cpBBk++, // copyNo
							 true); //check for overlaps
      new G4PVPlacement(fc3,G4ThreeVector( xpl*m,y,(-zpl)*m),"BeltBack",
							 belt,      //its logical volume   
							 fPhysOuterAir,           //its mother  volume
							 false,                 //no boolean operation
							 cpBBk++, // copyNo
							 true); //check for overlaps
      new G4PVPlacement(fc3,G4ThreeVector(-xpl*m,y,(+zpl)*m),"BeltFront",
							 belt,      //its logical volume   
							 fPhysOuterAir,           //its mother  volume
							 false,                 //no boolean operation
							 cpBF++, // copyNo
							 true); //check for overlaps
      new G4PVPlacement(fc3,G4ThreeVector( xpl*m,y,(+zpl)*m),"BeltFront",
							 belt,      //its logical volume   
							 fPhysOuterAir,           //its mother  volume
							 false,                 //no boolean operation
							 cpBF++, // copyNo
							 true); //check for overlaps
      
    }
    xpl+=zbsp;
  }



}

void DetectorConstruction::ShieldingWalls()
{

  double eps(0.215);
  const double ht =  fht ; //m
  const double st =  fst + 0.031 ; //m have to be out from the belts left, right

  // These are the y-locations of centers of shield panels, not the actual holes.
  const double yTopHole = -(fISideLength/2.+fIFlangeHeight/2.-5907./1000.-(4.0-1.)*fIPortSpacing + 9*fIPortHoleRad)/2.;
  const double yBotHole = -(fISideLength/2.+fIFlangeHeight/2.-5907./1000.-(0.0-1.)*fIPortSpacing + 9*fIPortHoleRad)/2.;
  const double y1HoleUp = -(fISideLength/2.+fIFlangeHeight/2.-5907./1000.-(2.0-1.)*fIPortSpacing + 9*fIPortHoleRad)/2.;;
  const double yTop = (ht-eps);
  
  const double BlockWidth = fSpacing-fIFlangeWaist-0.001 - 0.050;
  const double BlockHeight = fIPortSpacing - 0.050;
  const double BlockHeightTop =  BlockHeight*0.63 ;
  std::cout << "BlockHeight,Top are " << BlockHeight << ", " << BlockHeightTop << std::endl;
  const double BlockHeightBot = BlockHeight*0.15;

  double fBlockThickness(0.40);
  G4Box* ShieldBlock = new G4Box("ShieldBlock",(fBlockThickness/2.0)*m, (BlockHeight/2.)*m, (BlockWidth/2.)*m); // 20cm for the fBPSE density 1.60, 30cm for fBP density 1.0
  G4Box* ShieldBlockTop = new G4Box("ShieldBlockTop",(fBlockThickness/2.0)*m, (BlockHeightTop/2.)*m, (BlockWidth/2.)*m); // 20cm for the fBPSE density 1.60, 32cm for fBP density 1.0
  G4Box* ShieldBlockBot = new G4Box("ShieldBlockBot",(fBlockThickness/2.0)*m, (BlockHeightBot/2.)*m, (BlockWidth/2.)*m); // 20cm for the fBPSE density 1.60, 32cm for fBP density 1.0

  G4RotationMatrix* fc = new G4RotationMatrix();
  G4RotationMatrix* fc3 = new G4RotationMatrix();
  G4ThreeVector* axisfc = new G4ThreeVector(0.0,0.0,1.0);
  G4RotationMatrix* fc2 = new G4RotationMatrix();
  G4ThreeVector* axisfc2 = new G4ThreeVector(0.0,1.0,0.0);
  G4ThreeVector* axisfc3 = new G4ThreeVector(1.0,0.0,0.0);
  fc->rotate(CLHEP::pi/2.,axisfc);
  fc2->rotate(CLHEP::pi/2.,axisfc2);
  fc3->rotate(CLHEP::pi/2.,axisfc);
  fc3->rotate(CLHEP::pi/2.,axisfc3);

  G4LogicalVolume* ShieldBlockLog = new G4LogicalVolume(ShieldBlock, fBP, "ShieldBlockLog");
  G4VisAttributes* simpleShieldAtt= new G4VisAttributes(G4Colour::Blue());
  simpleShieldAtt->SetDaughtersInvisible(true);
  simpleShieldAtt->SetForceSolid(true);
  simpleShieldAtt->SetForceAuxEdgeVisible(true);
  ShieldBlockLog->SetVisAttributes(simpleShieldAtt);

  G4LogicalVolume* ShieldBlockTopLog = new G4LogicalVolume(ShieldBlockTop, fBP, "ShieldBlockTopLog");
  ShieldBlockTopLog->SetVisAttributes(simpleShieldAtt);
  G4LogicalVolume* ShieldBlockBotLog = new G4LogicalVolume(ShieldBlockBot, fBP, "ShieldBlockBotLog");
  ShieldBlockBotLog->SetVisAttributes(simpleShieldAtt);

  const double zbsp = fSpacing; //m
  double zpl(zbsp/2.);
  double xpl(0.);

  int cpIT(0), cpIB(0), cpIL(0), cpIR(0),cpBlt(0);  

  for (size_t ii=0;ii<=19;ii++) {

    // left and right sides. bot 3 y-levels just below the port hole heights, top one at the top
       for (size_t jj=0;jj<5;jj++) {
      double y;
      G4LogicalVolume* shield = ShieldBlockLog ;
     // loop on y bot to top for these  two side walls -- 2 of 3 is nohole
      if (jj==3) { // top 
	y = yTopHole + BlockHeight/2. + BlockHeightTop/2. + 0.050;
	shield = ShieldBlockTopLog;

      }
      if (jj==2) { // bot hole
	y = yBotHole;

      }
      if (jj==1) { // up 1 hole
	y = y1HoleUp;
      }
      if (jj==0) { // top hole
	y = yTopHole;
      }
      if (jj==4) { // bottom 
	y = yBotHole - BlockHeight/2. - 0.5*BlockHeightBot - 0.050 ;
	shield = ShieldBlockBotLog;

      }

      if (jj==0 || jj == 3) continue; // EC, 2-May-2025, drop upper shields.
      
      new G4PVPlacement(0,G4ThreeVector(-st*m,y*m,(-zpl)*m),"ShieldLeft",
							 shield,      //its logical volume   
							 fPhysOuterAir,           //its mother  volume
							 false,                 //no boolean operation
							 cpBlt++, // copyNo
							 true); //check for overlaps
      new G4PVPlacement(0,G4ThreeVector( st*m,y*m,(-zpl)*m),"ShieldRight",
							 shield,      //its logical volume   
							 fPhysOuterAir,           //its mother  volume
							 false,                 //no boolean operation
							 cpBlt++, // copyNo
							 true); //check for overlaps
      new G4PVPlacement(0,G4ThreeVector(-st*m,y*m,(+zpl)*m),"ShieldLeft",
							 shield,      //its logical volume   
							 fPhysOuterAir,           //its mother  volume
							 false,                 //no boolean operation
							 cpBlt++, // copyNo
							 true); //check for overlaps
      new G4PVPlacement(0,G4ThreeVector( st*m,y*m,(+zpl)*m),"ShieldRight",
							 shield,      //its logical volume   
							 fPhysOuterAir,           //its mother  volume
							 false,                 //no boolean operation
							 cpBlt++, // copyNo
							 true); //check for overlaps
      
    }
    zpl+=zbsp;
  }


  // Front face, back face
  int cpBF(0), cpBBk(0);
  xpl = zbsp/2.; //m
  zpl = fzpl/2.; //m
  for (size_t ii=0;ii<5;ii++) {

       for (size_t jj=0;jj<5;jj++) {
      double y;
      G4LogicalVolume* shield = ShieldBlockLog ;
     // loop on y bot to top for these  two side walls -- 2 of 3 is nohole
      if (jj==3) { // top 
	y = yTopHole + BlockHeight/2. + BlockHeightTop/2. + 0.050;
	shield = ShieldBlockTopLog;

      }
      if (jj==2) { // bot hole
	y = yBotHole;

      }
      if (jj==1) { // up 1 hole
	y = y1HoleUp;
      }
      if (jj==0) { // top hole
	y = yTopHole;
      }
      if (jj==4) { // bottom 
	y = yBotHole - BlockHeight/2. - 0.5*BlockHeightBot - 0.050 ;
	shield = ShieldBlockBotLog;
      }

      new G4PVPlacement(fc2,G4ThreeVector(-xpl*m,y*m,(-zpl)*m),"ShieldBack",
							 shield,      //its logical volume   
							 fPhysOuterAir,           //its mother  volume
							 false,                 //no boolean operation
							 cpBBk++, // copyNo
							 true); //check for overlaps
      new G4PVPlacement(fc2,G4ThreeVector( xpl*m,y*m,(-zpl)*m),"ShieldBack",
							 shield,      //its logical volume   
							 fPhysOuterAir,           //its mother  volume
							 false,                 //no boolean operation
							 cpBBk++, // copyNo
							 true); //check for overlaps
      new G4PVPlacement(fc2,G4ThreeVector(-xpl*m,y*m,(+zpl)*m),"ShieldFront",
							 shield,      //its logical volume   
							 fPhysOuterAir,           //its mother  volume
							 false,                 //no boolean operation
							 cpBF++, // copyNo
							 true); //check for overlaps
      new G4PVPlacement(fc2,G4ThreeVector( xpl*m,y*m,(+zpl)*m),"ShieldFront",
							 shield,      //its logical volume   
							 fPhysOuterAir,           //its mother  volume
							 false,                 //no boolean operation
							 cpBF++, // copyNo
							 true); //check for overlaps
      
    }
    
    xpl+=zbsp;
  }
  
}




void DetectorConstruction::DefineMaterials()
{  

  G4String name, symbol;             
  G4double density;            
  
  G4int natoms,nel;
  G4double z;
  
  // Define Elements   
  G4Element*   H  = new G4Element ("Hydrogen","H",1.,1.01*g/mole);
  G4Element*   C = new G4Element ("Carbon","C",6.,12.01*g/mole);
  G4Element*   O = new G4Element ("Oxygen","O",8.,16.0*g/mole);
  G4Element*  Al = new G4Element(name="Aluminium",symbol="Al",z=13.,26.98*g/mole);
  G4Element*  Fe = new G4Element(name="Iron",symbol="Fe",z=26.,55.85*g/mole);
  G4Element*  Ni = new G4Element(name="Niquel",symbol="Ni",z=28.,58.6934*g/mole);
  G4Element*  Si = new G4Element(name="Silicon",symbol="Si",z=14.,28.085*g/mole);
  G4Element*  Cr = new G4Element(name="Chromium",symbol="Cr",z=24.,51.9961*g/mole);
  G4Element* N  = new G4Element("Nitrogen", "N", 7, 14.01*g/mole);
  G4Element* Mn  = new G4Element("Manganese","Mn", 25, 54.94*g/mole);
  G4Element* Cu  = new G4Element("Copper","Cu", 29, 63.55*g/mole);
  G4Element* B10 = new G4Element (name="Boron10",symbol="B10",z=4.,10.00*g/mole);
  G4Element* B11 = new G4Element (name="Boron11",symbol="B11",z=4.,11.00*g/mole);

  G4NistManager * man = G4NistManager::Instance();

  G4Material* StainlessSteel = new G4Material(name="StainlessSteel",7.93*g/cm3,nel=4);//STEEL_STAINLESS_Fe7Cr2Ni
  StainlessSteel->AddElement(C, 0.0010);
  StainlessSteel->AddElement(Cr, 0.1792);
  StainlessSteel->AddElement(Fe, 0.7298);
  StainlessSteel->AddElement(Ni, 0.0900);

  G4int ncomponents;
  G4double fractionmass;
  fDUNESteel=new G4Material("DUNESteel"/*SS407L*/,7.93*g/cm3, ncomponents=7);
  fDUNESteel->AddElement(Fe,fractionmass=95.8/100.);
  fDUNESteel->AddElement(Mn,fractionmass=1.8/100.);
  fDUNESteel->AddElement(Ni,fractionmass=0.8/100.);
  fDUNESteel->AddElement(Si,fractionmass=0.6/100.);
  fDUNESteel->AddElement(Cu,fractionmass=0.5/100.);
  fDUNESteel->AddElement(Cr,fractionmass=0.3/100.);
  fDUNESteel->AddElement(C, fractionmass=0.2/100.);
  
  G4Material* G10 = new G4Material(name="G10",1.7*g/cm3,nel=4);
  G10->AddElement(Si, 0.2805);
  G10->AddElement(O, 0.3954);
  G10->AddElement(C, 0.2990);
  G10->AddElement(H, 0.0251);

  G4Material* Aluminium = new G4Material(name="Aluminium",z=13.,26.98*g/mole,2.7*g/cm3);
  G4Material* base_mat = man->FindOrBuildMaterial("G4_TEFLON");
  G4Material* env_mat = man->FindOrBuildMaterial("G4_lAr");
  G4Material* mAir = man->FindOrBuildMaterial("G4_AIR");

  G4Material* ptp_mat =  new G4Material(name = "ptp_mat", 1.079*g/cm3, nel = 2); //p-Terphenyl
  ptp_mat->AddElement (C, natoms=18);
  ptp_mat->AddElement (H, natoms=14);

  G4Material* Mylar =  new G4Material("Mylar", density= 1.40*g/cm3, ncomponents=3);
  Mylar->AddElement(H, natoms=4);
  Mylar->AddElement(C, natoms=5);
  Mylar->AddElement(O, natoms=2);

  //Foam. From https://indico.fnal.gov/event/20144/session/19/contribution/267/material/slides/1.pdf 

  G4int number_of_atoms;
  G4Material *foam=new G4Material("Foam",0.09*g/cm3, ncomponents=4);
  foam->AddElement(C,number_of_atoms=54);
  foam->AddElement(O,number_of_atoms=15);
  foam->AddElement(N,number_of_atoms=4);
  foam->AddElement(H,number_of_atoms=60);

  G4Material *wood=new G4Material("Wood",0.5*g/cm3, ncomponents=3);
  wood->AddElement(C,number_of_atoms=50);
  wood->AddElement(O,number_of_atoms=44);
  wood->AddElement(H,number_of_atoms=6);

  fShieldMater =  foam ;
  fWoodMater =  wood ;


  G4Material* H2O = new G4Material("Water",density= 1.0*g/cm3,ncomponents=2);
  H2O->AddElement(H, number_of_atoms=2);
  H2O->AddElement(O, number_of_atoms=1);
  // Borated-Poly SE self-extinguishing, https://johncaunt.com/products/jc207-hd-hd5/, EC, 7-Jan-2025.
  G4Material* fBP_SE = new G4Material (name="BP_SE", density= 1.60*g/cm3 /*0.95*/, ncomponents=5);
  /*
  "Borated Polyethyl."         5     0.95       #. 5%-borated polyethylene (BPE)
          "Hydrogen"               11.6       #. C.R.Wuest  SSCL-GEM TN 92-172.
          "Carbon"                 61.2       #.  composition of "Reactor Experiments, Inc."
          "Bor 11"                  4.0       #. Weight fraction               GMIX
          "Bor 10"                  1.0       #. Weight fraction               GMIX
          "Oxygen"                 22.2       #. Weight fraction               GMIX
  */
  fBP_SE->AddElement(H,fractionmass=6.6*perCent);
  fBP_SE->AddElement(C,fractionmass=66.5*perCent);
  fBP_SE->AddElement(B11,fractionmass=3.76*perCent);
  fBP_SE->AddElement(B10,fractionmass=0.94*perCent);
  fBP_SE->AddElement(O,fractionmass=22.2*perCent);

  G4Material* fBP_norm = new G4Material (name="BP_norm", density= 0.95*g/cm3 , ncomponents=5);  
  fBP_norm->AddElement(H,fractionmass=11.6*perCent);
  fBP_norm->AddElement(C,fractionmass=61.2*perCent);
  fBP_norm->AddElement(B11,fractionmass=4.0*perCent);
  fBP_norm->AddElement(B10,fractionmass=1.0*perCent);
  fBP_norm->AddElement(O,fractionmass=22.2*perCent);

  
  fBP = H2O; // H2O; //fBP_norm; //fBP_SE;
  
  /*const G4int nEntries = 6;
  G4double PhotonEnergy[nEntries] =
    { 2.0*eV, 2.341*eV, 2.757*eV, 3.353*eV, 4.136*eV, 10.0*eV };

  G4double l_lAr[nEntries] =
  {20*m,20*m,20*m,20*m,20*m,20*m}; */

    const G4int nEntries = 8;
  //_________ RELEVANT ENERGY VALUES Xe 175nm -> 7.08eV; Ar 128 -> 9.69eV_________
  G4double PhotonEnergy[nEntries] = {2.5*eV, 5.0*eV, 7.0*eV, 7.5*eV, 8.0*eV, 9.0*eV, 9.5*eV, 10.136*eV};
  G4double l_lAr[nEntries] =        {80*m,     80*m,   80*m,   80*m,   20*m,   20*m,   20*m,      20*m}; 
  //  G4double l_lAr[nEntries] =        {80*m,     80*m,   80*m,   80*m,   2000*m,   2000*m,   2000*m,      2000*m}; 


  //  G4double n_lAr[nEntries] = {1.3,1.3,1.3,1.3,1.3};
  G4double n_lAr[72] = {1.2310665394518083, 1.2311200318147684, 1.2312727482978456, 1.2313496348295065, 1.2313506914097514, 1.2314023012909403, 1.2321538166889967, 1.231909947460871, 1.2321384938953646, 1.2323923169803301, 1.232443926861519, 1.2328668458034384, 1.233095392237932, 1.233053736213583, 1.233509772502325, 1.2335361057330418, 1.2334866090123424, 1.2338668153496684, 1.2343481282888824, 1.2344502914710154, 1.2346282846045646, 1.2350337675923626, 1.2354898038811046, 1.2358874461725522, 1.2361238333033957, 1.2364787629902494, 1.2371038951359457, 1.2374666655191495, 1.237720488604115, 1.2384292913975776, 1.2388853276863196, 1.2393413639750617, 1.2402271033218288, 1.2409106294648191, 1.2417458155106422, 1.2425810015564651, 1.2435678475051206, 1.2446558000556642, 1.245844859208096, 1.2472614082147766, 1.2488296171242892, 1.2505242092861624, 1.2525221212537003, 1.2548739063278473, 1.257579564508603, 1.2607654790483278, 1.2644822032479663, 1.269108886864599, 1.2747719131505861, 1.2816987719601767, 1.2910690403154002, 1.302368759656658, 1.316713517166515, 1.3329411664124329, 1.3493804523854165, 1.365926455672306, 1.3820018039345536, 1.4010549636291998, 1.423580491960391, 1.442600693158671, 1.46103487273757, 1.4894895758980782, 1.508978836972031, 1.5300372520419365, 1.5496431051408859, 1.5737334642006038, 1.5996151986767673, 1.6187496758472137, 1.6361783106957755, 1.6540241603412935, 1.672490509546729, 1.69156768392026};

  G4double Energy_n_lar[72] = {1.88901692613609*eV, 1.915491549763*eV, 1.94434453495496*eV, 1.9740800365654*eV, 2.00473917327901*eV, 2.03636565851099*eV, 2.0690060083582*eV, 2.1027097698752*eV, 2.1375297720313*eV, 2.17352240202191*eV, 2.21074790997382*eV, 2.24927074550798*eV, 2.28915993011439*eV, 2.33048946986518*eV, 2.37333881365753*eV, 2.41779336295621*eV, 2.46394503991691*eV, 2.51189292184303*eV, 2.56174395119099*eV, 2.61361373183203*eV, 2.66762742404865*eV, 2.72392075285046*eV, 2.78264114670945*eV, 2.84394902682875*eV, 2.90801927068415*eV, 2.97504287795441*eV, 3.04522887226339*eV, 3.11880647861692*eV, 3.19602762431763*eV, 3.27716982084588*eV, 3.36253949617556*eV, 3.45247586185779*eV, 3.54735541774646*eV, 3.64759722049586*eV, 3.75366907130329*eV, 3.86609481562225*eV, 3.98546299517591*eV, 4.11243715385428*eV, 4.24776817847434*eV, 4.39230915909235*eV, 4.54703339014745*eV, 4.71305631519082*eV, 4.89166246132417*eV, 5.08433873911506*eV, 5.29281593505014*eV, 5.51912084854449*eV, 5.76564240171605*eV, 6.03521629502331*eV, 6.33123457629323*eV, 6.65778911807759*eV, 7.01986191167459*eV, 7.42358102512033*eV, 7.80768693883725*eV, 8.18166617653237*eV, 8.49578904787588*eV, 8.75179837892604*eV, 8.97145190403002*eV, 9.21112941189563*eV, 9.40734062594689*eV, 9.56190983757659*eV, 9.71915769371575*eV, 9.88023992289941*eV, 9.99798584772925*eV, 10.0930619775564*eV, 10.1644798618311*eV, 10.2587258039568*eV, 10.3321720936895*eV, 10.3923239488175*eV, 10.4560982471211*eV, 10.4827908882082*eV, 10.5288181909777*eV, 10.5783779318939*eV};

  G4double ray_e_lAr[21] = { 1.18626*eV, 1.68626*eV, 2.18626*eV, 2.68626*eV, 3.18626*eV, 3.68626*eV, 4.18626*eV, 4.68626*eV, 5.18626*eV, 5.68626*eV, 6.18626*eV, 6.68626*eV, 7.18626*eV, 7.68626*eV, 8.18626*eV, 8.68626*eV, 9.18626*eV, 9.68626*eV, 10.1863*eV, 10.6863*eV, 11.1863*eV};
  G4double ray_s_lAr[21] = { 1200800*cm, 390747*cm, 128633*cm, 54969.1*cm, 27191.8*cm, 14853.7*cm, 8716.9*cm, 5397.42*cm, 3481.37*cm, 2316.51*cm, 1577.63*cm, 1092.02*cm, 763.045*cm, 534.232*cm, 371.335*cm, 252.942*cm, 165.38*cm, 99.9003*cm, 51.2653*cm, 17.495*cm, 0.964341*cm };


  G4double RayleighEnergies[22] = {   2.80*eV,   3.00*eV,   3.50*eV,   4.00*eV,  5.00*eV,  6.00*eV,  7.00*eV,  8.00*eV,  8.50*eV,  9.00*eV,  9.20*eV,  9.40*eV,  9.50*eV,  9.60*eV,  9.70*eV,  9.80*eV,  9.90*eV,  10.0*eV,  10.2*eV,  10.4*eV,  10.6*eV, 10.8*eV };
    G4double RayleighSpectrum[22] = { 47923.*cm, 35981.*cm, 18825.*cm, 10653.*cm, 3972.*cm, 1681.*cm, 750.9*cm, 334.7*cm, 216.8*cm, 135.0*cm, 109.7*cm, 88.06*cm, 78.32*cm, 69.34*cm, 61.06*cm, 53.46*cm, 46.50*cm, 40.13*cm, 28.91*cm, 19.81*cm, 12.61*cm, 7.20*cm };

  G4MaterialPropertiesTable* lAr_pt = new G4MaterialPropertiesTable();
  lAr_pt->AddProperty("RINDEX", Energy_n_lar, n_lAr, 72);
  //lAr_pt->AddProperty("RAYLEIGH", ray_e_lAr, ray_s_lAr, 21);

  lAr_pt->AddProperty("RAYLEIGH", RayleighEnergies, RayleighSpectrum, 22);
  lAr_pt->AddProperty("ABSLENGTH", PhotonEnergy, l_lAr, nEntries);
  lAr_pt->AddConstProperty("SCINTILLATIONTIMECONSTANT1", 6. * ns);       // FASTTIMECONSTANT
  lAr_pt->AddConstProperty("SCINTILLATIONTIMECONSTANT2", 1590. * ns);//	  SLOWTIMECONSTANT

  std::vector<double> FastScintEnergies { 6.0*eV,  6.7*eV,  7.1*eV,  7.4*eV,  7.7*eV, 7.9*eV,  8.1*eV,  8.4*eV,  8.5*eV,  8.6*eV,  8.8*eV,  9.0*eV,  9.1*eV,  9.4*eV,  9.8*eV,  10.4*eV,  10.7*eV};
  std::vector<double> SlowScintEnergies { 6.0*eV,  6.7*eV,  7.1*eV,  7.4*eV,  7.7*eV, 7.9*eV,  8.1*eV,  8.4*eV,  8.5*eV,  8.6*eV,  8.8*eV,  9.0*eV,  9.1*eV,  9.4*eV,  9.8*eV,  10.4*eV,  10.7*eV};
  //  std::vector<double> FastScintSpectrumloc { 0.0,  0.04, 0.12, 0.27, 0.44, 0.62, 0.80, 0.91, 0.92, 0.85, 0.70, 0.50, 0.31, 0.13, 0.04,  0.01, 0.0};
  // std::vector<double> SlowScintSpectrumloc { 0.0,  0.04, 0.12, 0.27, 0.44, 0.62, 0.80, 0.91, 0.92, 0.85, 0.70, 0.50, 0.31, 0.13, 0.04,  0.01, 0.0};
  std::vector<double> FastScintSpectrumloc { 0.0,  0.0, 1., 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,  0.0, 0.0};
  std::vector<double> SlowScintSpectrumloc { 0.0,  0.0, 1., 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,  0.0, 0.0};
  lAr_pt->AddProperty("SCINTILLATIONCOMPONENT1", FastScintEnergies, FastScintSpectrumloc); // FASTCOMPNENT
  lAr_pt->AddProperty("SCINTILLATIONCOMPONENT2", SlowScintEnergies, SlowScintSpectrumloc); // SLOWCOMPONENT
  lAr_pt->AddConstProperty("SCINTILLATIONYIELD", 25000 / MeV );
  lAr_pt->AddConstProperty("SCINTILLATIONYIELD1", 0.3 ); // YIELDRATIO
  lAr_pt->AddConstProperty("RESOLUTIONSCALE", 1.0 );
  env_mat->GetIonisation()->SetBirksConstant(0.069 * cm / MeV);


  // By commenting out below, and running these lines instead I will enforce the G4_lAr properties set in MaterialPropertiesLoader.
  env_mat->SetMaterialPropertiesTable(lAr_pt); 
  std::cout << "Dumping G4_lAr properties .... " << std::endl;
  env_mat->GetMaterialPropertiesTable()->DumpTable();
  //fMPL->SetPropertiesFromServices();  // fills local LArprop class with hard-coded data cutnpasted from fcl file.
  //fMPL->GetPropertiesFromServices();  // Shoves these into local MaterialTables
  //fMPL->UpdateGeometry(G4LogicalVolumeStore::GetInstance()); // Finally, loads properties into G4MaterialProperties
  
  const G4int num = 3;
  G4double n_teflon[num] = {1.41, 1.41, 1.41};
  G4double l_teflon[num] = {3448*m,  4082*m, 9174*m};
  G4double pp[num] = {1.9587*eV, 7.5407*eV, 10.1328*eV};

  G4MaterialPropertiesTable*  Teflon_pt = new G4MaterialPropertiesTable();
  Teflon_pt->AddProperty("RINDEX", pp, n_teflon, num);
  Teflon_pt->AddProperty("ABSLENGTH", pp, l_teflon, num);

  base_mat->SetMaterialPropertiesTable(Teflon_pt);

  G4Material* acrylic = new G4Material("acrylic", 1.18*g/cm3,3); //acrylic
  acrylic->AddElement (C, 5);
  acrylic->AddElement (O, 2);
  acrylic->AddElement (H, 8);

  //__________________ PTP parametrisation __________________//
  
  const int ptp_e_entries = 22;
  const int ptp_a_entries = 15;
  const int ptp_fe_entries = 9;
  //  const int ptp_ab_entries = 3;

  // Set emission energy range and parameters

  G4double ptp_e_energy[ptp_e_entries] = 
    {3.06434*eV, 3.1098*eV,3.14086*eV, 3.18325*eV, 3.24903*eV, 
    3.31176*eV, 3.33517*eV, 3.37092*eV,3.4198*eV, 3.45108*eV, 
    3.47012*eV, 3.48937*eV, 3.51537*eV, 3.53513*eV, 3.55512*eV, 
    3.59577*eV, 3.62339*eV, 3.63736*eV, 3.65144*eV, 3.66562*eV, 
    3.67992*eV, 3.7235*eV};

  G4double ptp_e[ptp_e_entries] = 
    {0.016285, 0.0520565, 0.0788851, 0.132542, 0.275628, 
     0.409771, 0.543914, 0.678057, 0.722772, 0.821143, 
     0.955286, 1, 0.946343, 0.8122, 0.740657, 0.8122, 
     0.678057, 0.534971, 0.400828, 0.266685, 0.132542, 0.00734215};

  // Set absorption energy range and parameters
  G4double ptp_a_energy[ptp_a_entries] =
    {2.0*eV, 3.0*eV, 4.0*eV, 4.5*eV, 5.0*eV, 5.5*eV, 6.0*eV, 6.5*eV, 7.0*eV, 7.5*eV, 8.0*eV, 8.27*eV, 12.40*eV, 24.80*eV, 123.99*eV}; // many divisions to make long absorption length within emission range
   
  double c_ptp = 0.3*um; // value without absorption
  //double c_ptp = 2.45*um; // value with absorption

  double m_ptp = 1e6;

  G4double ptp_l[ptp_a_entries] =
    {m_ptp*c_ptp, m_ptp*c_ptp, m_ptp*c_ptp, m_ptp*c_ptp, m_ptp*c_ptp, m_ptp*c_ptp, c_ptp, c_ptp, c_ptp, c_ptp, c_ptp, c_ptp, c_ptp, c_ptp, c_ptp};

  // Set refraction and absorption energy range and parameters

  G4double ptp_fe_energy[ptp_fe_entries] = 
    {1*eV, 3.0*eV, 3.8*eV, 6.20*eV, 8.27*eV, 12.40*eV, 24.80*eV, 123.99*eV, 200*eV};
   
  G4double ptp_n[ptp_e_entries] =
    {1.35, 1.35, 1.35, 1.35, 1.35, 1.35, 1.35, 1.35, 1.35};

  //double ptp_yield = 1.0; // yield without absorption included
  double ptp_yield = 0.9; // yield without absorption included
    // double ptp_yield = 2.75; // yield value when including absorption

  G4MaterialPropertiesTable* ptp_pt = new G4MaterialPropertiesTable();
  ptp_pt->AddProperty("RINDEX", ptp_fe_energy, ptp_n, ptp_fe_entries);
  ptp_pt->AddProperty("WLSABSLENGTH", ptp_a_energy, ptp_l, ptp_a_entries);
  ptp_pt->AddProperty("WLSCOMPONENT", ptp_e_energy, ptp_e, ptp_e_entries);
  //ptp_pt->AddConstProperty("WLSMEANNUMBERPHOTONS",ptp_yield);
  ptp_pt->AddConstProperty("WLSTIMECONSTANT", 5*ns);
  ptp_mat->SetMaterialPropertiesTable(ptp_pt); 
  
  fDefaultMaterial = env_mat;
  fmAir = mAir;
  fBase = base_mat; 
  facrylic = acrylic;
  fPTP = ptp_mat;
  fSteel = StainlessSteel;
  fAluminium = Aluminium;
  fG10 = G10;
  fMylar = Mylar;
  // DISPLAY MATERIALS
  G4cout << G4endl << *(G4Material::GetMaterialTable()) << G4endl;
  G4cout << " " << G4endl;
}

G4VPhysicalVolume* DetectorConstruction::ConstructLine()
{
  // WORLD
  //  fWorldSizeXY  = 2*m;
  //  fWorldSizeZ   = 2*m;
     
  //*************
  // WORLD VOLUME
  //*************
  
  fSolidWorld = new G4Box("World",			         //its name
		  fWorldSizeX/2*m,fWorldSizeY/2*m,fWorldSizeZ/2*m);  //its size
  
  
  fLogicWorld = new G4LogicalVolume(fSolidWorld,	//its solid
				    fDefaultMaterial,	//its material
				    "World");		//its name
  
  fPhysiWorld = new G4PVPlacement(0,			//no rotation
  				 G4ThreeVector(),	//at (0,0,0)
                                 "World",		//its name
                                 fLogicWorld,		//its logical volume
                                 NULL,			//its mother  volume
                                 false,			//no boolean operation
                                 0);			//copy number

  // Beginning of cryostat construction: coldskin, wood, foam, nougat, etc.*cm, ...*cm, warmskin

  //coldskin
  fColdSkinThickness =  0.0018; //m
  G4double Offset(0.10); //m
  G4Box* fSolidCryostat = new G4Box("Cryostat",(fCryostat_x/2+Offset)*m, (fCryostat_y/2.0+Offset)*m,(fCryostat_z/2+Offset)*m); //make it a little bigger to avoid overlaps
  G4Box* ShellOut = new G4Box("ShellOut",(fCryostat_x/2+fColdSkinThickness+Offset)*m, (fCryostat_y/2.0+fColdSkinThickness+Offset)*m,(fCryostat_z/2+fColdSkinThickness+Offset)*m);
  G4SubtractionSolid* fShell = new G4SubtractionSolid("ColdSkin", ShellOut, fSolidCryostat);
  G4LogicalVolume* fLogicShell = new G4LogicalVolume(fShell,fDUNESteel,"ColdSkin");

  
  G4VPhysicalVolume* fPhysShell = new G4PVPlacement(0,G4ThreeVector(0,0,0),"ColdSkin",
                                 fLogicShell,     //its logical volume
                                 fPhysiWorld,    	//its mother  volume
                                 false,			//no boolean operation
  						    0, true);
  
						    
  // foam
  fShieldThickness = 0.776 ; //m 
  G4Box* sOutShield = new G4Box("InShield", (fCryostat_x/2+fShieldThickness+fColdSkinThickness+Offset)*m, (fCryostat_y/2+fShieldThickness+fColdSkinThickness+Offset)*m, (fCryostat_z/2.+fShieldThickness+fColdSkinThickness+Offset)*m);
  G4SubtractionSolid *sShield = new G4SubtractionSolid("Foam",sOutShield, ShellOut);  

  fLogicShield = new G4LogicalVolume(sShield,       //shape
                             fShieldMater,            //material
                             "Foam");               //name
                                  
           new G4PVPlacement(0,                         //no rotation
			     G4ThreeVector(0.,0.,0.),  // fWorldLength/2.-1*fDetectorLength/2.),             //at (0,0,0)
			     "Foam",                  //name
			     fLogicShield,              //logical volume
			     fPhysiWorld,                      //mother  volume
                           false,                       //no boolean operation
			     0, true);                          //copy number
  

// wood
  fWoodThickness = 0.048; //m
  G4Box* sOutWood = new G4Box("InWood", (fCryostat_x/2+fWoodThickness+fShieldThickness+fColdSkinThickness+Offset)*m, (fCryostat_y/2+fWoodThickness+fShieldThickness+fColdSkinThickness+Offset)*m, (fCryostat_z/2.+fWoodThickness+fShieldThickness+fColdSkinThickness+Offset)*m);
  G4SubtractionSolid *sWood = new G4SubtractionSolid("Wood",sOutWood, sOutShield);  
  
  fLogicWood = new G4LogicalVolume(sWood,       //shape
                             fWoodMater,            //material
                             "Wood");               //name
  
           new G4PVPlacement(0,                         //no rotation
			     G4ThreeVector(0.,0.,0.),  // fWorldLength/2.-1*fDetectorLength/2.),             //at (0,0,0)
			     "Wood",                  //name
			     fLogicWood,              //logical volume
			     fPhysiWorld,                      //mother  volume
			     false,                       //no boolean operation
			     0, true);                          //copy number
  

  // warmskin
  fWarmSkinThickness =  0.024 ; //m
  G4Box* ShellOutW = new G4Box("ShellOut",(fCryostat_x/2+fWarmSkinThickness+fWoodThickness+fShieldThickness+fColdSkinThickness+Offset)*m, (fCryostat_y/2.0+fWarmSkinThickness+fWoodThickness+fShieldThickness+fColdSkinThickness+Offset)*m,(fCryostat_z/2+fWarmSkinThickness+fWoodThickness+fShieldThickness+fColdSkinThickness+Offset)*m);
  G4SubtractionSolid* fShellW = new G4SubtractionSolid("WarmSkin", ShellOutW, sOutWood);
  G4LogicalVolume* fLogicShellW = new G4LogicalVolume(fShellW,fDUNESteel,"WarmSkin");

  
  G4VPhysicalVolume* fPhysShellW = new G4PVPlacement(0,G4ThreeVector(0,0,0),
						     "WarmSkin",
						     fLogicShellW,     //its logical volume
						     fPhysiWorld,    	//its mother  volume
						     false,			//no boolean operation
						     0, true);
  
  
  // Now, need Air volume outside cryo, inside World. Otherwise this space will be G4_lAr
  
  G4SubtractionSolid* fOuterAir = new G4SubtractionSolid("OuterAir", fSolidWorld, ShellOutW);
  G4LogicalVolume* fLogicOAir = new G4LogicalVolume(fOuterAir,fmAir,"OuterAir");
  fPhysOuterAir = new G4PVPlacement(0,G4ThreeVector(0,0,0),
				    "OuterAir",
				    fLogicOAir,     //its logical volume
				    fPhysiWorld,    	//its mother  volume
				    false,			//no boolean operation
				    0, true);
  
				    
  
  std::cout << "Checking units on warm CryoSkin. xout size [mm]: " << (fCryostat_x/2+fWarmSkinThickness+fWoodThickness+fShieldThickness+fColdSkinThickness+Offset)*m << std::endl;

  // Create and Place I-Beams and Belts and Shielding panels.
  // All of this must have mother volume fPhysOuterAir
  IBeams();
  Belts();
  ShieldingWalls();
  //  ShieldingFloor();
  
  //Bulk box for wls optical properties tests

  /*G4Box* bulk = new G4Box("bulk",1.0*um,0.5*m,0.5*m);
  G4LogicalVolume* lbulk = new G4LogicalVolume(bulk,fPTP,"bulk");
  G4VPhysicalVolume* pbulk;
  pbulk = new G4PVPlacement(0,G4ThreeVector(2.0*um,0.6*m,0.0*m),"bulk",
			    lbulk,
			    fPhysiWorld,
			    true,
			    0,
			    true);

  pbulk = new G4PVPlacement(0,G4ThreeVector(-1.0*um,-1.1*m,0.0*m),"bulk",
			    lbulk,
			    fPhysiWorld,
			    true,
			    1,
			    true);*/


  
  //############## Cathode Parameters ###############
  double heightCathode=0.04; //m
  double CathodeBorder=0.04; //m
  double widthCathode=fFC_x/4.0; //in meters
  double lengthCathode=fFC_z/20.0;
  double widthCathodeVoid=(widthCathode-5.*CathodeBorder)/4;
  double lengthCathodeVoid=(lengthCathode-5.*CathodeBorder)/4;
  
  
  G4Box* fCathodeBlock = new G4Box("CathodeBlock",(widthCathode/2.)*m,(heightCathode/2.)*m,(lengthCathode/2.)*m);
  G4Box* fCathodeVoid = new G4Box("CathodeVoid",(widthCathodeVoid/2.)*m,(heightCathode)*m,(lengthCathodeVoid/2.)*m);
  G4SubtractionSolid* Cathode1 = new G4SubtractionSolid("Cathode1", fCathodeBlock, fCathodeVoid, 0, G4ThreeVector((-1.5*widthCathodeVoid-2.*CathodeBorder)*m,0,(-1.5*lengthCathodeVoid-2.*CathodeBorder)*m));
  G4SubtractionSolid* Cathode2 = new G4SubtractionSolid("Cathode2", Cathode1, fCathodeVoid, 0, G4ThreeVector((-1.5*widthCathodeVoid-2.*CathodeBorder)*m,0,(-0.5*lengthCathodeVoid-1.*CathodeBorder)*m));
  G4SubtractionSolid* Cathode3 = new G4SubtractionSolid("Cathode3", Cathode2, fCathodeVoid, 0, G4ThreeVector((-1.5*widthCathodeVoid-2.*CathodeBorder)*m,0,(0.5*lengthCathodeVoid+1.*CathodeBorder)*m));
  G4SubtractionSolid* Cathode4 = new G4SubtractionSolid("Cathode4", Cathode3, fCathodeVoid, 0, G4ThreeVector((-1.5*widthCathodeVoid-2.*CathodeBorder)*m,0,(1.5*lengthCathodeVoid+2.*CathodeBorder)*m));
  G4SubtractionSolid* Cathode5 = new G4SubtractionSolid("Cathode5", Cathode4, fCathodeVoid, 0, G4ThreeVector((-0.5*widthCathodeVoid-1.*CathodeBorder)*m,0,(-1.5*lengthCathodeVoid-2.*CathodeBorder)*m));
  G4SubtractionSolid* Cathode6 = new G4SubtractionSolid("Cathode6", Cathode5, fCathodeVoid, 0, G4ThreeVector((-0.5*widthCathodeVoid-1.*CathodeBorder)*m,0,(-0.5*lengthCathodeVoid-1.*CathodeBorder)*m));
  G4SubtractionSolid* Cathode7 = new G4SubtractionSolid("Cathode7", Cathode6, fCathodeVoid, 0, G4ThreeVector((-0.5*widthCathodeVoid-1.*CathodeBorder)*m,0,(0.5*lengthCathodeVoid+1.*CathodeBorder)*m));
  G4SubtractionSolid* Cathode8 = new G4SubtractionSolid("Cathode8", Cathode7, fCathodeVoid, 0, G4ThreeVector((-0.5*widthCathodeVoid-1.*CathodeBorder)*m,0,(1.5*lengthCathodeVoid+2.*CathodeBorder)*m));
  G4SubtractionSolid* Cathode9 = new G4SubtractionSolid("Cathode9", Cathode8, fCathodeVoid, 0, G4ThreeVector((0.5*widthCathodeVoid+1.*CathodeBorder)*m,0,(-1.5*lengthCathodeVoid-2.*CathodeBorder)*m));
  G4SubtractionSolid* Cathode10 = new G4SubtractionSolid("Cathode10", Cathode9, fCathodeVoid, 0, G4ThreeVector((0.5*widthCathodeVoid+1.*CathodeBorder)*m,0,(-0.5*lengthCathodeVoid-1.*CathodeBorder)*m));
  G4SubtractionSolid* Cathode11 = new G4SubtractionSolid("Cathode11", Cathode10, fCathodeVoid, 0, G4ThreeVector((0.5*widthCathodeVoid+1.*CathodeBorder)*m,0,(0.5*lengthCathodeVoid+1.*CathodeBorder)*m));
  G4SubtractionSolid* Cathode12 = new G4SubtractionSolid("Cathode12", Cathode11, fCathodeVoid, 0, G4ThreeVector((0.5*widthCathodeVoid+1.*CathodeBorder)*m,0,(1.5*lengthCathodeVoid+2.*CathodeBorder)*m));
G4SubtractionSolid* Cathode13 = new G4SubtractionSolid("Cathode13", Cathode12, fCathodeVoid, 0, G4ThreeVector((1.5*widthCathodeVoid+2.*CathodeBorder)*m,0,(-1.5*lengthCathodeVoid-2.*CathodeBorder)*m));
  G4SubtractionSolid* Cathode14 = new G4SubtractionSolid("Cathode14", Cathode13, fCathodeVoid, 0, G4ThreeVector((1.5*widthCathodeVoid+2.*CathodeBorder)*m,0,(-0.5*lengthCathodeVoid-1.*CathodeBorder)*m));
  G4SubtractionSolid* Cathode15 = new G4SubtractionSolid("Cathode15", Cathode14, fCathodeVoid, 0, G4ThreeVector((1.5*widthCathodeVoid+2.*CathodeBorder)*m,0,(0.5*lengthCathodeVoid+1.*CathodeBorder)*m));
  G4SubtractionSolid* CathodeGrid = new G4SubtractionSolid("CathodeGrid", Cathode15, fCathodeVoid, 0, G4ThreeVector((1.5*widthCathodeVoid+2.*CathodeBorder)*m,0,(1.5*lengthCathodeVoid+2.*CathodeBorder)*m));
  G4LogicalVolume* fLogicCathode = new G4LogicalVolume(CathodeGrid,fSteel,"Cathode");
  

  double cxpos=-1.5*widthCathode;
  double czpos=-9.5*lengthCathode;

  int cnt=0;

  for(int i = 0; i<4; i++){ // copies of Cathode
    cxpos=(i-1.5)*widthCathode;
    for(int k = 0; k<20; k++){
      czpos=(k-9.5)*lengthCathode;
      G4PVPlacement* cathode_cp  = new G4PVPlacement(0,G4ThreeVector(cxpos*m,0,czpos*m), "Cathode", fLogicCathode, fPhysiWorld, false, cnt, true);
      cnt++;
    }
  }
  
  G4Box* fSolidAnode = new G4Box("Anode",(fCathode_x/2)*m,fthickness/2*m,(fCathode_z/2)*m);

  G4LogicalVolume* fLogicAnodeT = new G4LogicalVolume(fSolidAnode,fSteel,"AnodeT");
  G4VPhysicalVolume* fPhysiAnodeT = new G4PVPlacement(0,G4ThreeVector(0,(fCryostat_y/2.0+fthickness/2)*m,0),"AnodeT",
                                 fLogicAnodeT,     //its logical volume
                                 fPhysiWorld,    	//its mother  volume
                                 false,			//no boolean operation
				 0,
				 true); //check for overlaps

  G4LogicalVolume* fLogicAnodeB = new G4LogicalVolume(fSolidAnode,fSteel,"AnodeB");
  G4VPhysicalVolume* fPhysiAnodeB = new G4PVPlacement(0,G4ThreeVector(0,-(fCryostat_y/2.0+fthickness/2)*m,0),"AnodeB",
                                 fLogicAnodeB,     //its logical volume
                                 fPhysiWorld,    	//its mother  volume
                                 false,			//no boolean operation
				 0,
				 true); //check for overlaps   

  //FC Structure

  //26 c-shaped profiles per FC module/each modelled as eliptical cutted tubes
  //Longer lateral

  G4double ypos=-fCryostat_y/2.0 + fFCOut_y/2.0; //in m
  G4double zpos=-fCryostat_z/2.0 + fFCOut_z/2.0; //in m
  
  G4EllipticalTube* fFieldCageOut = new G4EllipticalTube("fFieldCageOut",fFCOut_x*m,fFCOut_y*m,(fFCOut_z-0.001)*m);
  G4EllipticalTube* fFieldCageIn = new G4EllipticalTube("fFieldCageIn",(fFCOut_x-0.002)*m,(fFCOut_y-0.002)*m,(fFCOut_z-0.001)*m);
  G4SubtractionSolid* fFCAux = new G4SubtractionSolid("FieldCageAux", fFieldCageOut, fFieldCageIn, 0, G4ThreeVector(0.,0.,0.));
  G4Box* fFCcutout = new G4Box("fFCcutout",0.02*m,fFCOut_y*m,(fFCOut_z+0.001)*m);
  G4SubtractionSolid* fSolidFC = new G4SubtractionSolid("FieldCage", fFCAux,  fFCcutout, 0, G4ThreeVector((-0.02-0.015)*m,0.,0.));
  G4LogicalVolume* fLogicFC = new G4LogicalVolume(fSolidFC,fAluminium,"FieldCage");

  cnt=0;
  int ncols=20, nrows=26;

  G4RotationMatrix* fcMirror = new G4RotationMatrix();
  G4ThreeVector* axisfcMirror = new G4ThreeVector(0.0,0.0,1.0);
  fcMirror->rotate(CLHEP::pi,axisfcMirror);

  
  for(int i=0; i<ncols; i++){    
    for(int j=0; j<nrows; j++){
      G4cout << "i: " << i << " j: " << j << G4endl;
      ypos=-fCryostat_y/2.0+(j+0.5)*(fCryostat_y/2.0/nrows);
      zpos=-fCryostat_z/2.0+(i+0.5)*2*fFCOut_z;
      G4PVPlacement* ph_cp = new G4PVPlacement(0,G4ThreeVector((fFC_x/2.-fAra_offset+fFCOut_x-0.015)*m,ypos*m,zpos*m),
						"FieldCage", fLogicFC, fPhysiWorld, false, cnt, true);
      cnt++;
      ypos=fCryostat_y/2.0-(j+0.5)*(fCryostat_y/2.0/nrows);
      zpos=-fCryostat_z/2.0+(i+0.5)*2*fFCOut_z;
      G4PVPlacement* ph_cp2 = new G4PVPlacement(0,G4ThreeVector((fFC_x/2.-fAra_offset+fFCOut_x-0.015)*m,ypos*m,zpos*m),
						"FieldCage", fLogicFC, fPhysiWorld, false, cnt, true);
      G4cout << "ypos: " << ypos << " zpos: " << zpos << G4endl;
      cnt++;

      ypos=-fCryostat_y/2.0+(j+0.5)*(fCryostat_y/2.0/nrows);
      zpos=-fCryostat_z/2.0+(i+0.5)*2*fFCOut_z;
      G4PVPlacement* ph_cp3 = new G4PVPlacement(fcMirror,G4ThreeVector(-(fFC_x/2.-fAra_offset+fFCOut_x-0.015)*m,ypos*m,zpos*m),
						"FieldCage", fLogicFC, fPhysiWorld, false, cnt, true);
      cnt++;
      ypos=fCryostat_y/2.0-(j+0.5)*(fCryostat_y/2.0/nrows);
      zpos=-fCryostat_z/2.0+(i+0.5)*2*fFCOut_z;
      G4PVPlacement* ph_cp4 = new G4PVPlacement(fcMirror,G4ThreeVector(-(fFC_x/2.-fAra_offset+fFCOut_x-0.015)*m,ypos*m,zpos*m),
						"FieldCage", fLogicFC, fPhysiWorld, false, cnt, true);
      G4cout << "ypos: " << ypos << " zpos: " << zpos << G4endl;
      cnt++;   
    }
  }


  ypos = 0.0;
  zpos = 0.0;
  G4double xpos = 0.0;
  ncols=4;

  //shorter laterals
  G4RotationMatrix* rSh = new G4RotationMatrix();
  G4ThreeVector* axisSh = new G4ThreeVector(0.0,1.0,0.0);
  rSh->rotate(CLHEP::pi/2,axisSh);

  G4RotationMatrix* lSh = new G4RotationMatrix();
  lSh->rotate(-CLHEP::pi/2,axisSh);

  for(int i=0; i<ncols; i++){    
    for(int j=0; j<nrows; j++){
      ypos=-fCryostat_y/2.0+(j+0.5)*(fCryostat_y/2.0/nrows);
      xpos=-fFC_x/2.0+fFCOut_z/2+(i+0.5)*2*fFCOut_z;
      zpos=fCryostat_z/2.0-fAra_offset+fFCOut_x-0.015;
      G4PVPlacement* ph_cp5 = new G4PVPlacement(rSh,G4ThreeVector(xpos*m,ypos*m,zpos*m),
						"FieldCageS", fLogicFC, fPhysiWorld, false, cnt, true);
      cnt++;
      ypos=-(-fCryostat_y/2.0+(j+0.5)*(fCryostat_y/2.0/nrows));
      xpos=-fFC_x/2.0+fFCOut_z/2+(i+0.5)*2*fFCOut_z;
      zpos=fCryostat_z/2.0-fAra_offset+fFCOut_x-0.015;
      G4PVPlacement* ph_cp6 = new G4PVPlacement(rSh,G4ThreeVector(xpos*m,ypos*m,zpos*m),

						"FieldCageS", fLogicFC, fPhysiWorld, false, cnt, true);
      cnt++;      
      ypos=-fCryostat_y/2.0+(j+0.5)*(fCryostat_y/2.0/nrows);
      xpos=-fFC_x/2.0+fFCOut_z/2+(i+0.5)*2*fFCOut_z;
      zpos=-(fCryostat_z/2.0-fAra_offset+fFCOut_x-0.015);
      G4PVPlacement* ph_cp7 = new G4PVPlacement(lSh,G4ThreeVector(xpos*m,ypos*m,zpos*m),
						"FieldCageS", fLogicFC, fPhysiWorld, false, cnt, true);
      cnt++;      
      ypos=-(-fCryostat_y/2.0+(j+0.5)*(fCryostat_y/2.0/nrows));
      xpos=-fFC_x/2.0+fFCOut_z/2+(i+0.5)*2*fFCOut_z;
      zpos=-(fCryostat_z/2.0-fAra_offset+fFCOut_x-0.015);
      G4PVPlacement* ph_cp8 = new G4PVPlacement(lSh,G4ThreeVector(xpos*m,ypos*m,zpos*m),
						"FieldCageS", fLogicFC, fPhysiWorld, false, cnt, true);
      cnt++;
      
     }
  }
   
   //ARAPUCAs
   ncols=120, nrows=12;
   zpos=0.0;
   G4double darapuca = fAra_z;
   G4Box* Arapuca = new G4Box("Arapuca",fAra_x/2*m,fAra_y/2*m,(darapuca/2-1e-3)*m);
   G4LogicalVolume* fLogicArapuca = new G4LogicalVolume(Arapuca,facrylic,"Arapuca");

   G4Box* ptp_film = new G4Box("PTP_film",fptp_width/2.0,fAra_y/2*m,(darapuca/2-1e-3)*m);
   G4LogicalVolume* fLogicPTP = new G4LogicalVolume(ptp_film,fPTP,"PTP_film");
   G4LogicalVolume* fLogicAB = new G4LogicalVolume(ptp_film,fMylar,"ArapucaBackCoating");
   
  int ptp_cnt = 0; G4VPhysicalVolume* ptp_phys;
  int myl_cnt = 0; G4VPhysicalVolume* myl_phys;
  G4String name, physname, name2, physname2;
  for(int i=0; i<ncols; i++){
    ypos=(-fCryostat_y/2.0) + (fCryostat_y/2.0 - (nrows*fAra_y+(nrows-1)*fAras_yspacing))/2 + fAra_y/2;
    //lower edge of cryostat + space between closest edges of cryostat & lower XARAPUCA + y_arapuca_side/2 in m
    
    zpos=(-fCryostat_z/2.0) + (i+0.5)*darapuca; //in m;
    for(int j=0; j<nrows;j++){
      std::cout << "ypos: " << ypos << " zpos: " << zpos << std::endl;
      
      ptp_phys = new G4PVPlacement(0,G4ThreeVector((-fFC_x/2.+fAra_offset)*m+fAra_x*m+fptp_width/2,ypos*m,zpos*m),"PTP_film", fLogicPTP, fPhysiWorld, false,ptp_cnt, true);
      ptp_cnt++;
      ptp_phys = new G4PVPlacement(0,G4ThreeVector((fFC_x/2.-fAra_offset)*m-fAra_x*m-fptp_width/2,ypos*m,zpos*m),"PTP_film", fLogicPTP, fPhysiWorld, false,ptp_cnt, true);
      ptp_cnt++;
      
      name = "ArapucaR"; name.append(std::to_string(i+1)); name.append("_"); name.append(std::to_string(j+1));
      physname = "fPhysArapucaR"; physname.append(std::to_string(i+1)); physname.append("_"); physname.append(std::to_string(j+1));
      new G4PVPlacement(0,G4ThreeVector((fFC_x/2.-fAra_offset-fAra_x/2)*m,ypos*m,zpos*m),name.c_str(), fLogicArapuca, fPhysiWorld, false,0, true);
      name2 = "ArapucaL"; name2.append(std::to_string(i+1)); name2.append("_"); name2.append(std::to_string(j+1));
      physname2 = "fPhysArapucaL"; physname2.append(std::to_string(i+1)); physname2.append("_"); physname2.append(std::to_string(j+1));
      new G4PVPlacement(0,G4ThreeVector((-fFC_x/2.+fAra_offset+fAra_x/2)*m,ypos*m,zpos*m),name2.c_str(), fLogicArapuca, fPhysiWorld, false,0, true);

      myl_phys = new G4PVPlacement(0,G4ThreeVector((-fFC_x/2.+fAra_offset)*m-fptp_width/2,ypos*m,zpos*m),"MYL_film", fLogicAB, fPhysiWorld, false,myl_cnt, true);
      myl_cnt++;
      myl_phys = new G4PVPlacement(0,G4ThreeVector((fFC_x/2.-fAra_offset)*m+fptp_width/2,ypos*m,zpos*m),"MYL_film", fLogicAB, fPhysiWorld, false,myl_cnt, true);
      myl_cnt++;
      
      ypos+=(fAra_y+fAras_yspacing);
    }
    
    ypos=(fCryostat_y/2.0) - (fCryostat_y/2.0 - (nrows*fAra_y+(nrows-1)*fAras_yspacing))/2 - fAra_y/2; //in m
    
    for(int j=nrows; j<2*nrows;j++){
      std::cout << "ypos: " << ypos << " zpos: " << zpos << std::endl;
      
      ptp_phys = new G4PVPlacement(0,G4ThreeVector((-fFC_x/2.+fAra_offset)*m+fAra_x*m+fptp_width/2,ypos*m,zpos*m),"PTP_film", fLogicPTP, fPhysiWorld, false,ptp_cnt, true);
      ptp_cnt++;
      ptp_phys = new G4PVPlacement(0,G4ThreeVector((fFC_x/2.-fAra_offset)*m-fAra_x*m-fptp_width/2,ypos*m,zpos*m),"PTP_film", fLogicPTP, fPhysiWorld, false,ptp_cnt, true);
      ptp_cnt++;
      
      name = "ArapucaR"; name.append(std::to_string(i+1)); name.append("_"); name.append(std::to_string(j+1));
      physname = "fPhysArapucaR"; physname.append(std::to_string(i+1)); physname.append("_"); physname.append(std::to_string(j+1));
      new G4PVPlacement(0,G4ThreeVector((fFC_x/2.-fAra_offset-fAra_x/2)*m,ypos*m,zpos*m),name.c_str(), fLogicArapuca, fPhysiWorld, false,0, true);
      name2 = "ArapucaL"; name2.append(std::to_string(i+1)); name2.append("_"); name2.append(std::to_string(j+1));
      physname2 = "fPhysArapucaL"; physname2.append(std::to_string(i+1)); physname2.append("_"); physname2.append(std::to_string(j+1));
      new G4PVPlacement(0,G4ThreeVector((-fFC_x/2.+fAra_offset+fAra_x/2)*m,ypos*m,zpos*m),name2.c_str(), fLogicArapuca, fPhysiWorld, false,0, true);

      myl_phys = new G4PVPlacement(0,G4ThreeVector((-fFC_x/2.+fAra_offset)*m-fptp_width/2,ypos*m,zpos*m),"MYL_film", fLogicAB, fPhysiWorld, false,myl_cnt, true);
      myl_cnt++;
      myl_phys = new G4PVPlacement(0,G4ThreeVector((fFC_x/2.-fAra_offset)*m+fptp_width/2,ypos*m,zpos*m),"MYL_film", fLogicAB, fPhysiWorld, false,myl_cnt, true);
      myl_cnt++;
      
      ypos-=(fAra_y+fAras_yspacing);
    }
    
  }
  
  //Shorter side
  G4Box* ArapucaShort = new G4Box("ArapucaShort",fAra_x/2*m,fAra_y/2*m,(darapuca/2-1e-3)*m);
  G4LogicalVolume* fLogicArapucaShort = new G4LogicalVolume(ArapucaShort,facrylic,"ArapucaShort");
  
  G4Box* ptps_film = new G4Box("PTPs_film",fptp_width/2.0,fAra_y/2*m,(darapuca/2-1e-3)*m);
  G4LogicalVolume* fLogicPTPs = new G4LogicalVolume(ptps_film,fPTP,"PTPs_film");
  G4LogicalVolume* fLogicABs = new G4LogicalVolume(ptps_film,fMylar,"ArapucaBackCoatingShort");
   
  ncols=24, nrows=12;
  int ptps_cnt = 0; G4VPhysicalVolume* ptps_phys;
  int myls_cnt = 0; G4VPhysicalVolume* myls_phys;
  
  //G4double xpos;
  xpos=0;
  
  for(int i=0; i<ncols; i++){
    ypos=-fCryostat_y/2.0 + (fCryostat_y/2.0 - (nrows*fAra_y+(nrows-1)*fAras_yspacing))/2 + fAra_y/2; //in m
    xpos = -fFC_x/2+(fFC_x-(ncols*darapuca+(ncols-1)*fAras_yspacing))/2.0 + darapuca/2 + i*(darapuca+fAras_yspacing);
    std::cout << xpos << std::endl;
    for(int j=0; j<nrows;j++){
      
      ptps_phys = new G4PVPlacement(rSh,G4ThreeVector(xpos*m,ypos*m,(fFC_z/2.-fAra_offset)*m-fAra_x*m-fptp_width/2.0),"PTPs_film", fLogicPTPs, fPhysiWorld, false,ptps_cnt, true);
      ptps_cnt++;
      ptps_phys = new G4PVPlacement(rSh,G4ThreeVector(xpos*m,ypos*m,(-fFC_z/2.+fAra_offset)*m+fAra_x*m+fptp_width/2.0),"PTPs_film", fLogicPTPs, fPhysiWorld, false,ptps_cnt, true);
      ptps_cnt++;
      
      name = "ArapucaF"; name.append(std::to_string(i+1)); name.append("_"); name.append(std::to_string(j+1));
      physname = "fPhysArapucaF"; physname.append(std::to_string(i+1)); physname.append("_"); physname.append(std::to_string(j+1));
      new G4PVPlacement(rSh,G4ThreeVector(xpos*m,ypos*m,(fFC_z/2.-fAra_offset-fAra_x/2)*m),name.c_str(), fLogicArapucaShort, fPhysiWorld, false,0, true);
      name2 = "ArapucaB"; name2.append(std::to_string(i+1)); name2.append("_"); name2.append(std::to_string(j+1));
      physname2 = "fPhysArapucaB"; physname2.append(std::to_string(i+1)); physname2.append("_"); physname2.append(std::to_string(j+1));
      new G4PVPlacement(rSh,G4ThreeVector(xpos*m,ypos*m,(-fFC_z/2.+fAra_offset+fAra_x/2)*m),name2.c_str(), fLogicArapucaShort, fPhysiWorld, false,0, true);

      myls_phys = new G4PVPlacement(rSh,G4ThreeVector(xpos*m,ypos*m,(fFC_z/2.-fAra_offset+fptp_width/2.0)*m),"MYL_film", fLogicABs, fPhysiWorld, false,myls_cnt, true);
      myls_cnt++;
      myls_phys = new G4PVPlacement(rSh,G4ThreeVector(xpos*m,ypos*m,(-fFC_z/2.+fAra_offset-fptp_width/2.0)*m),"MYL_film", fLogicABs, fPhysiWorld, false,myls_cnt, true);
      myls_cnt++;
      
      ypos+=(fAra_y+fAras_yspacing);
    }
    
    ypos=fCryostat_y/2.0 - (fCryostat_y/2.0 - (nrows*fAra_y+(nrows-1)*fAras_yspacing))/2 - fAra_y/2; //in m
    for(int j=nrows; j<2*nrows;j++){
      
      ptps_phys = new G4PVPlacement(rSh,G4ThreeVector(xpos*m,ypos*m,(fFC_z/2.-fAra_offset)*m-fAra_x*m-fptp_width/2.0),"PTPs_film", fLogicPTPs, fPhysiWorld, false,ptps_cnt, true);
      ptps_cnt++;
      ptps_phys = new G4PVPlacement(rSh,G4ThreeVector(xpos*m,ypos*m,(-fFC_z/2.+fAra_offset)*m+fAra_x*m+fptp_width/2.0),"PTPs_film", fLogicPTPs, fPhysiWorld, false,ptps_cnt, true);
      ptps_cnt++;
      
      name = "ArapucaF"; name.append(std::to_string(i+1)); name.append("_"); name.append(std::to_string(j+1));
      physname = "fPhysArapucaF"; physname.append(std::to_string(i+1)); physname.append("_"); physname.append(std::to_string(j+1));
      new G4PVPlacement(rSh,G4ThreeVector(xpos*m,ypos*m,(fFC_z/2.-fAra_offset-fAra_x/2)*m),name.c_str(), fLogicArapucaShort, fPhysiWorld, false,0, true);
      name2 = "ArapucaB"; name2.append(std::to_string(i+1)); name2.append("_"); name2.append(std::to_string(j+1));
      physname2 = "fPhysArapucaB"; physname2.append(std::to_string(i+1)); physname2.append("_"); physname2.append(std::to_string(j+1));
      new G4PVPlacement(rSh,G4ThreeVector(xpos*m,ypos*m,(-fFC_z/2.+fAra_offset+fAra_x/2)*m),name2.c_str(), fLogicArapucaShort, fPhysiWorld, false,0, true);

      myls_phys = new G4PVPlacement(rSh,G4ThreeVector(xpos*m,ypos*m,(fFC_z/2.-fAra_offset+fptp_width/2.0)*m),"MYL_film", fLogicABs, fPhysiWorld, false,myls_cnt, true);
      myls_cnt++;
      myls_phys = new G4PVPlacement(rSh,G4ThreeVector(xpos*m,ypos*m,(-fFC_z/2.+fAra_offset-fptp_width/2.0)*m),"MYL_film", fLogicABs, fPhysiWorld, false,myls_cnt, true);
      myls_cnt++;
      
      ypos-=(fAra_y+fAras_yspacing);
    }
    
  }
  
  //vertical bars on longer laterals
  ncols = 60;
  G4Box* VerticalBar = new G4Box("VerticalBar",fvert_bar_x/2*m,fvert_bar_y/2*m,fvert_bar_z/2*m);
  G4LogicalVolume* fLogicVerticalBar = new G4LogicalVolume(VerticalBar,facrylic,"VerticalBar");

  int vbar_cnt = 0; G4VPhysicalVolume* vbar_phys;
  
  for(int i=0; i<ncols; i++){
    zpos = 2*fAra_z*(i+0.5) -fFC_z/2;
    //std::cout << (fFC_x/2.-fAra_offset-fAra_x-fvert_bar_x/2) <<" "<< zpos << std::endl;    
    //vbar_phys = new G4PVPlacement(0,G4ThreeVector((fFC_x/2.-fAra_offset-fAra_x-fvert_bar_x/2)*m-fptp_width,0*m,zpos*m),"vbar", fLogicVerticalBar, fPhysiWorld, false, vbar_cnt, true);
    vbar_phys = new G4PVPlacement(0,G4ThreeVector((fFC_x/2.-fAra_offset-fAra_x-fvert_bar_x/2)*m-fptp_width,fFC_y/4*m,zpos*m),"vbar", fLogicVerticalBar, fPhysiWorld, false, vbar_cnt, true);
    vbar_cnt++;
    vbar_phys = new G4PVPlacement(0,G4ThreeVector((-fFC_x/2.+fAra_offset+fAra_x+fvert_bar_x/2)*m+fptp_width,fFC_y/4*m,zpos*m),"vbar", fLogicVerticalBar, fPhysiWorld, false, vbar_cnt, true);
    vbar_cnt++;
    vbar_phys = new G4PVPlacement(0,G4ThreeVector((fFC_x/2.-fAra_offset-fAra_x-fvert_bar_x/2)*m-fptp_width,-fFC_y/4*m,zpos*m),"vbar", fLogicVerticalBar, fPhysiWorld, false, vbar_cnt, true);
    vbar_cnt++;
    vbar_phys = new G4PVPlacement(0,G4ThreeVector((-fFC_x/2.+fAra_offset+fAra_x+fvert_bar_x/2)*m+fptp_width,-fFC_y/4*m,zpos*m),"vbar", fLogicVerticalBar, fPhysiWorld, false, vbar_cnt, true);
    vbar_cnt++;
  }

  //vertical bars on shorter laterals
  ncols = 12;
  int ncolx = 24;

  G4Box* VerticalBarS = new G4Box("VerticalBarS",fvert_bar_x/2*m,fvert_bar_y/2*m,fvert_bar_z/2*m);
  G4LogicalVolume* fLogicVerticalBarS = new G4LogicalVolume(VerticalBarS,facrylic,"VerticalBarS");

  vbar_cnt = 0;
  
  for(int i=0; i<ncols; i++){

    //zpos;
    xpos = -fFC_x/2+(fFC_x-(ncolx*darapuca+(ncolx-1)*fAras_yspacing))/2.0 + darapuca + fAras_yspacing/2 + i*2*(darapuca + fAras_yspacing);
    
    //std::cout << xpos  <<" "<< zpos << std::endl;    
    
    vbar_phys = new G4PVPlacement(rSh,G4ThreeVector(xpos*m,fFC_y/4*m,(fFC_z/2.-fAra_offset-fAra_x-fvert_bar_x/2)*m-fptp_width),"vbars", fLogicVerticalBarS, fPhysiWorld, false, vbar_cnt, true);
    vbar_cnt++;
    vbar_phys = new G4PVPlacement(rSh,G4ThreeVector(xpos*m,fFC_y/4*m,(-fFC_z/2.+fAra_offset+fAra_x+fvert_bar_x/2)*m+fptp_width),"vbars", fLogicVerticalBarS, fPhysiWorld, false, vbar_cnt, true);
    vbar_cnt++;
    vbar_phys = new G4PVPlacement(rSh,G4ThreeVector(xpos*m,-fFC_y/4*m,(fFC_z/2.-fAra_offset-fAra_x-fvert_bar_x/2)*m-fptp_width),"vbars", fLogicVerticalBarS, fPhysiWorld, false, vbar_cnt, true);
    vbar_cnt++;
    vbar_phys = new G4PVPlacement(rSh,G4ThreeVector(xpos*m,-fFC_y/4*m,(-fFC_z/2.+fAra_offset+fAra_x+fvert_bar_x/2)*m+fptp_width),"vbars", fLogicVerticalBarS, fPhysiWorld, false, vbar_cnt, true);
    vbar_cnt++;
    
  }
  
    
  
  //Surfaces setup

  //lAr-Anode interface

  G4OpticalSurface* AnodeSurfaceT = new G4OpticalSurface("AnodeSurfaceT");
  AnodeSurfaceT->SetType(dielectric_metal);
  AnodeSurfaceT->SetModel(unified);
  AnodeSurfaceT->SetFinish(ground);
  AnodeSurfaceT->SetSigmaAlpha(0.0*deg); // for vikuit

  G4OpticalSurface* AnodeSurfaceB = new G4OpticalSurface("AnodeSurfaceB");
  AnodeSurfaceB->SetType(dielectric_metal);
  AnodeSurfaceB->SetModel(unified);
  AnodeSurfaceB->SetFinish(ground);
  AnodeSurfaceB->SetSigmaAlpha(0.0*deg); // for vikuit
  
  const G4int nEntries = 8;
  G4double PhotonEnergy[nEntries] =
    { 2.5*eV, 5.0*eV, 7.0*eV, 7.5*eV, 8.0*eV, 9.0*eV, 9.5*eV, 10.136*eV};

  G4double Anode_r[nEntries] = {0.20, 0.20, 0.20, 0.20, 0.0, 0.0, 0.0, 0.0}; //reflection coef for base
  //G4double Anode_r[nEntries] = {0.5, 0.5, 0.5, 0.5, 0.0, 0.0, 0.0, 0.0}; //reflection coef for base
  G4double Anode_e[nEntries] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0}; //absorption coefficient

  G4MaterialPropertiesTable* AnodeSurface_pt = new G4MaterialPropertiesTable();

  AnodeSurface_pt->AddProperty("REFLECTIVITY", PhotonEnergy, Anode_r, nEntries);
  AnodeSurface_pt->AddProperty("EFFICIENCY", PhotonEnergy, Anode_e, nEntries);

  AnodeSurfaceT->SetMaterialPropertiesTable(AnodeSurface_pt);
  new G4LogicalSkinSurface("AnodeSurface", fLogicAnodeT, AnodeSurfaceT);

  AnodeSurfaceB->SetMaterialPropertiesTable(AnodeSurface_pt);
  new G4LogicalSkinSurface("AnodeSurface", fLogicAnodeB, AnodeSurfaceB);

  //_____________________________FC REFLECTIVITY CHANGE 0.2 -> 0.7 ______________________________
  G4OpticalSurface* FCSurface = new G4OpticalSurface("FCSurface");
  FCSurface->SetType(dielectric_metal);
  FCSurface->SetModel(unified);
  FCSurface->SetFinish(ground);
  FCSurface->SetSigmaAlpha(0.0*deg); // for vikuit

  G4double FC_r[nEntries] = {0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7}; //reflection coef for base
  G4double FC_e[nEntries] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0}; //absorption coefficient
  //G4double FC_r[nEntries] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}; //reflection coef for base
  
  G4MaterialPropertiesTable* FCSurface_pt = new G4MaterialPropertiesTable();
  
  FCSurface_pt->AddProperty("REFLECTIVITY", PhotonEnergy, FC_r, nEntries);
  FCSurface_pt->AddProperty("EFFICIENCY", PhotonEnergy, FC_e, nEntries);
  
  FCSurface->SetMaterialPropertiesTable(FCSurface_pt);
  new G4LogicalSkinSurface("FCSurface",fLogicFC,FCSurface);
  //  new G4LogicalSkinSurface("FCSurfaceShort",fLogicFCShort,FCSurface);

  //_____________________________Mylar REFLECTIVITY 1.0 ______________________________
  G4OpticalSurface* ArapucaBackSurface = new G4OpticalSurface("ArapucaBackSurface");
  ArapucaBackSurface->SetType(dielectric_metal);
  ArapucaBackSurface->SetModel(unified);
  ArapucaBackSurface->SetFinish(ground);
  ArapucaBackSurface->SetSigmaAlpha(0.0*deg); // for vikuit

  G4double AB_r[nEntries] = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0}; //reflection coef for base
  G4double AB_e[nEntries] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0}; //absorption coefficient
  
  G4MaterialPropertiesTable* ArapucaBackSurface_pt = new G4MaterialPropertiesTable();
  
  ArapucaBackSurface_pt->AddProperty("REFLECTIVITY", PhotonEnergy, AB_r, nEntries);
  ArapucaBackSurface_pt->AddProperty("EFFICIENCY", PhotonEnergy, AB_e, nEntries);
  
  ArapucaBackSurface->SetMaterialPropertiesTable(ArapucaBackSurface_pt);
  new G4LogicalSkinSurface("ArapucaBackSurface",fLogicAB,ArapucaBackSurface);
  new G4LogicalSkinSurface("ArapucaBackSurfaceShort",fLogicABs,ArapucaBackSurface);


  // _________________ CHANGE NEW CRYOSTAT PROPERTIES _____________________________________
  G4OpticalSurface* CryostatSurface = new G4OpticalSurface("CryostatSurface");
  CryostatSurface->SetType(dielectric_metal);
  CryostatSurface->SetModel(unified);
  CryostatSurface->SetFinish(ground);
  CryostatSurface->SetSigmaAlpha(0.0*deg); // for vikuit

  G4double Cryo_r[nEntries] = {0.4, 0.4, 0.4, 0.4, 0.3, 0.3, 0.3, 0.3}; //reflection coef for base
  G4double Cryo_e[nEntries] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0}; //absorption coefficient
  
  G4MaterialPropertiesTable* CryostatSurface_pt = new G4MaterialPropertiesTable();

  CryostatSurface_pt->AddProperty("REFLECTIVITY", PhotonEnergy, Cryo_r, nEntries);
  CryostatSurface_pt->AddProperty("EFFICIENCY", PhotonEnergy, Cryo_e, nEntries);
  
  CryostatSurface->SetMaterialPropertiesTable(CryostatSurface_pt);
  new G4LogicalSkinSurface("CryostatSurface", fLogicShell, CryostatSurface);
  
  G4VisAttributes* simpleWorldVisAtt= new G4VisAttributes(G4Colour(1.0,1.0,1.0)); //White
  simpleWorldVisAtt->SetVisibility(true);
  
  //G4VisAttributes* simplePlain= new G4VisAttributes(G4Colour(1.0,1.0,1.0,0.4));   
  G4VisAttributes* simplePlain= new G4VisAttributes(G4Colour(1.0,1.0,1.0,0.5)); //White
  simplePlain->SetVisibility(true);
  simplePlain->SetForceSolid(true);
  simplePlain->SetForceAuxEdgeVisible(true);
   
  //G4VisAttributes* simpleBoxAtt= new G4VisAttributes(G4Colour(1.0,1.0,0.0,0.3));
  G4VisAttributes* simpleBoxAtt= new G4VisAttributes(G4Colour(1.0,1.0,0.0,0.5));
  simpleBoxAtt->SetDaughtersInvisible(true);
  simpleBoxAtt->SetForceSolid(true);
  simpleBoxAtt->SetForceAuxEdgeVisible(true);
  
  G4VisAttributes* simpleBoxAttKGM= new G4VisAttributes(G4Colour(0.0,0.0,1.0,1.0));
  simpleBoxAttKGM->SetVisibility(true);
  simpleBoxAttKGM->SetForceWireframe(true);
  simpleBoxAttKGM->SetDaughtersInvisible(false);
  simpleBoxAttKGM->SetForceSolid(true);

  G4VisAttributes* BoxAtt= new G4VisAttributes(G4Colour(1.0,0.0,0.0,1.0));
  BoxAtt->SetDaughtersInvisible(true);
  BoxAtt->SetForceSolid(true);
  BoxAtt->SetForceAuxEdgeVisible(true);

  //  fLogicCathode->SetVisAttributes(simpleBoxAtt);
  fLogicArapuca->SetVisAttributes(simpleBoxAttKGM);
  fLogicArapucaShort->SetVisAttributes(simpleBoxAttKGM);
  //  fLogicFCShort->SetVisAttributes(simpleBoxAtt);
  fLogicAnodeT->SetVisAttributes(simplePlain);
  fLogicAnodeB->SetVisAttributes(simplePlain);
  fLogicFC->SetVisAttributes(simpleBoxAtt);

  
  return fPhysiWorld;
}


void DetectorConstruction::ShieldingFloor()
{

  const double ht =  fht ; //m                                                                                                                                                                          
  const double zbsp = fSpacing; //m                                                                                                                                                                      
  double zpl(zbsp/2.);
  double xpl(0.);
  double eps(0.215);
  int cpIT(0), cpIB(0), cpIL(0), cpIR(0),cpBlt(0);
  double BlockThickness(0.40);
  const double BlockWidth = fSpacing-fIFlangeWaist-0.001 - 0.050;

  G4Box* ShieldBlock = new G4Box("ShieldBlock", (zbsp-0.050)/2.*m, (BlockThickness/2.0)*m, (BlockWidth/2.)*m); // 20cm for the fBPSE density 1.60, 30cm for fBP density 1.0
  G4LogicalVolume* ShieldBlockLog = new G4LogicalVolume(ShieldBlock, fBP, "ShieldBlockLog");
  G4VisAttributes* simpleShieldAtt= new G4VisAttributes(G4Colour::Blue());
  simpleShieldAtt->SetDaughtersInvisible(true);
  simpleShieldAtt->SetForceSolid(true);
  simpleShieldAtt->SetForceAuxEdgeVisible(true);
  ShieldBlockLog->SetVisAttributes(simpleShieldAtt);
  
  // Top, bottom  
  for (size_t ii=0;ii<=19 /*20*/;ii++) {
    // loop on x for top and bottom                                                                                                                                                                     
    for (int jj=-6;jj<6;jj++) {

      new G4PVPlacement(0,G4ThreeVector((jj)*zbsp*m,(-ht+eps)*m,(zpl)*m),"ShieldBot",
							 ShieldBlockLog,      //its logical volume   
							 fPhysOuterAir,           //its mother  volume
							 false,                 //no boolean operation
							 cpIB++, // copyNo
							 true); //check for overlaps
      new G4PVPlacement(0,G4ThreeVector((jj)*zbsp*m,(-ht+eps)*m,(-zpl)*m),"ShieldBot",
							 ShieldBlockLog,      //its logical volume   
							 fPhysOuterAir,           //its mother  volume
							 false,                 //no boolean operation
							 cpIB++, // copyNo
							 true); //check for overlaps
    }
    zpl+=zbsp;
  }

}
