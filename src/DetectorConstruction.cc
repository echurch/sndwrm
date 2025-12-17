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
#include "G4RunManager.hh"
#include "G4GDMLParser.hh"

#include <string>

#include "globals.hh"
#include <CLHEP/Geometry/Transform3D.h>


DetectorConstruction::DetectorConstruction()
  
  :fDefaultMaterial(NULL),
   fPhysiWorld(NULL),fLogicWorld(NULL),fSolidWorld(NULL),
   fPhysiVol(NULL),fLogicVol(NULL),fSolidVol(NULL)
  
{

  // World/Envelope Dimensions
  fWorldSizeX = 50.0 * m; 
  fWorldSizeY = 50.0 * m; 
  fWorldSizeZ = 75.0 * m;
   
  // Cryostat Dimensions
  fCryostat_x = 14.8 * m; 
  fCryostat_y = 13.0 * m; 
  fCryostat_z = 62.0 * m; 
  
  // Field Cage (FC) Dimensions
  fFC_x = 13.5 * m; 
  fFC_y = fCryostat_y; // Retains unit from fCryostat_y
  fFC_z = 60.0*m; // Retains unit from fCryostat_z
   
  // Cathode Dimensions
  fCathode_x = 13.5 * m; 
  fCathode_z = 60.0 * m; 
  
  // APA internal size (LatY was commented as 'cm' but value is 6.5, which is large for cm. Assuming m.)
  fLatY = 6.5 * m; 
  fLatZ = 60.0 * m;
  
  // Thickness
  fthickness = 0.10 * m; 
  
  // Field Cage Outer profiles (FCOut)
  fFCOut_x = 0.0345 * m; 
  fFCOut_y = 0.04 * m; 
  fFCOut_z = 3.0/2 * m; 
  
  // Arapuca window size
  fAra_x = 0.007 * m; 
  fAra_y = 0.50 * m; 
  fAra_z = 0.50 * m; 
  
  // Arapuca spacing/offset
  fAra_offset = 0.07 * m; 
  fAras_yspacing = 0.005 * m; 
  
  // PTP (Wavelength Shifting) strip width
  fptp_width = 2.0 * um; // Explicitly kept in micrometers (um) as it was small
  //fptp_width = 10 * cm; // Alternative (commented out)
  
  // Vertical support bar dimensions
  fvert_bar_x = 0.075 * m; 
  // fvert_bar_y depends on fFC_y, retains unit
  fvert_bar_y = fFC_y/2.0 - 0.1 * m; 
  fvert_bar_z = 0.075 * m;
  
 // IBeams (corrected)
  fIFlangeWidth  = 0.402 * m;      // flange width (across X for top beams)
  fIFlangeThick  = 0.040 * m;      // flange thickness (along Y)
  fIFlangeWaist  = 0.022 * m;      // web (waist) thickness (across X in center box)
                                    // NOTE: fIFlangeWaist is web thickness, not used to compute flange height
  
  // Total I-beam overall depth (web height + 2*flange_thickness) is 1.108 m
  // So web (clear distance between flanges) = total_depth - 2*flange_thickness
  fIFlangeHeight = 1.108 * m - 2.0 * fIFlangeThick; // web clear height (used below as fIFlangeHeight)
  
  fITopLength    = 18.940 * m;       // top I-beam length (along Z in your top placement)
  fISideLength   = 17.840 * m;       // side beam length — DO NOT subtract cross-section dims from length
  
  // IBeam Port Holes/Locations
  fIPortHoleRad  = 0.80/2.0 * m;    // 0.8 m radius 
  
  // Keep port spacing and nominal locations; we will compute symmetric positions in IBeams()
  fISidePortLoc  = 5.907 * m;        // distance from a beam end to first port (useful to compute symmetric offsets)
  fIPortSpacing  = 4.0 * m;          // spacing between ports along the beam
  fIBotPortLoc   = 5.0 * m;          // bottom beam port location if needed
  
  // Placement distances (unchanged)
  fht = 4.0 * m; 
  fst = 16.732/2.0 * m;
  fzpl = 64.732/2.0 * m;
  fSpacing = 64.732/41.0 * m;
   

  fDetectorMessenger = new DetectorMessenger(this); // re-insert this to allow to set fFidVolume,fFloorShield
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

  
  if (!GetGDMLfile().length())
    {
      DefineMaterials();
      return ConstructLine();
    }
  
  G4GDMLParser* parser = new G4GDMLParser();
  parser->Read(GetGDMLfile(), false);  
  fPhysiWorld = parser->GetWorldVolume();
  return fPhysiWorld;
  
}


// https://indico.cern.ch/event/698002/contributions/2868259/attachments/1591642/2519098/AC_introduction_to_the_cryostat_design_warm_vessel.pdf
//

void DetectorConstruction::IBeams()
{
    //
    // === Corrected dimensions ===
    //
    const G4double flangeWidth   = fIFlangeWidth;       // 0.402 m
    const G4double flangeThick   = fIFlangeThick;       // 0.040 m
    const G4double webThick      = fIFlangeWaist;       // 0.022 m
    const G4double webHeight     = fIFlangeHeight;      // 1.028 m
    const G4double beamDepth     = flangeThick*2. + webHeight;  // = 1.108 m

    const G4double topLength  = fITopLength;    // 18.940 m
    const G4double sideLength = fISideLength - 2*beamDepth;   // 17.840 m
    const G4double holeR      = fIPortHoleRad;  // user-defined
    const G4double holeSpace  = fIPortSpacing;  // 4.0 m

    const G4double halfWebH   = webHeight / 2.;
    const G4double halfFlange = flangeThick / 2.;

    // Vertical frame half-height for placement
    const G4double fst_local = 16.732 * m / 2.;     // = 8.920 m
    const G4double sidePosX   = 17.832 * m / 2.;     // = half width
    const G4double sidePosZ   = 65.84 * m / 2. - beamDepth/2 + 25*cm;     // need ~25 cm extra since last belt is a bit different from others
    const G4double PitchIBeams = 1.6*m;
    const G4int nIBeamsLongSide = 39;
    const G4int nIBeamsShortSide = 9;
    G4double origin_z = -19*PitchIBeams ;
    G4double origin_x = -4*PitchIBeams;

    G4cout << "Corrected I-beam: webHeight=" << webHeight/m
           << " flangeThick=" << flangeThick/m
           << " totalDepth=" << beamDepth/m << G4endl;

    // === Basic solids (web + flanges) ===
    // Top/Bottom beams (long ones)
    G4Box* boxWebTop = new G4Box("IWebTop", webThick/2., webHeight/2., topLength/2.);
    G4Box* boxFlangeTop = new G4Box("IFlangeTop", flangeWidth/2., flangeThick/2., topLength/2.);

    // Side beams (vertical ones)
    G4Box* boxWebSide = new G4Box("IWebSide", webThick/2., webHeight/2., sideLength/2.);
    G4Box* boxFlangeSide = new G4Box("IFlangeSide", flangeWidth/2., flangeThick/2., sideLength/2.);

    // === Hole tube (for both top and side beams) ===
    G4Tubs* portTube = new G4Tubs("IBeamHole", 0., holeR, flangeThick/2., 0., 2.*CLHEP::pi);
    // === Hole rotation (long beam: hole axis along beam Z) ===
    G4RotationMatrix* rotHole = new G4RotationMatrix();
    rotHole->rotateY(90.*deg);

    // === Build TOP beam (with two holes) ===
    G4SubtractionSolid* topWithHole1 = new G4SubtractionSolid("TopWebHole1",
                               boxWebTop, portTube,
                               rotHole,
                               G4ThreeVector(0,0, holeSpace/2.));

    G4SubtractionSolid* topWithHole2 = new G4SubtractionSolid("TopWeb",
                               topWithHole1, portTube,
                               rotHole,
                               G4ThreeVector(0,0,-holeSpace/2.));

    // Assemble 3-piece cross-section (web + 2 flanges)
    G4MultiUnion* muTop = new G4MultiUnion("IBeamTopVol");

    muTop->AddNode(boxWebTop, HepGeom::Transform3D());
    muTop->AddNode(boxFlangeTop, HepGeom::TranslateY3D( halfWebH + halfFlange ));
    muTop->AddNode(boxFlangeTop, HepGeom::TranslateY3D(-halfWebH - halfFlange ));
    muTop->Voxelize();

    G4LogicalVolume* logTop = new G4LogicalVolume(muTop, fDUNESteel, "IBeamTopLog");

    // === Build BOTTOM beam identically without holes===
    G4MultiUnion* muBottom = new G4MultiUnion("IBeamBottomVol");

    muBottom->AddNode(topWithHole2, HepGeom::Transform3D());
    muBottom->AddNode(boxFlangeTop, HepGeom::TranslateY3D( halfWebH + halfFlange ));
    muBottom->AddNode(boxFlangeTop, HepGeom::TranslateY3D(-halfWebH - halfFlange ));
    muBottom->Voxelize();

    G4LogicalVolume* logBottom = new G4LogicalVolume(muBottom, fDUNESteel, "IBeamBottomLog");

    // === Build SIDE beams (3-hole pattern) ===
    // Three hole Z-positions:
    const G4double baseZ = sideLength/2. - 5.907*m; // original logic refined

    G4SubtractionSolid* side1 = new G4SubtractionSolid("SideWeb1", boxWebSide, portTube,
                               rotHole, G4ThreeVector(0,0, baseZ));

    G4SubtractionSolid* side2 = new G4SubtractionSolid("SideWeb2", side1, portTube,
                               rotHole, G4ThreeVector(0,0, baseZ - holeSpace));

    G4SubtractionSolid* side3 = new G4SubtractionSolid("SideWeb",  side2, portTube,
                               rotHole, G4ThreeVector(0,0, baseZ - 2.*holeSpace));

    G4MultiUnion* muSide = new G4MultiUnion("IBeamSideVol");

    muSide->AddNode(side3, HepGeom::Transform3D());
    muSide->AddNode(boxFlangeSide, HepGeom::TranslateY3D( halfWebH + halfFlange ));
    muSide->AddNode(boxFlangeSide, HepGeom::TranslateY3D(-halfWebH - halfFlange ));
    muSide->Voxelize();

    G4LogicalVolume* logSide = new G4LogicalVolume(muSide, fDUNESteel, "IBeamSideLog");

    //
    // === Vis attributes ===
    //
    auto vis = new G4VisAttributes(G4Colour::Green());
    vis->SetForceSolid(true);
    vis->SetForceAuxEdgeVisible(true);
    vis->SetDaughtersInvisible(true);

    logTop->SetVisAttributes(vis);
    logBottom->SetVisAttributes(vis);
    logSide->SetVisAttributes(vis);

    //
    // === Rotations for placement ===
    //
    G4RotationMatrix* rotTopBottom = new G4RotationMatrix();
    rotTopBottom->rotateY(90.*deg);

    G4RotationMatrix* rotLongSide = new G4RotationMatrix();  // identity
    rotLongSide->rotateZ(90.*deg); 
    rotLongSide->rotateY(-90.*deg);

    G4RotationMatrix* rotShortSide = new G4RotationMatrix();  // identity
    rotShortSide->rotateX(-90.*deg);

    G4double zpos = 0. *m;
    G4double xpos = 0. *m;

    // === Placements IBeams short side===
    for(int i = 1 ; i <= nIBeamsShortSide ; i++){
	    xpos = origin_x + (i-1)*PitchIBeams;
	    new G4PVPlacement(rotShortSide,
                      G4ThreeVector( xpos, 0., sidePosZ),
                      "IBeamZ+",
                      logSide,
                      fPhysOuterAir,
                      false, 0, true);
	    new G4PVPlacement(rotShortSide,
                      G4ThreeVector( xpos, 0., -sidePosZ),
                      "IBeamZ-",
                      logSide,
                      fPhysOuterAir,
                      false, 0, true);
    }


    // === Placements IBeams long side===
    for(int i = 1 ; i <= nIBeamsLongSide ; i++){
	    zpos = origin_z + (i-1)*PitchIBeams;
	    //Top
    	    new G4PVPlacement(rotTopBottom,
                      G4ThreeVector(0.,  fst_local, zpos),
                      "IBeamTop",
                      logTop,
                      fPhysOuterAir,
                      false, 0, true);
	    //Bottom
	    new G4PVPlacement(rotTopBottom,
                      G4ThreeVector(0., -fst_local, zpos),
                      "IBeamBottom",
                      logBottom,
                      fPhysOuterAir,
                      false, 0, true);

	    // Sides (X+ and X-)
	    new G4PVPlacement(rotLongSide,
                      G4ThreeVector( sidePosX, 0., zpos),
                      "IBeamX+",
                      logSide,
                      fPhysOuterAir,
                      false, 0, true);

	    new G4PVPlacement(rotLongSide,
                      G4ThreeVector(-sidePosX, 0., zpos),
                      "IBeamX-",
                      logSide,
                      fPhysOuterAir,
                      false, 0, true);
    }
}

void DetectorConstruction::Belts()
{
    const double clearance = 1.0*mm;

    // These must be already defined consistently:
    // fIFlangeHeightInside
    // fIFlangeWaist
    // fIFlangeThick
    // fSpacing  (Z-spacing between belts)
    // fBeltFlangeBotWidth (should be ~18.136m)
    const G4double flangeThick   = fIFlangeThick;       // 0.040 m
    const G4double flangeWidth   = fIFlangeWidth;       // 0.402 m
    const G4double webHeight     = fIFlangeHeight;      // 1.028 m
    const G4double holeR      = fIPortHoleRad;  // user-defined
    const G4double beamDepth     = flangeThick*2. + webHeight;  // = 1.108 m
    const G4double fst_local = 16.732 * m / 2.;     // = 8.920 m
    const G4double sidePosX   = 17.832 * m / 2.;     // = half width
    const G4double sidePosZ   = 65.84 * m / 2. - beamDepth/2 + 25*cm;     // need ~25 cm extra since last belt is a bit different from others
    const G4double PitchIBeams = 1.6*m;
    const G4int nIBeamsLongSide = 39;
    const G4int nIBeamsShortSide = 9;
    G4double origin_z = -19*PitchIBeams;
    G4double origin_x = -4*PitchIBeams;


    const double halfSpacingZ = PitchIBeams/2.0 - fIFlangeWaist/2.0;

    // -------------------------
    // SOLIDS
    // -------------------------

    // Belt center piece (inside I-beam gap)
    G4Box* BeltMid = new G4Box("BeltMid", fIFlangeWaist/2.0, webHeight/2.0, halfSpacingZ);
    G4Box* BeltMidTop = new G4Box("BeltMidTop", fIFlangeWaist/2.0, webHeight/3.0, halfSpacingZ);
    // Belt flanges (top and bottom)
    G4Box* BeltFlange = new G4Box("BeltFlange", flangeWidth/2.0, fIFlangeWaist/2.0, halfSpacingZ - flangeWidth/2.0);
    G4Box* BeltFlangeTop = new G4Box("BeltFlangeTop", flangeWidth/3.0, fIFlangeWaist/2.0, halfSpacingZ);
    // Hole (one side belts)
    G4Tubs* BeltPort =  new G4Tubs("BeltPortHole", 0.0, holeR, fIFlangeThick/2.0, 0.0, 2.0*CLHEP::pi);

    // Rotate hole so it's aligned along Z
    G4RotationMatrix* rotY = new G4RotationMatrix();
    rotY->rotateY(90*deg);

    G4RotationMatrix* rotZ = new G4RotationMatrix(0,0,0);
    rotZ->rotateZ(90*deg);

    G4RotationMatrix* rotX = new G4RotationMatrix(0,0,0);
    rotX->rotateX(90*deg);

    G4RotationMatrix* rotXY = new G4RotationMatrix(0,0,0);
    rotXY->rotateX(90*deg);
    rotXY->rotateY(90*deg);

    G4SubtractionSolid* BeltHole = new G4SubtractionSolid("BeltHole",
                               BeltMid,
                               BeltPort,
                               rotY,
                               G4ThreeVector());

    // -------------------------
    // TRANSFORMS FOR UNION
    // -------------------------

    HepGeom::Transform3D Tmid = HepGeom::TranslateY3D(0);
    HepGeom::Transform3D Ttop = HepGeom::TranslateY3D(+(fIFlangeHeight/2.0 + fIFlangeThick/2.0));
    HepGeom::Transform3D Tbot = HepGeom::TranslateY3D(-(fIFlangeHeight/2.0 + fIFlangeThick/2.0));
    //Top belt which is different from other sides
    HepGeom::Transform3D TFtop = HepGeom::TranslateY3D(+(webHeight/3.0 + fIFlangeThick/2.0));
    HepGeom::Transform3D TFbot = HepGeom::TranslateY3D(-(webHeight/3.0 + fIFlangeThick/2.0));

    // -------------------------
    // BUILD MULTIUNIONS
    // -------------------------
    // Belt for top (uses top-mid + small flanges)
    G4MultiUnion* BeltTop = new G4MultiUnion("BeltWithHole");
    BeltTop->AddNode(BeltMidTop, Tmid);
    BeltTop->AddNode(BeltFlangeTop, TFtop);
    BeltTop->AddNode(BeltFlangeTop, TFbot);
    BeltTop->Voxelize();
    // Belt with hole
    G4MultiUnion* BeltWithHole = new G4MultiUnion("BeltWithHole");
    BeltWithHole->AddNode(BeltHole, Tmid);
    BeltWithHole->AddNode(BeltFlange, Ttop);
    BeltWithHole->AddNode(BeltFlange, Tbot);
    BeltWithHole->Voxelize();
    // Belt without hole
    G4MultiUnion* BeltWithoutHole = new G4MultiUnion("BeltWithHole");
    BeltWithoutHole->AddNode(BeltMid, Tmid);
    BeltWithoutHole->AddNode(BeltFlange, Ttop);
    BeltWithoutHole->AddNode(BeltFlange, Tbot);
    BeltWithoutHole->Voxelize();
    // -------------------------
    // LOGICAL VOLUMES
    // -------------------------
    G4LogicalVolume *fBeltTopLog = new G4LogicalVolume(BeltTop, fDUNESteel, "BeltTopLog");
    G4LogicalVolume *fBeltWithHoleLog = new G4LogicalVolume(BeltWithHole, fDUNESteel, "BeltWithHoleLog");
    G4LogicalVolume *fBeltWithoutHoleLog = new G4LogicalVolume(BeltWithoutHole, fDUNESteel, "BeltWithHoleLog");


    // -------------------------
    // VISUALIZATION
    // -------------------------
    G4VisAttributes* beltTopVis = new G4VisAttributes(G4Colour::Green());
    G4VisAttributes* beltVis = new G4VisAttributes(G4Colour::Cyan());
    G4VisAttributes* holeVis = new G4VisAttributes(G4Colour::Brown());
    beltTopVis->SetForceSolid(true);
    beltVis->SetForceSolid(true);
    holeVis->SetForceSolid(true);

    fBeltTopLog->SetVisAttributes(beltTopVis);
    fBeltWithHoleLog->SetVisAttributes(holeVis);
    fBeltWithoutHoleLog->SetVisAttributes(beltVis);

    
    // === Placements Belts long side===
    G4double xpos = 0. *m;
    G4double ypos = 0. *m;
    G4double zpos = 0. *m;
    for(int i = 0 ; i <= nIBeamsLongSide  ; i++){
	    zpos = origin_z + (i-1)*PitchIBeams + PitchIBeams/2;
    	    for(int j = 1 ; j <= nIBeamsShortSide ; j++){
		    xpos = origin_x + (j-1)*PitchIBeams;
		    //Bottom belts
	    	    new G4PVPlacement(0,
                      G4ThreeVector(xpos,  -fst_local, zpos),
                      "IBeltBottom",
                      fBeltWithHoleLog,
                      fPhysOuterAir,
                      false, 0, true);
		    //Top belts
	    	    new G4PVPlacement(0,
                      G4ThreeVector(xpos,  fst_local, zpos),
                      "IBeltTop",
                      fBeltTopLog,
                      fPhysOuterAir,
                      false, 0, true);
	    }
	    //Lateral belts
	    for(int k = 1 ; k <= 4; k++){
		        if (k == 1) { 
				ypos = fst_local - (k-1)*(fht);
			}else{
				ypos = fst_local - (k-1)*(fht) + fht/2 - 2*holeR;
			}

			if ( (i+k) % 3 == 0 || k == 4){
				//x side +
				new G4PVPlacement(
             			rotZ,
                		G4ThreeVector(sidePosX, ypos, zpos),
        	        	"BeltX+",
	                	fBeltWithHoleLog,
                		fPhysOuterAir,
                		false,
	                	1,
        	        	true);
				//x side -
				new G4PVPlacement(
	             		rotZ,
        	        	G4ThreeVector(-sidePosX, ypos, zpos),
                		"BeltX-",
                		fBeltWithHoleLog,
                		fPhysOuterAir,
                		false,
                		1,
                		true);
			}else{
				//x side +
				new G4PVPlacement(
             			rotZ,
                		G4ThreeVector(sidePosX, ypos, zpos),
        	        	"BeltX+",
	                	fBeltWithoutHoleLog,
                		fPhysOuterAir,
                		false,
	                	1,
        	        	true);
				//x side -
				new G4PVPlacement(
	             		rotZ,
        	        	G4ThreeVector(-sidePosX, ypos, zpos),
                		"BeltX-",
                		fBeltWithoutHoleLog,
                		fPhysOuterAir,
                		false,
                		1,
                		true);

			}
		}
    }
    //Short side
    for(int i = 0 ; i <= nIBeamsShortSide + 2; i++){
		xpos = origin_x + (i -1)*PitchIBeams - PitchIBeams/2;
		for(int j = 1 ; j <= 4; j++){
		        if (j == 1) { 
				ypos = fst_local - (j-1)*(fht);
			}else{
				ypos = fst_local - (j-1)*(fht) + fht/2 - 2*holeR;
			}
			if ( (i+j) % 3 == 0 || j == 1){
				//z plus
				new G4PVPlacement(
             			rotXY,
	                	G4ThreeVector(xpos, ypos, sidePosZ),
                		"BeltZ+",
        	        	fBeltWithHoleLog,
                		fPhysOuterAir,
	                	false,
        	        	1,
                		true);
				//z minus
				new G4PVPlacement(
        	     		rotXY,
                		G4ThreeVector(xpos, ypos, -sidePosZ),
	                	"BeltZ-",
                		fBeltWithHoleLog,
        	        	fPhysOuterAir,
                		false,
                		1,
	                	true);
			}else{
				//z plus
				new G4PVPlacement(
             			rotXY,
	                	G4ThreeVector(xpos, ypos, sidePosZ),
                		"BeltZ+",
        	        	fBeltWithoutHoleLog,
                		fPhysOuterAir,
	                	false,
        	        	1,
                		true);
				//z minus
				new G4PVPlacement(
        	     		rotXY,
                		G4ThreeVector(xpos, ypos, -sidePosZ),
	                	"BeltZ-",
                		fBeltWithoutHoleLog,
        	        	fPhysOuterAir,
                		false,
                		1,
	                	true);

			}
		}
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

  G4Element* Mg  = new G4Element("Magnesium","Mg", 12, 24.3*g/mole);
  G4Element* Na  = new G4Element("Sodium","Na", 11, 22.99*g/mole);
  G4Element* Ca  = new G4Element("Calcium","Ca", 20, 40.08*g/mole);
  G4Element* S  = new G4Element("Sulphur","S",  16, 32.06*g/mole);
  
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
  fWater = man->FindOrBuildMaterial("G4_WATER");
  //  fPb = man->FindOrBuildMaterial("G4_LEAD");
  fPb = new G4Material(name="Lead", z=82., 207*g/mole,11.348*g/cm3);

  
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


  G4Material* H2O = new G4Material("Water",density= 1.0*g/cm3,ncomponents=2, kStateSolid, 293.15*kelvin);
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

  G4Material* fBP_norm = new G4Material (name="BP_norm", density= 0.95*g/cm3 , ncomponents=5, kStateSolid, 293.15*kelvin);  
  fBP_norm->AddElement(H,fractionmass=11.6*perCent);
  fBP_norm->AddElement(C,fractionmass=61.2*perCent);
  fBP_norm->AddElement(B11,fractionmass=4.0*perCent);
  fBP_norm->AddElement(B10,fractionmass=1.0*perCent);
  fBP_norm->AddElement(O,fractionmass=22.2*perCent);


  G4Material* SiO2 = new G4Material(name="SiO2", density = 2.2*g/cm3, ncomponents=2);
  SiO2->AddElement(Si,number_of_atoms=1);
  SiO2->AddElement(O,number_of_atoms=2);
  G4Material* FeO = new G4Material(name="FeO", density = 5.745*g/cm3, ncomponents=2);
  FeO->AddElement(Fe,number_of_atoms=1);
  FeO->AddElement(O,number_of_atoms=1);
  G4Material* Al2O3 = new G4Material(name="Al2O3", density = 3.97*g/cm3, ncomponents=2);
  Al2O3->AddElement(Al,number_of_atoms=2);
  Al2O3->AddElement(O,number_of_atoms=3);
  G4Material* MgO = new G4Material(name="MgO", density = 3.58 *g/cm3, ncomponents=2);
  MgO->AddElement(Mg,number_of_atoms=1);
  MgO->AddElement(O,number_of_atoms=1);
  G4Material* CO2 = new G4Material(name="CO2", density = 1.562 *g/cm3, ncomponents=2);
  CO2->AddElement(C,number_of_atoms=1);
  CO2->AddElement(O,number_of_atoms=2);
  G4Material* CaO = new G4Material(name="CaO", density = 3.35 *g/cm3, ncomponents=2);
  CaO->AddElement(Ca,number_of_atoms=1);
  CaO->AddElement(O,number_of_atoms=1);
  G4Material* Na2O = new G4Material(name="Na2O", density = 2.27 *g/cm3, ncomponents=2);
  Na2O->AddElement(Na,number_of_atoms=2);
  Na2O->AddElement(O,number_of_atoms=1);

  G4Material* DUNERock = new G4Material("DUNERock", density = 2.82 *g/cm3, ncomponents=9);
  DUNERock->AddMaterial(SiO2,fractionmass=0.5267);
  DUNERock->AddMaterial(FeO,fractionmass=0.1174);
  DUNERock->AddMaterial(Al2O3,fractionmass=0.1025);
  DUNERock->AddMaterial(MgO,fractionmass=0.0473);
  DUNERock->AddMaterial(CO2,fractionmass=0.0422);
  DUNERock->AddMaterial(CaO,fractionmass=0.0382);
  DUNERock->AddElement(C,fractionmass=0.0240);
  DUNERock->AddElement(S,fractionmass=0.0186);
  DUNERock->AddMaterial(Na2O,fractionmass=0.0053);
  fRock = DUNERock;
		       
  fBP = fBP_norm; // H2O; //fBP_norm; //fBP_SE; // Lead
  
  G4Material* acrylic = new G4Material("acrylic", 1.18*g/cm3,3); //acrylic
  acrylic->AddElement (C, 5);
  acrylic->AddElement (O, 2);
  acrylic->AddElement (H, 8);
   
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
		  fWorldSizeX/2,fWorldSizeY/2,fWorldSizeZ/2);  //its size
  
  
  fLogicWorld = new G4LogicalVolume(fSolidWorld,	//its solid
				    fRock,	//its material
				    "World");		//its name
  
  fPhysiWorld = new G4PVPlacement(0,			//no rotation
  				 G4ThreeVector(),	//at (0,0,0)
                                 "World",		//its name
                                 fLogicWorld,		//its logical volume
                                 NULL,			//its mother  volume
                                 false,			//no boolean operation
                                 0);			//copy number

  // Beginning of cryostat construction: coldskin, wood, foam, nougat, etc.*cm, ...*cm, warmskin

    // Now, need Air volume outside cryo, inside World. Otherwise this space will be G4_lAr
 
  G4Box* boxAir = new G4Box("Cavern",fCryostat_x/2 + 6*m, fCryostat_y/2 + 6*m, fCryostat_z/2 +6*m); 
  G4LogicalVolume* fLogicOAir = new G4LogicalVolume(boxAir,fmAir,"OuterAir");
  fPhysOuterAir = new G4PVPlacement(0,G4ThreeVector(0,0,0),
				    "OuterAir",
				    fLogicOAir,     //its logical volume
				    fPhysiWorld,    	//its mother  volume
				    false,			//no boolean operation
				    0, true);

				    
  

  // Create and Place I-Beams and Belts and Shielding panels.
  // All of this must have mother volume fPhysOuterAir
  IBeams();
  Belts();
  ShieldingFloor();
  //if (GetFloorShield()>0.)
  ShieldingWalls();
  //if (GetFloorShield()>0.)
  
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


  
    
  return fPhysiWorld;
}

void DetectorConstruction::ShieldingFloor()
{
    // === PARAMETERS ===
    const double ht = fht;                  // detector half-height (meters)
    const G4double fst_local = 16.732 * m / 2.;     // = 8.920 m
    const G4double webHeight     = fIFlangeHeight;      // 1.028 m
    const double eps = 0.215;               // floor clearance
    const double BlockThickness = 30*cm;   // mm -> m
    const double BlockThicknessPb = 2.5*cm;                //  
    const G4double PitchIBeams = 1.6*m;
    const G4int nIBeamsLongSide = 39;
    const G4int nIBeamsShortSide = 9;

    const int nZ = nIBeamsLongSide;         // number of rows in z
    const int nX = nIBeamsShortSide;        // number of columns in x

    // block width between beams
    const double clearance = 2. *mm; //assume there will be always a small space between shielding and beams/belts
    const double BlockWidth = PitchIBeams - fIFlangeWaist - clearance;

    std::cout << "ShieldingFloor(): ShieldBlockNeutrons thickness [m]: "
              << BlockThickness << std::endl;

    // === SOLIDS ===
    G4Box* ShieldBlockNeutron = new G4Box("ShieldBlockNeutron", BlockWidth / 2.0, BlockThickness / 2.0, BlockWidth / 2.0);

    G4Box* ShieldBlockPb = new G4Box("ShieldBlockPb", BlockWidth / 2.0, BlockThicknessPb / 2.0, BlockWidth / 2.0);

    // === LOGICAL VOLUMES ===
    G4LogicalVolume* ShieldBlockNeutronLog = new G4LogicalVolume(ShieldBlockNeutron, fBP, "ShieldBlockNeutronLog");

    G4LogicalVolume* ShieldBlockPbLog = new G4LogicalVolume(ShieldBlockPb, fPb, "ShieldBlockPbLog");

    // === VISUALIZATION ===
    auto* visShieldNeutron = new G4VisAttributes(G4Colour::Blue());
    visShieldNeutron->SetDaughtersInvisible(true);
    visShieldNeutron->SetForceSolid(true);
    visShieldNeutron->SetForceAuxEdgeVisible(true);
    ShieldBlockNeutronLog->SetVisAttributes(visShieldNeutron);
    auto* visShieldGammas = new G4VisAttributes(G4Colour::Grey());
    visShieldGammas->SetDaughtersInvisible(true);
    visShieldGammas->SetForceSolid(true);
    visShieldGammas->SetForceAuxEdgeVisible(true);
    ShieldBlockPbLog->SetVisAttributes(visShieldGammas);

    // === PLACEMENT SHIELDING===
    // Y-positions
    const double origin_z = -19 * PitchIBeams + PitchIBeams/2;    // same logic as belts
    const double origin_x = -4 * PitchIBeams + PitchIBeams/2;
    double yPbBlock  = -fst_local - webHeight/2 + BlockThicknessPb/2.0;
    double yNeutronBlock = yPbBlock + BlockThicknessPb + BlockThickness/2;
    int cpSteel = 0, cpLead = 0;

    for (int i = 0; i <= nZ; i++)
    {
        double zpos = origin_z + (i - 1) * PitchIBeams;

        for (int j = 0; j <= nX; j++)
        {
            double xpos = origin_x + (j - 1) * PitchIBeams ;

            // === BLOCKS (bottom rows) ===
            new G4PVPlacement(0,
                G4ThreeVector(xpos, yNeutronBlock, zpos),
                "ShieldBotNeutron",
                ShieldBlockNeutronLog,
                fPhysOuterAir,
                false, cpSteel++, true);

            new G4PVPlacement(0,
                G4ThreeVector(xpos, yPbBlock, zpos),
                "ShieldBotPb",
                ShieldBlockPbLog,
                fPhysOuterAir,
                false, cpSteel++, true);

        }
    }
}

void DetectorConstruction::ShieldingWalls()
{
    // === PARAMETERS ===
    const double ht = fht;                  // detector half-height (meters)
    const G4double flangeThick   = fIFlangeThick;       // 0.040 m
    const G4double fst_local = 16.732 * m / 2.;     // = 8.920 m
    const G4double webHeight     = fIFlangeHeight;      // 1.028 m
    const G4double beamDepth     = flangeThick*2. + webHeight;  // = 1.108 m
    const G4double topLength  = fITopLength;    // 18.940 m
    const double BlockThickness = 23*cm;   // mm -> m
    const double BlockHeight = 100*cm;                //  
    const double BlockWidth = 100*cm;                //  
    const G4double ContainerThickness = 2.0*mm;
    const double WaterThickness = BlockThickness - 2*ContainerThickness;
    const double WaterHeight = BlockHeight - 2*ContainerThickness;
    const double WaterWidth = BlockWidth - 2*ContainerThickness;
    const G4int nBlocksLongSide = 64;
    const G4int nBlocksShortSide = 14;
    const G4int nBlocksStack = 10;

    const int nZ = nBlocksLongSide;         // number of rows in z
    const int nX = nBlocksShortSide;        // number of columns in x
    const int nY = nBlocksStack;        // number of columns in x


    std::cout << "ShieldingFloor(): ShieldBlockNeutrons thickness [m]: "
              << BlockThickness << std::endl;

    // === SOLIDS ===
    G4Box* ShieldBlockContainer = new G4Box("ShieldBlockContainer", BlockThickness / 2.0, BlockHeight / 2.0, BlockWidth / 2.0);

    G4Box* ShieldBlockWater = new G4Box("ShieldBlockWater", WaterThickness / 2.0, WaterHeight / 2.0, WaterWidth / 2.0);

    // === LOGICAL VOLUMES ===
    G4LogicalVolume* ShieldBlockContainerLog = new G4LogicalVolume(ShieldBlockContainer, fBP, "ShieldBlockContainerLog");

    G4LogicalVolume* ShieldBlockWaterLog = new G4LogicalVolume(ShieldBlockWater, fWater, "ShieldBlockWaterLog");

    //place the water inside the container  
    new G4PVPlacement(0,
                G4ThreeVector(0, 0, 0),
                ShieldBlockWaterLog,
                "ShieldLatWater",
                ShieldBlockContainerLog,
                false, 0, true);
              
    // === VISUALIZATION ===
    auto* visShieldNeutron = new G4VisAttributes(G4Colour::Blue());
    visShieldNeutron->SetDaughtersInvisible(true);
    visShieldNeutron->SetForceSolid(true);
    visShieldNeutron->SetForceAuxEdgeVisible(true);
    ShieldBlockContainerLog->SetVisAttributes(visShieldNeutron);
    //auto* visShieldGammas = new G4VisAttributes(G4Colour::Grey());
    //visShieldGammas->SetDaughtersInvisible(true);
    //visShieldGammas->SetForceSolid(true);
    //visShieldGammas->SetForceAuxEdgeVisible(true);
    //ShieldBlockWaterLog->SetVisAttributes(visShieldGammas);

    // === PLACEMENT SHIELDING===
    // Y-positions
    const double origin_z = -nZ * BlockWidth/2 + BlockWidth/2;    // same logic as belts
    const double origin_x = -4 * BlockWidth + BlockWidth/2;
    const double origin_y = -fst_local - beamDepth + BlockHeight/2;
    const double xLatWall = topLength/2 + BlockThickness/2;
    int cpContainer = 0, cpWater = 0;

    for (int i = 0; i <= nZ; i++)
    {
        double zpos = origin_z + (i - 1) * BlockWidth;

        for (int j = 0; j <= nY; j++)
        {
            double ypos = origin_y + (j - 1) * BlockHeight ;

            // === BLOCKS (+x) ===
            new G4PVPlacement(0,
                G4ThreeVector(xLatWall, ypos, zpos),
                "ShieldLatContainerX+",
                ShieldBlockContainerLog,
                fPhysOuterAir,
                false, cpContainer++, true);

            // === BLOCKS (-x) ===
            new G4PVPlacement(0,
                G4ThreeVector(-xLatWall, ypos, zpos),
                "ShieldLatContainerX-",
                ShieldBlockContainerLog,
                fPhysOuterAir,
                false, cpContainer++, true);

        }
    }
}


void DetectorConstruction::SetFidVolume(G4ThreeVector value)
{
  fFidVol = value;
  // Do not need to re-initialize geom. This is strictly for analysis sake. EC, 8-May-2025.
  //  G4RunManager::GetRunManager()->ReinitializeGeometry();
}

void DetectorConstruction::SetFloorShield(G4double value)
{
  fFloorShield = value;
}

void DetectorConstruction::SetGDMLfile(G4String value)
{
  fGDMLfile = value;
}
