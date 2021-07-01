////////////////////////////////////////////////////////////////////////
/// \file MaterialPropertyLoader.cxx
//
/// \author  bjpjones@mit.edu
////////////////////////////////////////////////////////////////////////
// Class to set material properties for different materials in
// the detector. Currently mainly used to set optical properties
// for LAr and other optical components
//

// TODO convert tabs into spaces

// TODO verify the inclusion list
#include "MaterialPropertyLoader.hh"

//#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"

#include "G4LogicalVolume.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4Material.hh"
#include "G4MaterialPropertiesTable.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4OpticalSurface.hh"


//#include "messagefacility/MessageLogger/MessageLogger.h"



  //----------------------------------------------
  void MaterialPropertyLoader::SetMaterialProperty(std::string Material,
                                                   std::string Property,
                                                   std::map<double, double> PropertyVector,
                                                   double Unit)
  {
    std::map<double,double> PropVectorWithUnit;
    for(std::map<double,double>::const_iterator it=PropertyVector.begin();
        it!=PropertyVector.end();
        it++)
      {
        PropVectorWithUnit[it->first*CLHEP::eV]=it->second*Unit;
      }
    fPropertyList[Material][Property]=PropVectorWithUnit;
    // replace with MF_LOGDEBUG()
    std::cout << "MaterialPropertyLoader" <<"Added property "
	      << Material<< "  "
	      << Property << std::endl;
  }

  //----------------------------------------------
  void MaterialPropertyLoader::SetMaterialConstProperty(std::string Material,
                                                        std::string Property,
                                                        double PropertyValue,
                                                        double Unit)
  {
    fConstPropertyList[Material][Property]=PropertyValue*Unit;
    // replace with MF_LOGDEBUG()
    std::cout << "MaterialPropertyLoader" << "Added const property "
	      << Material << "  "
	      << Property << " = " << PropertyValue << std::endl;
  }

  //----------------------------------------------
  void MaterialPropertyLoader::SetBirksConstant(std::string Material,
                                                double PropertyValue,
                                                double Unit)
  {
    fBirksConstants[Material]=PropertyValue*Unit;
    // replace with MF_LOGDEBUG()
    std::cout << "MaterialPropertyLoader" << "Set Birks constant "
	      << Material << std::endl;
  }

  //----------------------------------------------
  void MaterialPropertyLoader::UpdateGeometry(G4LogicalVolumeStore * lvs)
  {
    std::map<std::string,G4MaterialPropertiesTable*> MaterialTables;
    std::map<std::string,bool> MaterialsSet;

    // TODO replace console output with messagefacility output
    std::cout <<"MaterialPropertyLoader" << "UPDATING GEOMETRY" << std::endl;

    // Loop over each material with a property vector and create a new material table for it
    for(std::map<std::string,std::map<std::string,std::map<double,double> > >::const_iterator i=fPropertyList.begin(); i!=fPropertyList.end(); i++){
      std::string Material=i->first;
      MaterialsSet[Material]=true;

      std::cout << "MatPropLoad:UpdateGeom(): building new G4MatPropTable for Material: " << Material << std::endl;
      MaterialTables[Material]=new G4MaterialPropertiesTable;
    }

    // Loop over each material with a const property,
    // if material table does not exist, create one
    std::cout << "MatPropLoad:UpdateGeom(): fPropList.size = " << fPropertyList.size() << std::endl;
    std::cout << "MatPropLoad:UpdateGeom(): fConstPropList.size = " << fConstPropertyList.size() << std::endl;
    for(std::map<std::string,std::map<std::string,double> >::const_iterator i=fConstPropertyList.begin(); i!=fConstPropertyList.end(); i++){
      std::string Material=i->first;
      if(!MaterialsSet[Material]){
        MaterialsSet[Material]=true;
	std::cout << "MatPropLoad:UpdateGeom(): building new Const G4MatPropTable for Material: " << Material << std::endl;
        MaterialTables[Material]=new G4MaterialPropertiesTable;
      }
    }

    // For each property vector, convert to an array of g4doubles and
    // feed to materials table Lots of firsts and seconds!  See annotation
    // in MaterialPropertyLoader.h to follow what each element is

    for(std::map<std::string,std::map<std::string,std::map<double,double> > >::const_iterator i=fPropertyList.begin(); i!=fPropertyList.end(); i++){
      std::string Material=i->first;
      for(std::map<std::string,std::map<double,double> >::const_iterator j = i->second.begin(); j!=i->second.end(); j++){
        std::string Property=j->first;
        std::vector<G4double> g4MomentumVector;
        std::vector<G4double> g4PropertyVector;

        for(std::map<double,double>::const_iterator k=j->second.begin(); k!=j->second.end(); k++){
          g4MomentumVector.push_back(k->first);
          g4PropertyVector.push_back(k->second);
        }
        int NoOfElements=g4MomentumVector.size();
        MaterialTables[Material]->AddProperty(Property.c_str(),&g4MomentumVector[0], &g4PropertyVector[0],NoOfElements);
        // replace with mf::LogVerbatim()
	std::cout << "MaterialPropertyLoader" << "Added property "
		  <<Property 
		  <<" to material table "
		  << Material << std::endl;
      }
    }

    //Add each const property element
    for(std::map<std::string,std::map<std::string,double > >::const_iterator i = fConstPropertyList.begin(); i!=fConstPropertyList.end(); i++){
      std::string Material=i->first;
      for(std::map<std::string,double>::const_iterator j = i->second.begin(); j!=i->second.end(); j++){
        std::string Property=j->first;
        G4double PropertyValue=j->second;
        MaterialTables[Material]->AddConstProperty(Property.c_str(), PropertyValue);
        // replace with mf::LogVerbatim()
	std::cout << "MaterialPropertyLoader" << "Added const property "
		  <<Property<<" = " << PropertyValue
		  <<" to material table "
		  << Material << std::endl;
      }
    }

    //Loop through geometry elements and apply relevant material table where materials match
    for ( G4LogicalVolumeStore::iterator i = lvs->begin(); i != lvs->end(); ++i ){
      G4LogicalVolume* volume = (*i);
      G4Material* TheMaterial = volume->GetMaterial();
      std::string Material = TheMaterial->GetName();
      std::string volLower = volume->GetName();


      //
      // create reflective surfaces corresponding to the volumes made of some
      // selected materials
      //

      //--------------------------> FIXME <-----------------(parameters from fcl files(?))

      G4MaterialPropertyVector* PropertyPointer = 0;
      if(MaterialTables[Material])
        PropertyPointer = MaterialTables[Material]->GetProperty("REFLECTIVITY");

      if(Material=="Copper"){
        std::cout<< "copper foil surface set "<<volume->GetName()<<std::endl;
        if(PropertyPointer) {
          std::cout<< "defining Copper optical boundary "<<std::endl;
          G4OpticalSurface* refl_opsurfc = new G4OpticalSurface("Surface copper",glisur,ground,dielectric_metal);
          refl_opsurfc->SetMaterialPropertiesTable(MaterialTables[Material]);
          refl_opsurfc->SetPolish(0.2);
          new G4LogicalSkinSurface("refl_surfacec",volume, refl_opsurfc);
        }
        else
          std::cout<< "Warning: Copper surface in the geometry without REFLECTIVITY assigned"<<std::endl;
      }


      /* 
	 Note that I have a G10 clause below. EC, 12-May-2021.
      if(Material=="G10"){
        std::cout<< "G10 surface set "<<volume->GetName()<<std::endl;
        if(PropertyPointer) {
          std::cout<< "defining G10 optical boundary "<<std::endl;
          G4OpticalSurface* refl_opsurfg = new G4OpticalSurface("g10 Surface",glisur,ground,dielectric_metal);
          refl_opsurfg->SetMaterialPropertiesTable(MaterialTables[Material]);
          refl_opsurfg->SetPolish(0.1);
          new G4LogicalSkinSurface("refl_surfaceg",volume, refl_opsurfg);
        }
        else
          std::cout<< "Warning: G10 surface in the geometry without REFLECTIVITY assigned"<<std::endl;
      }
      */

      if(Material=="vm2000"){
        std::cout<< "vm2000 surface set "<<volume->GetName()<<std::endl;
        if(PropertyPointer) {
          std::cout<< "defining vm2000 optical boundary "<<std::endl;
          G4OpticalSurface* refl_opsurf = new G4OpticalSurface("Reflector Surface",unified,groundfrontpainted,dielectric_dielectric);
          refl_opsurf->SetMaterialPropertiesTable(MaterialTables[Material]);
          G4double sigma_alpha = 0.8;
          refl_opsurf->SetSigmaAlpha(sigma_alpha);
          new G4LogicalSkinSurface("refl_surface",volume, refl_opsurf);
        }
        else
          std::cout<< "Warning: vm2000 surface in the geometry without REFLECTIVITY assigned"<<std::endl;
      }
      if(Material=="ALUMINUM_Al"){
        std::cout<< "ALUMINUM_Al surface set "<<volume->GetName()<<std::endl;
        if(PropertyPointer) {
          std::cout<< "defining ALUMINUM_Al optical boundary "<<std::endl;
          G4OpticalSurface* refl_opsurfs = new G4OpticalSurface("Surface Aluminum",glisur,ground,dielectric_metal);
          refl_opsurfs->SetMaterialPropertiesTable(MaterialTables[Material]);
          refl_opsurfs->SetPolish(0.5);
          new G4LogicalSkinSurface("refl_surfaces",volume, refl_opsurfs);
        }
        else
          std::cout<< "Warning: ALUMINUM_Al surface in the geometry without REFLECTIVITY assigned"<<std::endl;
      }
      if(Material=="STEEL_STAINLESS_Fe7Cr2Ni"){
        std::cout<< "STEEL_STAINLESS_Fe7Cr2Ni surface set "<<volume->GetName()<<std::endl;
        if(PropertyPointer) {
          std::cout<< "defining STEEL_STAINLESS_Fe7Cr2Ni optical boundary "<<std::endl;
          G4OpticalSurface* refl_opsurfs = new G4OpticalSurface("Surface Steel",glisur,ground,dielectric_metal);
          refl_opsurfs->SetMaterialPropertiesTable(MaterialTables[Material]);
          refl_opsurfs->SetPolish(0.5);
          new G4LogicalSkinSurface("refl_surfaces",volume, refl_opsurfs);
        }
        else
          std::cout<< "Warning: STEEL_STAINLESS_Fe7Cr2Ni surface in the geometry without REFLECTIVITY assigned"<<std::endl;
      }
      
      // EC: adding Acrylic
      if(Material=="Acrylic"){
        std::cout<< "Acrylic surface set "<<volume->GetName()<<std::endl;
	G4MaterialPropertyVector* PropertyPointer2 = 0;
        PropertyPointer2 = MaterialTables[Material]->GetProperty("REFLECTANCE_Acrylic");
        if(PropertyPointer || PropertyPointer2 ) {
          std::cout<< "defining Acrylic optical boundary "<<std::endl;
          G4OpticalSurface* refl_opsurfs = new G4OpticalSurface("Surface Acrylic",glisur,ground,dielectric_metal); //dielectric_dielectric);
          refl_opsurfs->SetMaterialPropertiesTable(MaterialTables[Material]);
          refl_opsurfs->SetPolish(0.95);
          new G4LogicalSkinSurface("refl_surfaces",volume, refl_opsurfs);
        }
        else
          std::cout<< "Warning: Acrylic surface in the geometry without REFLECTIVITY assigned"<<std::endl;
      }      // EC: adding G10



      if(Material=="G10"){
        std::cout<< "G10 surface set "<<volume->GetName()<<std::endl;
	G4MaterialPropertyVector* PropertyPointer2 = 0;
        PropertyPointer2 = MaterialTables[Material]->GetProperty("REFLECTANCE_G10");
        if(PropertyPointer || PropertyPointer2 ) {
          std::cout<< "defining G10 optical boundary "<<std::endl;
          G4OpticalSurface* refl_opsurfs = new G4OpticalSurface("Surface G10",glisur,ground,dielectric_metal); //dielectric_metal);
          refl_opsurfs->SetMaterialPropertiesTable(MaterialTables[Material]);
          refl_opsurfs->SetPolish(0.95);
          new G4LogicalSkinSurface("refl_surfaces",volume, refl_opsurfs);
        }
        else
          std::cout<< "Warning: G10 surface in the geometry without REFLECTIVITY assigned"<<std::endl;
      }
      //-----------------------------------------------------------------------------

      // EC: adding SiPM
      if(Material=="SiPM"){
        std::cout<< "SiPM surface set "<<volume->GetName()<<std::endl;
	G4MaterialPropertyVector* PropertyPointer3 = 0;
        PropertyPointer3 = MaterialTables[Material]->GetProperty("REFLECTANCE_SiPM");
        if(PropertyPointer || PropertyPointer3 ) {
          std::cout<< "defining SiPM optical boundary "<<std::endl;
          G4OpticalSurface* refl_opsurfs = new G4OpticalSurface("Surface SiPM",glisur,ground,dielectric_dielectric);  
          refl_opsurfs->SetMaterialPropertiesTable(MaterialTables[Material]);
          refl_opsurfs->SetPolish(0.0); 
          new G4LogicalSkinSurface("refl_surfaces",volume, refl_opsurfs);
        }
        else
          std::cout<< "Warning: SiPM surface in the geometry without REFLECTIVITY assigned"<<std::endl;
      }
      //-----------------------------------------------------------------------------

      //
      // apply the remaining material properties
      //
      for(std::map<std::string,G4MaterialPropertiesTable*>::const_iterator j=MaterialTables.begin(); j!=MaterialTables.end(); j++){
        if(Material==j->first){
          TheMaterial->SetMaterialPropertiesTable(j->second);
	  std::cout << "MLP::UpdateGeometry(): Material/Pointer" << Material << "/" << TheMaterial << std::endl;
          //Birks Constant, for some reason, must be set separately
          if(fBirksConstants[Material]!=0) {
	    std::cout << "MLP::UpdateGeometry(): Setting Birks Constant: " << fBirksConstants[Material]  << std::endl;
            TheMaterial->GetIonisation()->SetBirksConstant(fBirksConstants[Material]);
	  }
          volume->SetMaterial(TheMaterial);
        }
      }
    }
  }


void MaterialPropertyLoader::SetReflectances(std::string Material, std::map<std::string,std::map<double, double> > Reflectances,  std::map<std::string,std::map<double, double> >  DiffuseFractions)
  {
    std::map<double, double> ReflectanceToStore;
    std::map<double, double> DiffuseToStore;

    for(std::map<std::string,std::map<double,double> >::const_iterator itMat=Reflectances.begin();
        itMat!=Reflectances.end();
        ++itMat)
      {
        std::string ReflectancePropName = std::string("REFLECTANCE_") + itMat->first;
        ReflectanceToStore.clear();
        for(std::map<double,double>::const_iterator itEn=itMat->second.begin();
            itEn!=itMat->second.end();
            ++itEn)
          {
            ReflectanceToStore[itEn->first]=itEn->second;
          }
        SetMaterialProperty(Material, ReflectancePropName, ReflectanceToStore,1);
      }

    for(std::map<std::string,std::map<double,double> >::const_iterator itMat=DiffuseFractions.begin();
        itMat!=DiffuseFractions.end();
        ++itMat)
      {
        std::string DiffusePropName = std::string("DIFFUSE_REFLECTANCE_FRACTION_") + itMat->first;
        DiffuseToStore.clear();
        for(std::map<double,double>::const_iterator itEn=itMat->second.begin();
            itEn!=itMat->second.end();
            ++itEn)
          {
            DiffuseToStore[itEn->first]=itEn->second;
          }
        SetMaterialProperty(Material, DiffusePropName, DiffuseToStore,1);
      }

  }

// To work over all Materials for which we've set ReflectanceFractions/DiffuseFractions EC, 11-Feb-2021.
void MaterialPropertyLoader::SetReflectances(std::map<std::string,std::map<double, double> > Reflectances,  std::map<std::string,std::map<double, double> >  DiffuseFractions)
  {
    std::map<double, double> ReflectanceToStore;
    std::map<double, double> DiffuseToStore;

    for(std::map<std::string,std::map<double,double> >::const_iterator itMat=Reflectances.begin();
        itMat!=Reflectances.end();
        ++itMat)
      {
        std::string ReflectancePropName = std::string("REFLECTANCE_") + itMat->first;
        ReflectanceToStore.clear();
        for(std::map<double,double>::const_iterator itEn=itMat->second.begin();
            itEn!=itMat->second.end();
            ++itEn)
          {
            ReflectanceToStore[itEn->first]=itEn->second;
          }
        SetMaterialProperty(itMat->first, ReflectancePropName, ReflectanceToStore,1);
      }

    for(std::map<std::string,std::map<double,double> >::const_iterator itMat=DiffuseFractions.begin();
        itMat!=DiffuseFractions.end();
        ++itMat)
      {
        std::string DiffusePropName = std::string("DIFFUSE_REFLECTANCE_FRACTION_") + itMat->first;
        DiffuseToStore.clear();
        for(std::map<double,double>::const_iterator itEn=itMat->second.begin();
            itEn!=itMat->second.end();
            ++itEn)
          {
            DiffuseToStore[itEn->first]=itEn->second;
          }
        SetMaterialProperty(itMat->first, DiffusePropName, DiffuseToStore,1);
      }

  }


  void MaterialPropertyLoader::SetReflectances(std::map<std::string,std::map<double, double> > Reflectances)
  {
    std::map<double, double> ReflectanceToStore;

    for(std::map<std::string,std::map<double,double> >::const_iterator itMat=Reflectances.begin();
        itMat!=Reflectances.end();
        ++itMat)
      {
        ReflectanceToStore.clear();
        for(std::map<double,double>::const_iterator itEn=itMat->second.begin();
            itEn!=itMat->second.end();
            ++itEn)
          {
            ReflectanceToStore[itEn->first]=itEn->second;
          }
        SetMaterialProperty(itMat->first, "REFLECTIVITY", ReflectanceToStore,1);
      }
  }

 

  void MaterialPropertyLoader::SetPropertiesFromServices()
  {
    // Ripped out of LArSoft here, we have no service to read the fhicl parameters that ordinarily will have set these properties.
    // Hence we must hand-set them first. EC, 5-Jan-2021.
    // Start by cutnpasting in /dune/products/lardataalg/v08_09_02/job/larproperties.fcl
    // Perhaps this should really be done with a G4Messenger and the .mac run-time file...
    


//For following parameters, see http://pdg.lbl.gov/AtomicNuclearProperties/

    LarProp = new LArPropertiesStandard();

    double RadiationLength( 19.55);   //g/cm^2
    LarProp->SetRadiationLength(RadiationLength);


    double AtomicNumber(     18.);      //Ar atomic number.
    LarProp->SetAtomicNumber(AtomicNumber);
    double AtomicMass(     39.948);      //Ar atomic number.
    LarProp->SetAtomicNumber(AtomicMass);


    double ExcitationEnergy (188.0);   //Ar mean excitation energy (eV).
    LarProp->SetMeanExcitationEnergy(ExcitationEnergy);
    double Argon39DecayRate( 0.00141); //decays per cm^3 per second.  Assumes 1.01 Bq/kg and a density of 1.396 g/cc
    LarProp->SetArgon39DecayRate(Argon39DecayRate);

      // Optical properties
//Fast and slow scintillation emission spectra, from [J Chem Phys vol 91 (1989) 1469]
    std::vector<double> FastScintEnergies { 6.0,  6.7,  7.1,  7.4,  7.7, 7.9,  8.1,  8.4,  8.5,  8.6,  8.8,  9.0,  9.1,  9.4,  9.8,  10.4,  10.7};
    LarProp->SetFastScintEnergies(FastScintEnergies);
    std::vector<double> SlowScintEnergies { 6.0,  6.7,  7.1,  7.4,  7.7, 7.9,  8.1,  8.4,  8.5,  8.6,  8.8,  9.0,  9.1,  9.4,  9.8,  10.4,  10.7};
    LarProp->SetSlowScintEnergies(SlowScintEnergies);
    std::vector<double> FastScintSpectrumloc { 0.0,  0.04, 0.12, 0.27, 0.44, 0.62, 0.80, 0.91, 0.92, 0.85, 0.70, 0.50, 0.31, 0.13, 0.04,  0.01, 0.0};
    LarProp->SetFastScintSpectrum(FastScintSpectrumloc);
    std::vector<double> SlowScintSpectrumloc { 0.0,  0.04, 0.12, 0.27, 0.44, 0.62, 0.80, 0.91, 0.92, 0.85, 0.70, 0.50, 0.31, 0.13, 0.04,  0.01, 0.0};
    LarProp->SetSlowScintSpectrum(SlowScintSpectrumloc);
  
    double ScintResolutionScale ( 1.);     //resolution factor used by G4 scintillation
    LarProp->SetScintResolutionScale(ScintResolutionScale);
    double ScintFastTimeConst (  6.);     //fast scintillation time constant (ns)
    LarProp->SetScintFastTimeConst(ScintFastTimeConst);
    double ScintSlowTimeConst (  1590.);  //slow scintillation time constant (ns)
    LarProp->SetScintSlowTimeConst( ScintSlowTimeConst );
    double ScintBirksConstant (  0.069);  //birks constant for LAr (1/MeV cm)
    LarProp->SetScintBirksConstant (ScintBirksConstant);
    double ScintYield       ( 24000.); //total scintillation yield (ph/Mev)         
    LarProp->SetScintYield(ScintYield);
    //    double ScintPreScale      (0.03);   //Scale factor to reduce number of photons simulated
    double ScintPreScale      (1.);   //Scale factor to reduce number of photons simulated
    LarProp->SetScintPreScale(ScintPreScale);
                              //Later QE should be corrected for this scale
    double ScintYieldRatio   (  0.3);    //fast / slow scint ratio (needs revisitting)
    LarProp->SetScintYieldRatio(ScintYieldRatio);

                              //quenching per particle in fast op sim
    bool EnableCerenkovLight( true);    //whether to switch on cerenkov light (slow)
    LarProp->SetEnableCerenkovLight(EnableCerenkovLight);

    bool ScintByParticleType(true);
    //    bool ScintByParticleType(false);
    LarProp->SetScintByParticleType(ScintByParticleType);
 //Scintillation yields and fast/slow ratios per particle type 
    double MuonScintYield          (24000);
    LarProp->SetMuonScintYield(MuonScintYield);
    double MuonScintYieldRatio     (0.23);
    LarProp->SetMuonScintYieldRatio(MuonScintYieldRatio);
    double PionScintYield          (24000);
    LarProp->SetPionScintYield(PionScintYield);
    double PionScintYieldRatio     (0.23 );
    LarProp->SetPionScintYieldRatio(PionScintYieldRatio);
    double ElectronScintYield      (20000);
    LarProp->SetElectronScintYield (ElectronScintYield );
    double ElectronScintYieldRatio (0.27);
    LarProp->SetElectronScintYieldRatio(ElectronScintYieldRatio);
    double KaonScintYield          (24000);
    LarProp->SetKaonScintYield (KaonScintYield );
    double KaonScintYieldRatio     (0.23);
    LarProp->SetKaonScintYieldRatio(KaonScintYieldRatio);
    double ProtonScintYield        (19200);
    LarProp->SetProtonScintYield(ProtonScintYield);
    double ProtonScintYieldRatio   (0.29);
    LarProp->SetProtonScintYieldRatio(ProtonScintYieldRatio);
    double AlphaScintYield         (16800);
    LarProp->SetAlphaScintYield(AlphaScintYield);
    double AlphaScintYieldRatio    (0.56);
    LarProp->SetAlphaScintYieldRatio(AlphaScintYieldRatio);
    double Ar40ScintYield         (16800);
    LarProp->SetAr40ScintYield(Ar40ScintYield);
    double Ar40ScintYieldRatio    (0.90); // EC made-up!
    LarProp->SetAr40ScintYieldRatio(Ar40ScintYieldRatio);


//Refractive index as a function of energy (eV) from arXiv:1502.04213v1
  std::vector<double> RIndexEnergies { 1.900,  2.934,  3.592,  5.566,  6.694,  7.540,  8.574,  9.044,  9.232,  9.420,  9.514,  9.608,  9.702,  9.796,  9.890,  9.984,  10.08,  10.27,  10.45,  10.74,  10.92 };
  LarProp->SetRIndexEnergies(RIndexEnergies);
  std::vector<double> RIndexSpectrum { 1.232,  1.236,  1.240,  1.261,  1.282,  1.306,  1.353,  1.387,  1.404,  1.423,  1.434,  1.446,  1.459,  1.473,  1.488,  1.505,  1.524,  1.569,  1.627,  1.751,  1.879 };
  LarProp->SetRIndexSpectrum(RIndexSpectrum);

 //absorption length as function of energy
  std::vector<double> AbsLengthEnergies { 4,     5,     6,     7,     8,     9,     10,    11   };
  LarProp->SetAbsLengthEnergies(AbsLengthEnergies);
  std::vector<double> AbsLengthSpectrum { 5000., 5000., 5000., 5000., 5000., 5000., 5000., 5000.};  // EC, 11-May-2021. Change back to 2000?!
  LarProp->SetAbsLengthSpectrum(AbsLengthSpectrum);

//Rayleigh scattering length (cm) @ 90K as a function of energy (eV) from arXiv:1502.04213
  std::vector<double> RayleighEnergies {   2.80,   3.00,   3.50,   4.00,  5.00,  6.00,  7.00,  8.00,  8.50,  9.00,  9.20,  9.40,  9.50,  9.60,  9.70,  9.80,  9.90,  10.0,  10.2,  10.4,  10.6, 10.8 };
  LarProp->SetRayleighEnergies(RayleighEnergies);
  std::vector<double> RayleighSpectrum { 47923., 35981., 18825., 10653., 3972., 1681., 750.9, 334.7, 216.8, 135.0, 109.7, 88.06, 78.32, 69.34, 61.06, 53.46, 46.50, 40.13, 28.91, 19.81, 12.61, 7.20 };
  LarProp->SetRayleighSpectrum(RayleighSpectrum);

 //Surface reflectivity data - vector of energy spectrum per
 //  surface type
  std::vector<double> ReflectiveSurfaceEnergies {  7, 9, 10 };
  LarProp->SetReflectiveSurfaceEnergies(ReflectiveSurfaceEnergies);
  std::vector<std::string> ReflectiveSurfaceNames {  "Acrylic" , "G10", "SiPM"};
  LarProp->SetReflectiveSurfaceNames (ReflectiveSurfaceNames );

  // 44% is the fraction of G10 vs 56% holes in the VD PCB top layer.
  std::vector<std::vector<double>> ReflectiveSurfaceReflectances   { {1., 1., 1. }, {1., 1., 1.0 }, {0., 0., 0.} }; // 0.44

  std::cout << "Replacing default Acrylic ReflectiveSurfaceReflectances " << ReflectiveSurfaceReflectances.at(0).at(0) << " with " << GetMaterialAcrylicSpecRef() << std::endl;
  ReflectiveSurfaceReflectances.at(0).at(0) = GetMaterialAcrylicSpecRef();
  ReflectiveSurfaceReflectances.at(0).at(1) = GetMaterialAcrylicSpecRef();
  ReflectiveSurfaceReflectances.at(0).at(2) = GetMaterialAcrylicSpecRef();
  std::cout << "Replacing default G10 ReflectiveSurfaceReflectances " << ReflectiveSurfaceReflectances.at(1).at(0) << " with " << GetMaterialG10SpecRef() << std::endl;
  ReflectiveSurfaceReflectances.at(1).at(0) = GetMaterialG10SpecRef();
  ReflectiveSurfaceReflectances.at(1).at(1) = GetMaterialG10SpecRef();
  ReflectiveSurfaceReflectances.at(1).at(2) = GetMaterialG10SpecRef();
  std::cout << "Replacing default LAr Absorption Length " << AbsLengthSpectrum.at(0) << " (0th element) with flat value " << GetMaterialLArAbsLength() << " [m]" << std::endl;  

  for (auto &it : AbsLengthSpectrum) 
    {it = GetMaterialLArAbsLength();}
  std::cout << "New LArAbsLength vector" ;
  for (auto &it : AbsLengthSpectrum) 
    {std::cout << it << ", "; }
  std::cout << std::endl;


  LarProp->SetAbsLengthSpectrum(AbsLengthSpectrum);

  LarProp->SetReflectiveSurfaceReflectances(ReflectiveSurfaceReflectances);
  std::vector<std::vector<double>> ReflectiveSurfaceDiffuseFractions { { 0.5,  0.5,  0.5  }, { 0.,  0.,  0.  }, {0., 0., 0.} };
  LarProp->SetReflectiveSurfaceDiffuseFractions(ReflectiveSurfaceDiffuseFractions);

//Information related with the simulation of the Wavelength Shifter (TPB) 
  bool LoadExtraMatProperties(true);    
  LarProp->SetExtraMatProperties(LoadExtraMatProperties);

 //TPB - WLS
  double TpbTimeConstant (2.5); // wls time constant in s J. Lumin 81(1999) 285
  LarProp->SetTpbTimeConstant(TpbTimeConstant);

  //WLS - TPB properties original tpb [0.0, 0.0, 0.0, 0.0588,0.235, 0.853, 1.0,1.0,0.9259,0.704,0.0296,0.011, 0.0,0.0, 0.]
  std::vector<double> TpbEmissionEnergies {0.05,1.0,1.5, 2.25, 2.481, 2.819, 2.952,2.988,3.024, 3.1, 3.14,3.1807, 3.54, 5.5, 50.39};
  LarProp->SetTpbEmissionEnergies(TpbEmissionEnergies);
  std::vector<double> TpbEmissionSpectrum {0.0, 0.0, 0.0, 0.0588,0.235, 0.853, 1.0,1.0,0.9259,0.704,0.0296,0.011, 0.0,0.0, 0.};
  LarProp->SetTpbEmissionSpectrum(TpbEmissionSpectrum);
  std::vector<double> TpbAbsorptionEnergies {0.05,1.77,2.0675, 7.42, 7.75, 8.16, 8.73, 9.78,10.69, 50.39};
  LarProp->SetTpbAbsorptionEnergies(TpbAbsorptionEnergies);
  std::vector<double> TpbAbsorptionSpectrum {100000.0,100000.0, 100000.0,0.001,0.00000000001,0.00000000001, 0.00000000001, 0.00000000001, 0.00000000001, 0.00000000001};
  LarProp->SetTpbAbsorptionSpectrum(TpbAbsorptionSpectrum);

  }


  void MaterialPropertyLoader::GetPropertiesFromServices()
  {


    //    const detinfo::DetectorProperties* DetProp = lar::providerFrom<detinfo::DetectorPropertiesService>();

    // wavelength dependent quantities

    SetMaterialProperty( "Argon", "FASTCOMPONENT", LarProp->FastScintSpectrum(), 1  );
    SetMaterialProperty( "Argon", "SLOWCOMPONENT", LarProp->SlowScintSpectrum(), 1  );
    SetMaterialProperty( "Argon", "RINDEX",        LarProp->RIndexSpectrum(),    1  );
    SetMaterialProperty( "Argon", "ABSLENGTH",     LarProp->AbsLengthSpectrum(), CLHEP::cm );
    SetMaterialProperty( "Argon", "RAYLEIGH",      LarProp->RayleighSpectrum(),  CLHEP::cm );

    // Just use Argon optphoton properties for Acrylic
    SetMaterialProperty( "Acrylic", "RINDEX",        LarProp->RIndexSpectrum(),    1  );
    SetMaterialProperty( "Acrylic", "ABSLENGTH",     LarProp->AbsLengthSpectrum(), CLHEP::cm );
    SetMaterialProperty( "Acrylic", "RAYLEIGH",      LarProp->RayleighSpectrum(),  CLHEP::cm );

    // Just use Argon optphoton properties for G10
    SetMaterialProperty( "G10", "RINDEX",        LarProp->RIndexSpectrum(),    1  );
    SetMaterialProperty( "G10", "ABSLENGTH",     LarProp->AbsLengthSpectrum(), CLHEP::cm );
    SetMaterialProperty( "G10", "RAYLEIGH",      LarProp->RayleighSpectrum(),  CLHEP::cm );

    auto VShortAbsLength = LarProp->AbsLengthSpectrum();
    for (auto& kv : VShortAbsLength)
      {
	// aMap has key value pairs kv.first and kv.second
	kv.second = 0.01 ;
      }

    // Just try to get all optphotons to enter SiPM and then die, so that the Termination Volume shows they entered. EC, 12-Feb-2021
    std::vector<double> VShortAbsLengthSpectrum { 0.00001, 0.00001, 0.00001, 0.00001, 0.000001, 0.00001, 0.00001, 0.00001};
    SetMaterialProperty( "SiPM", "RINDEX",        LarProp->RIndexSpectrum(),    1  );
    SetMaterialProperty( "SiPM", "ABSLENGTH",     VShortAbsLength, CLHEP::cm );
    SetMaterialProperty( "SiPM", "RAYLEIGH",      LarProp->RayleighSpectrum(),  CLHEP::cm );


    // scalar properties

    SetMaterialConstProperty("Argon", "SCINTILLATIONYIELD",  LarProp->ScintYield(true),       1/CLHEP::MeV ); // true = scaled down by prescale in larproperties
    SetMaterialConstProperty("Argon", "RESOLUTIONSCALE",     LarProp->ScintResolutionScale(), 1);
    SetMaterialConstProperty("Argon", "FASTTIMECONSTANT",    LarProp->ScintFastTimeConst(),   CLHEP::ns);
    SetMaterialConstProperty("Argon", "SLOWTIMECONSTANT",    LarProp->ScintSlowTimeConst(),   CLHEP::ns);
    SetMaterialConstProperty("Argon", "YIELDRATIO",          LarProp->ScintYieldRatio(),      1);
    //    SetMaterialConstProperty("Argon", "ELECTRICFIELD",       DetProp->Efield(),               CLHEP::kilovolt/CLHEP::cm);

    SetBirksConstant("Argon",LarProp->ScintBirksConstant(), CLHEP::cm/CLHEP::MeV);
    //    if(DetProp->SimpleBoundary())
    SetReflectances(LarProp->SurfaceReflectances()); // not else, but Also. EC, 12-Feb-2021.
    //    SetReflectances("Argon", LarProp->SurfaceReflectances(), LarProp->SurfaceReflectanceDiffuseFractions());
    // SetReflectances("Acrylic", LarProp->SurfaceReflectances(), LarProp->SurfaceReflectanceDiffuseFractions());
    SetReflectances(LarProp->SurfaceReflectances(), LarProp->SurfaceReflectanceDiffuseFractions());
      //    else
    //    SetReflectances(LarProp->SurfaceReflectances()); // not else, but Also. EC, 12-Feb-2021.

    // If we are using scint by particle type, load these

    if(LarProp->ScintByParticleType())
      {
        // true = scaled down by prescale in larproperties
        SetMaterialConstProperty("Argon", "PROTONSCINTILLATIONYIELD",  LarProp->ProtonScintYield(true),    1./CLHEP::MeV );
        SetMaterialConstProperty("Argon", "PROTONYIELDRATIO",          LarProp->ProtonScintYieldRatio(),   1.);
        SetMaterialConstProperty("Argon", "MUONSCINTILLATIONYIELD",    LarProp->MuonScintYield(true),      1./CLHEP::MeV );
        SetMaterialConstProperty("Argon", "MUONYIELDRATIO",            LarProp->MuonScintYieldRatio(),     1.);
        SetMaterialConstProperty("Argon", "KAONSCINTILLATIONYIELD",    LarProp->KaonScintYield(true),      1./CLHEP::MeV );
        SetMaterialConstProperty("Argon", "KAONYIELDRATIO",            LarProp->KaonScintYieldRatio(),     1.);
        SetMaterialConstProperty("Argon", "PIONSCINTILLATIONYIELD",    LarProp->PionScintYield(true),      1./CLHEP::MeV );
        SetMaterialConstProperty("Argon", "PIONYIELDRATIO",            LarProp->PionScintYieldRatio(),     1.);
        SetMaterialConstProperty("Argon", "ELECTRONSCINTILLATIONYIELD",LarProp->ElectronScintYield(true),  1./CLHEP::MeV );
        SetMaterialConstProperty("Argon", "ELECTRONYIELDRATIO",        LarProp->ElectronScintYieldRatio(), 1.);
        SetMaterialConstProperty("Argon", "ALPHASCINTILLATIONYIELD",   LarProp->AlphaScintYield(true),     1./CLHEP::MeV );
        SetMaterialConstProperty("Argon", "ALPHAYIELDRATIO",           LarProp->AlphaScintYieldRatio(),    1.);
        SetMaterialConstProperty("Argon", "AR40SCINTILLATIONYIELD",   LarProp->Ar40ScintYield(true),     1./CLHEP::MeV );
        SetMaterialConstProperty("Argon", "AR40YIELDRATIO",           LarProp->Ar40ScintYieldRatio(),    1.);
      }

    // If we are simulating the TPB load this

    if(LarProp->ExtraMatProperties())
      {
        SetMaterialProperty("TPB", "RINDEX",       LarProp->RIndexSpectrum(), 1 );
        SetMaterialProperty("TPB", "WLSABSLENGTH", LarProp->TpbAbs(),         CLHEP::m );
        SetMaterialProperty("TPB", "WLSCOMPONENT", LarProp->TpbEm(),          1 );

        SetMaterialConstProperty("TPB", "WLSTIMECONSTANT", LarProp->TpbTimeConstant(), CLHEP::ns );
      }

  }


