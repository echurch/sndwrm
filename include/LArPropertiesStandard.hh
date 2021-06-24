/*
 * @file    LArPropertiesStandard.h
 * @brief   Service provider with utility LAr functions
 * @see     LArPropertiesStandard.cxx LArPropertiesStandardTestHelpers.h
 *
 * The provider detinfo::LArProperiesStandard supports simple setup for testing
 * environment, by including in your test:
 *
 *     #include "lardata/DetectorInfo/LArPropertiesStandardTestHelpers.h"
 *
 */
////////////////////////////////////////////////////////////////////////
// Historical authors:
//
// maddalena.antonello@lngs.infn.it
// ornella.palamara@lngs.infn.it
// msoderbe@syr.edu
// joshua.spitz@yale.edu
//
// Optical Properties:
// bjpjones@mit.edu
//
// Separation of service from Detector info class:
// jpaley@fnal.gov
//
// Test system:
// petrillo@fnal.gov
//
////////////////////////////////////////////////////////////////////////

#ifndef DETECTORINFO_LARPROPERTIESSTANDARD_H
#define DETECTORINFO_LARPROPERTIESSTANDARD_H



// C/C++ standard libraries
#include <string>
#include <vector>
#include <map>
#include <set>
#include <iostream>

  /**
   * @brief Properties related to liquid argon environment in the detector
   *
   * This class can access databases via DatabaseUtil service.
   *
   * @note Some of the database connection properties are established before
   * the beginning of the job and if they change this service will not be
   * aware of it. These properties petrain, so far, only the connection mode
   * and not any content of the databases themselves.
   * @note 2: the database connection features for this base class have been removed
   */
class LArPropertiesStandard { //: public LArProperties {
      public:
    LArPropertiesStandard() {};
  //    explicit LArPropertiesStandard    ();
    //LArPropertiesStandard(LArPropertiesStandard const&) = delete;
  //  virtual ~LArPropertiesStandard(){} = default;
   ~LArPropertiesStandard(){} ;

    /**
     * @brief Configures the provider
     * @param pset configuration parameter set
     * @param ignore_params unknown parameters to be tolerated (optional)
     *
     * This method will validate the parameter set (except for the parameters
     * it's explicitly told to ignore) and extract the useful information out
     * of it.
     */


    double RadiationLength()  	     const { return fRadiationLength; } ///< g/cm^2

    double Argon39DecayRate()              const { return fArgon39DecayRate; }  // decays per cm^3 per second

    /// Ar atomic number
    double AtomicNumber() const { return fZ; }
    /// Ar atomic mass (g/mol)
    double AtomicMass() const { return fA; }
    /// Ar mean excitation energy (eV)
    double ExcitationEnergy() const { return fI; }

    double ScintResolutionScale() const { return fScintResolutionScale; }
    double ScintFastTimeConst()   const { return fScintFastTimeConst;   }
    double ScintSlowTimeConst()   const { return fScintSlowTimeConst;   }
    double ScintBirksConstant()   const { return fScintBirksConstant;   }

    bool ScintByParticleType()    const { return fScintByParticleType;  }

    double ScintYield(bool prescale = false)         const { return fScintYield * ScintPreScale(prescale);}
    double ScintPreScale(bool prescale = true)       const { return (prescale ? fScintPreScale : 1);      }
    double ScintYieldRatio()                         const { return fScintYieldRatio;                     }

    double ProtonScintYield(bool prescale = false)   const { return fProtonScintYield * ScintPreScale(prescale);  }
    double ProtonScintYieldRatio()                   const { return fProtonScintYieldRatio;                       }
    double MuonScintYield(bool prescale = false)     const { return fMuonScintYield * ScintPreScale(prescale);    }
    double MuonScintYieldRatio()                     const { return fMuonScintYieldRatio;                         }
    double KaonScintYield(bool prescale = false)     const { return fKaonScintYield * ScintPreScale(prescale);    }
    double KaonScintYieldRatio()                     const { return fKaonScintYieldRatio;                         }
    double PionScintYield(bool prescale = false)     const { return fPionScintYield * ScintPreScale(prescale);    }
    double PionScintYieldRatio()                     const { return fPionScintYieldRatio;                         }
    double ElectronScintYield(bool prescale = false) const { return fElectronScintYield * ScintPreScale(prescale);}
    double ElectronScintYieldRatio()                 const { return fElectronScintYieldRatio;                     }
    double AlphaScintYield(bool prescale = false)    const { return fAlphaScintYield * ScintPreScale(prescale);   }
    double AlphaScintYieldRatio()                    const { return fAlphaScintYieldRatio;                        }
    double Ar40ScintYield(bool prescale = false)    const { return fAr40ScintYield * ScintPreScale(prescale);   }
    double Ar40ScintYieldRatio()                    const { return fAr40ScintYieldRatio;                        }
    bool CerenkovLightEnabled()                      const { return fEnableCerenkovLight;                         }


  std::map<double, double> SlowScintSpectrum() const ;
  std::map<double, double> FastScintSpectrum() const ;
  std::map<double, double> RIndexSpectrum() const ;
  std::map<double, double> AbsLengthSpectrum() const ;
  std::map<double, double> RayleighSpectrum() const ;
  std::vector<std::string> ReflectiveSurfaceNames() const {return fReflectiveSurfaceNames;} 

    std::map<std::string, std::map<double, double> > SurfaceReflectances() const;
    std::map<std::string, std::map<double, double> > SurfaceReflectanceDiffuseFractions() const;
  

    void SetRadiationLength(double rl) { fRadiationLength = rl;}
    void SetArgon39DecayRate(double r) { fArgon39DecayRate = r;}
    void SetAtomicNumber(double z) { fZ = z;}
    void SetAtomicMass(double a) { fA = a;}
    void SetMeanExcitationEnergy(double e) { fI = e;}

    void SetFastScintSpectrum(std::vector<double> s) { fFastScintSpectrum = s;}
    void SetFastScintEnergies(std::vector<double> s) { fFastScintEnergies = s;}
    void SetSlowScintSpectrum(std::vector<double> s) { fSlowScintSpectrum = s;}
    void SetSlowScintEnergies(std::vector<double> s) { fSlowScintEnergies = s;}
    void SetRIndexSpectrum(std::vector<double> s)    { fRIndexSpectrum = s;}
    void SetRIndexEnergies(std::vector<double> s)    { fRIndexEnergies = s;}
    void SetAbsLengthSpectrum(std::vector<double> s) { fAbsLengthSpectrum = s;}
    void SetAbsLengthEnergies(std::vector<double> s) { fAbsLengthEnergies = s;}
    void SetRayleighSpectrum(std::vector<double> s)  { fRayleighSpectrum = s;}
    void SetRayleighEnergies(std::vector<double> s)  { fRayleighEnergies = s;}

    void SetScintByParticleType(bool l)        { fScintByParticleType = l;}
    void SetProtonScintYield(double y)         { fProtonScintYield = y;}
    void SetProtonScintYieldRatio(double r)    { fProtonScintYieldRatio = r;}
    void SetMuonScintYield(double y)           { fMuonScintYield = y;}
    void SetMuonScintYieldRatio(double r)      { fMuonScintYieldRatio = r;}
    void SetPionScintYield(double y)           { fPionScintYield = y;}
    void SetPionScintYieldRatio(double r)      { fPionScintYieldRatio = r;}
    void SetKaonScintYield(double y)           { fKaonScintYield = y;}
    void SetKaonScintYieldRatio(double r)      { fKaonScintYieldRatio = r;}
    void SetElectronScintYield(double y)       { fElectronScintYield = y;}
    void SetElectronScintYieldRatio(double r)  { fElectronScintYieldRatio = r;}
    void SetAlphaScintYield(double y)          { fAlphaScintYield = y;}
    void SetAlphaScintYieldRatio(double r)     { fAlphaScintYieldRatio = r;}
    void SetAr40ScintYield(double y)          { fAr40ScintYield = y;}
    void SetAr40ScintYieldRatio(double r)     { fAr40ScintYieldRatio = r;}

    void SetScintYield(double y)               { fScintYield = y;}
    void SetScintPreScale(double s)            { fScintPreScale = s;}
    void SetScintResolutionScale(double r)     { fScintResolutionScale = r; }
    void SetScintFastTimeConst(double t)       { fScintFastTimeConst = t;}
    void SetScintSlowTimeConst(double t)       { fScintSlowTimeConst = t;}
    void SetScintYieldRatio(double r)          { fScintYieldRatio = r;}
    void SetScintBirksConstant(double kb)      { fScintBirksConstant = kb;}
    void SetEnableCerenkovLight(bool f)        { fEnableCerenkovLight = f; }

    void SetReflectiveSurfaceNames(std::vector<std::string> n) { fReflectiveSurfaceNames = n;}
    void SetReflectiveSurfaceEnergies(std::vector<double> e)   { fReflectiveSurfaceEnergies = e;}
    void SetReflectiveSurfaceReflectances(std::vector<std::vector<double> > r) { fReflectiveSurfaceReflectances = r;}
    void SetReflectiveSurfaceDiffuseFractions(std::vector<std::vector<double> > f) { fReflectiveSurfaceDiffuseFractions = f;}

    void SetExtraMatProperties(bool l)        { fExtraMatProperties = l;}
    bool ExtraMatProperties() const { return fExtraMatProperties; }
    double TpbTimeConstant()  const { return fTpbTimeConstant;     }

    std::map<double, double>  TpbAbs() const;
    std::map<double, double>  TpbEm() const;

    void SetTpbTimeConstant(double y)         { fTpbTimeConstant = y;}

    void SetTpbEmissionEnergies(std::vector<double> s) { fTpbEmissionEnergies = s;}
    void SetTpbEmissionSpectrum(std::vector<double> s) { fTpbEmissionSpectrum = s;}
    void SetTpbAbsorptionEnergies(std::vector<double> s) { fTpbAbsorptionEnergies = s;}
    void SetTpbAbsorptionSpectrum(std::vector<double> s) { fTpbAbsorptionSpectrum = s;}


  private:
  protected:


    bool fIsConfigured;

    double                         fRadiationLength;  ///< g/cm^2
    double                         fArgon39DecayRate; ///<  decays per cm^3 per second

    // Following parameters are for use in Bethe-Bloch formula for dE/dx.

    double fZ;                ///< Ar atomic number
    double fA;                ///< Ar atomic mass (g/mol)
    double fI;                ///< Ar mean excitation energy (eV)


    // Optical parameters for LAr

    std::vector<double> fFastScintSpectrum;
    std::vector<double> fFastScintEnergies;
    std::vector<double> fSlowScintSpectrum;
    std::vector<double> fSlowScintEnergies;
    std::vector<double> fRIndexSpectrum;
    std::vector<double> fRIndexEnergies;
    std::vector<double> fAbsLengthSpectrum;
    std::vector<double> fAbsLengthEnergies;
    std::vector<double> fRayleighSpectrum;
    std::vector<double> fRayleighEnergies;

    bool fScintByParticleType;

    double fProtonScintYield;
    double fProtonScintYieldRatio;
    double fMuonScintYield;
    double fMuonScintYieldRatio;
    double fPionScintYield;
    double fPionScintYieldRatio;
    double fKaonScintYield;
    double fKaonScintYieldRatio;
    double fElectronScintYield;
    double fElectronScintYieldRatio;
    double fAlphaScintYield;
    double fAlphaScintYieldRatio;
    double fAr40ScintYield;
    double fAr40ScintYieldRatio;

    double fScintYield;
    double fScintPreScale;
    double fScintResolutionScale;
    double fScintFastTimeConst;
    double fScintSlowTimeConst;
    double fScintYieldRatio;
    double fScintBirksConstant;

    bool fEnableCerenkovLight;

    std::vector<std::string>          fReflectiveSurfaceNames;
    std::vector<double>               fReflectiveSurfaceEnergies;
    std::vector<std::vector<double> > fReflectiveSurfaceReflectances;
    std::vector<std::vector<double> > fReflectiveSurfaceDiffuseFractions;

    bool fExtraMatProperties;
    double fTpbTimeConstant;
    std::vector<double>               fTpbEmissionEnergies;
    std::vector<double>               fTpbEmissionSpectrum;
    std::vector<double>               fTpbAbsorptionEnergies;
    std::vector<double>               fTpbAbsorptionSpectrum;


  }; // class LArPropertiesStandard

#endif // LARPROPERTIES_H
