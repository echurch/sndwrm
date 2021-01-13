
#include "LArPropertiesStandard.hh"



std::map<double, double> LArPropertiesStandard::SlowScintSpectrum() const {
   std::map <double,double> m;
   for (size_t ii=0; ii<fSlowScintSpectrum.size(); ++ii) {
     m.insert(std::pair<double,double> (fSlowScintEnergies[ii],fSlowScintSpectrum[ii]));
   }
   return m;
}

std::map<double, double> LArPropertiesStandard::FastScintSpectrum() const {
   std::map <double,double> m;

   for (size_t ii=0; ii<fFastScintSpectrum.size(); ++ii) {
     //     m[fFastScintEnergies[ii]] = m[fFastScintSpectrum[ii]]; 
     m.insert(std::pair<double,double> (fFastScintEnergies[ii],fFastScintSpectrum[ii]));
   }
   return m;
}


std::map<double, double> LArPropertiesStandard::RayleighSpectrum() const {
   std::map <double,double> m;
   for (size_t ii=0; ii<fRayleighSpectrum.size(); ++ii) {
     m.insert(std::pair<double,double> (fRayleighEnergies[ii],fRayleighSpectrum[ii]));
   }
   return m;
}

std::map<double, double> LArPropertiesStandard::RIndexSpectrum() const {
   std::map <double,double> m;
   for (size_t ii=0; ii<fRIndexSpectrum.size(); ++ii) {
     m.insert(std::pair<double,double> (fRIndexEnergies[ii],fRIndexSpectrum[ii]));
   }
   return m;
}

std::map<double, double> LArPropertiesStandard::AbsLengthSpectrum() const {
   std::map <double,double> m;
   for (size_t ii=0; ii<fAbsLengthSpectrum.size(); ++ii) {
     m.insert(std::pair<double,double> (fAbsLengthEnergies[ii],fAbsLengthSpectrum[ii]));
   }
   return m;
}


std::map<std::string,std::map<double, double>> LArPropertiesStandard::SurfaceReflectances() const {
  std::map<std::string, std::map <double,double>> m;

  for (size_t ii=0; ii<fReflectiveSurfaceNames.size(); ++ii) {
    std::map<double,double> md;
    for (size_t jj=0; jj<fReflectiveSurfaceReflectances[ii].size(); ++jj) {
      md.insert(std::pair<double,double>(fReflectiveSurfaceEnergies[jj],fReflectiveSurfaceReflectances[ii][jj]));
     }
    m.insert(std::pair<std::string,std::map<double,double>>(fReflectiveSurfaceNames[ii],md));
   }
   return m;
}

std::map<std::string,std::map<double, double>> LArPropertiesStandard::SurfaceReflectanceDiffuseFractions() const {
  std::map<std::string,std::map <double,double>> m;

   for (size_t ii=0; ii<fReflectiveSurfaceNames.size(); ++ii) {
     std::map<double,double> md;
     for (size_t jj=0; jj<fReflectiveSurfaceDiffuseFractions[ii].size(); ++jj) {
      md.insert(std::pair<double,double>(fReflectiveSurfaceEnergies[jj],fReflectiveSurfaceDiffuseFractions[ii][jj]));

     }
      m.insert(std::pair<std::string,std::map<double,double>>(fReflectiveSurfaceNames[ii],md));
   }
   return m;
}


std::map<double, double> LArPropertiesStandard::TpbAbs() const {
   std::map <double,double> m;
   for (size_t ii=0; ii<fTpbAbsorptionSpectrum.size(); ++ii) {
     m.insert(std::pair<double,double>(fTpbAbsorptionEnergies[ii],fTpbAbsorptionSpectrum[ii]) );
   }
   return m;
}

std::map<double, double> LArPropertiesStandard::TpbEm() const {
   std::map <double,double> m;
   for (size_t ii=0; ii<fTpbEmissionSpectrum.size(); ++ii) {
     m.insert(std::pair<double,double>(fTpbEmissionEnergies[ii],fTpbEmissionSpectrum[ii]) );
   }
   return m;
}
