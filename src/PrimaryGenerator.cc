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
/// \file PrimaryGenerator.cc
/// \brief Implementation of the PrimaryGenerator1 class
//
//
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "PrimaryGenerator.hh"

#include "G4Event.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4PrimaryParticle.hh"
#include "G4PrimaryVertex.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGenerator::PrimaryGenerator()
  : G4VPrimaryGenerator(), fFSNeutrino(false)
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGenerator::~PrimaryGenerator()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGenerator::GeneratePrimaryVertexOpt(G4Event* event, std::vector<double> &xyzb)
{
  //vertex A uniform on a cylinder
  //

  const G4int n_particle = 1250; // from our DkMatter paper, via SCENE. 1250 photons per 100 keV n.r.

  const G4double x = xyzb.at(0)*(G4UniformRand()-0.5) ;  
  const G4double y = xyzb.at(1)*(G4UniformRand()-0.5) ;  
  const G4double z = xyzb.at(2)*(G4UniformRand()-0.5) ; 
  //
  G4ThreeVector positionA(x,y,z);

  for (int ii =0; ii< n_particle; ii++) 
    {
      G4double alpha = pi*2.*(G4UniformRand()-0.5);     //alpha uniform in (-pi, pi)
      G4double beta = pi*G4UniformRand();     //alpha uniform in (0, pi)
      G4double ux = std::cos(alpha)*std::sin(beta);
      G4double uy = std::sin(alpha)*std::sin(beta);
      G4double uz = std::cos(beta);

      G4double timeA = 0*s;
  // 
      G4PrimaryVertex* vertexA = new G4PrimaryVertex(positionA, timeA);
  
  //particle 1 at vertex A
  //
      G4ParticleDefinition* particleDefinition
	= G4ParticleTable::GetParticleTable()->FindParticle("opticalphoton");
      G4PrimaryParticle* particle1 = new G4PrimaryParticle(particleDefinition);
      particle1->SetMomentumDirection(G4ThreeVector(ux,uy,uz));    
      //      std::cout << "PrimaryGenerator: Added " << n_particle << "  opticalphoton w ux,uy,uz " << ux <<", " << uy << ", " << uz << "."  << std::endl;
      particle1->SetKineticEnergy(9.686 * eV); // 128nm
  //
      vertexA->SetPrimary(particle1);
      event->AddPrimaryVertex(vertexA);
    }

  SetParticlePosition(positionA);
  std::cout << "PrimaryGenerator: Added " << n_particle << " isotropic opticalphotons as primaries at " << x/1000. <<", " << y/1000. << ", " << z/1000. << " [m]. "  << std::endl;
}

void PrimaryGenerator::GeneratePrimaryVertexMarley(G4Event* event, std::vector<double> &xyzb, marley::Event marlev) 
{
  // Loop over each of the final particles in the MARLEY event
  const G4double x = xyzb.at(0) *2*(G4UniformRand()-0.5) ;
  const G4double y = xyzb.at(1) *2*(G4UniformRand()-0.5) ;
  const G4double z = xyzb.at(2) *2*(G4UniformRand()-0.5) ;

  const G4int n_particle = 1250; // from our DkMatter paper, via SCENE. 1250 photons per 100 keV n.r.                                                           

  G4ThreeVector positionA(x,y,z);
  G4double timeA = 0*s;
  G4PrimaryVertex* vertex = new G4PrimaryVertex(positionA, timeA);
  fFSNeutrino = false;
      
  for ( const auto& fp : marlev.get_final_particles() ) {

    // Convert each one from a marley::Particle into a G4PrimaryParticle.
    // Do this by first setting the PDG code and the 4-momentum components.

    //    std::cout << "GenPrimVertMarl:  final state particle PID is: " << fp->pdg_code() << " of energy " << fp->kinetic_energy() << std::endl;
    // If this is an Ar nuclear recoil, let's place opticalphotons on stack, instead.

    if( fp->pdg_code() == 12)
      fFSNeutrino = true;

    if (fp->pdg_code() == 1000180400) {
      int Nphotons_stack = fp->kinetic_energy() * 1000. * n_particle/100.;
      G4ParticleDefinition* particleDefinition
	= G4ParticleTable::GetParticleTable()->FindParticle("opticalphoton");

      for (int phot=0; phot<Nphotons_stack; phot++)
	{
	  G4double alpha = pi*2.*(G4UniformRand()-0.5);     //alpha uniform in (-pi, pi)
	  G4double beta = pi*G4UniformRand();     //alpha uniform in (0, pi)
	  G4double ux = std::cos(alpha)*std::sin(beta);
	  G4double uy = std::sin(alpha)*std::sin(beta);
	  G4double uz = std::cos(beta);

	  G4PrimaryParticle* particle = new G4PrimaryParticle(particleDefinition);
	  particle->SetMomentumDirection(G4ThreeVector(ux,uy,uz));
	  particle->SetKineticEnergy(9.686 * eV); // 128nm
	  vertex->SetPrimary(particle);
	}
      if (Nphotons_stack)
	std::cout << "GeneratePrimaryVertexMarley: Put " << int(Nphotons_stack) << " 9.7 eV photons on stack in lieu of Ar nucleus of KE " << fp->kinetic_energy()*1000. << " keV" << std::endl;

    }
    else {
      G4PrimaryParticle* particle = new G4PrimaryParticle( fp->pdg_code(),
							   fp->px(), fp->py(), fp->pz(), fp->total_energy() );

      // Also set the charge of the G4PrimaryParticle appropriately
      particle->SetCharge( fp->charge() );

      // Add the fully-initialized G4PrimaryParticle to the primary vertex
      vertex->SetPrimary( particle );
    }

  } // end of loop on particles

  // The primary vertex has been fully populated with all final-state particles
  // from the MARLEY event. Add it to the G4Event object so that Geant4 can
  // begin tracking the particles through the simulated geometry.

  SetParticlePosition(positionA);
  event->AddPrimaryVertex( vertex );

}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


