// source/processes/hadronic/util/include/G4Nucleus.hh

// Class Description:
// This class knows how to describe a nucleus;
// to be used in your physics implementation (not physics list) in case you need this physics.
#pragma once
#include "isotope.h"
#include "material.h"
#include "dynamicParticle.h"
#include "debug.h"
#include "type.h"
#include "threeVector.h"
#include "reactionProduct.h"
//#include "G4ReactionProductVector.hh"
#include "lorentzVector.h"

typedef struct nucleus {
	int  theA;
	int  theZ;
	dREAL aEff;  // effective atomic weight
	dREAL zEff;  // effective atomic number

	const isotope* fIsotope;

	dREAL pnBlackTrackEnergy;  // the kinetic energy available for
								  // proton/neutron black track particles
	dREAL dtaBlackTrackEnergy; // the kinetic energy available for
								  // deuteron/triton/alpha particles
	dREAL pnBlackTrackEnergyfromAnnihilation;
					 // kinetic energy available for proton/neutron black
					 // track particles based on baryon annihilation
	dREAL dtaBlackTrackEnergyfromAnnihilation;
					 // kinetic energy available for deuteron/triton/alpha
					 // black track particles based on baryon annihilation


	// ************************** member variables by ChV *******************
  // Excitation Energy leading to evaporation or deexcitation.
	dREAL  excitationEnergy;

  // Momentum, accumulated by absorbing Particles
	dREAL momentum;

  // Fermi Gas model: at present, we assume constant nucleon density for all
  // nuclei. The radius of a nucleon is taken to be 1 fm.
  // see for example S.Fl"ugge, Encyclopedia of Physics, Vol XXXIX,
  // Structure of Atomic Nuclei (Berlin-Gottingen-Heidelberg, 1957) page 426.

  // maximum momentum possible from fermi gas model:
	dREAL fermiMomentum;
	dREAL theTemp; // temperature
	// ****************************** end ChV ******************************

} nucleus;

extern void constructNucleus(nucleus* nuc);
extern void constructNucleus_AZ(nucleus* nuc, const dREAL A, const dREAL Z);
//extern void constructNucleus(nucleus nuc, const int A, const int Z);
extern void constructNucleus_Material(nucleus* nuc, const material aMaterial);

//~G4Nucleus();
//inline G4Nucleus( const G4Nucleus &right )
//{ *this = right; }
//inline G4Nucleus& operator = (const G4Nucleus& right){
//  if (this != &right) {
//	theA=right.theA;
//	theZ=right.theZ;
//	aEff=right.aEff;
//	zEff=right.zEff;
//	fIsotope = right.fIsotope;
//	pnBlackTrackEnergy=right.pnBlackTrackEnergy;
//	dtaBlackTrackEnergy=right.dtaBlackTrackEnergy;
//	pnBlackTrackEnergyfromAnnihilation =
//				 right.pnBlackTrackEnergyfromAnnihilation;
//	dtaBlackTrackEnergyfromAnnihilation =
//				 right.dtaBlackTrackEnergyfromAnnihilation;
//	theTemp = right.theTemp;
//	excitationEnergy = right.excitationEnergy;
//	momentum = right.momentum;
//	fermiMomentum = right.fermiMomentum;
//  }
//  return *this;
//}
//inline G4bool operator==( const G4Nucleus &right ) const
//{ return ( this == (G4Nucleus *) &right ); }
//inline G4bool operator!=( const G4Nucleus &right ) const
//{ return ( this != (G4Nucleus *) &right ); }

extern void nucleus_ChooseParameters(nucleus* nuc, const material aMaterial );
extern void nucleus_SetParameters(nucleus* nuc, const dREAL A, const dREAL Z );
//void SetParameters( const G4int A, const G4int Z );

/*
#ifndef G4Hadr_Nucleus_IntegerAZ
//deprecated Jan 2010, GF
inline G4double GetN() const
{ return aEff; }

inline G4double GetZ() const
{ return zEff; }
#endif
//to be replaced by new
*/

extern int nucleus_GetA_asInt(nucleus* nuc);
extern int nucleus_GetN_asInt(nucleus* nuc);
extern int nucleus_GetZ_asInt(nucleus* nuc);
extern const isotope* nucleus_GetIsotope(nucleus* nuc);
extern void nucleus_SetIsotope(nucleus* nuc, const isotope* iso);

extern dREAL nucleus_GetPNBlackTrackEnergy(nucleus* nuc);
extern dREAL nucleus_GetDTABlackTrackEnergy(nucleus* nuc);
extern dREAL nucleus_GetAnnihilationPNBlackTrackEnergy(nucleus* nuc);
extern dREAL nucleus_GetAnnihilationDTABlackTrackEnergy(nucleus* nuc);

extern dREAL nucleus_Cinema(nucleus* nuc, dREAL kineticEnergy );
extern dREAL nucleus_EvaporationEffects(nucleus* nuc, dREAL kineticEnergy );
extern dREAL nucleus_AnnihilationEvaporationEffects(nucleus* nuc, dREAL kineticEnergy, dREAL ekOrg);

extern dREAL nucleus_AtomicMass(nucleus* nuc, const dREAL A, const dREAL Z );
extern dREAL nucleus_AtomicMass_int(nucleus* nuc, const int A, const int Z );
extern dREAL nucleus_GetThermalPz(nucleus* nuc, const dREAL mass, const dREAL temp );

extern dynamicParticle *nucleus_ReturnTargetParticle(nucleus* nuc);

//FIXME: G4DynamicParticle, G4ReactionProduct, G4ThreeVector

extern reactionProduct nucleus_GetThermalNucleus(nucleus* nuc, dREAL aMass, dREAL temp/*=-1*/);
extern reactionProduct nucleus_GetBiasedThermalNucleus(nucleus* nuc, dREAL aMass, threeVector aVelocity, dREAL temp/*=-1*/);

// ******************  methods introduced by ChV ***********************
// return excitation Energy
extern dREAL nucleus_GetEnergyDeposit(nucleus* nuc);
// excitation Energy...
extern void nucleus_AddExcitationEnergy(nucleus* nuc, dREAL anEnergy);


//FIXME: G4ReactionProductVector, G4ThreeVector

// return fermi momentum
extern threeVector nucleus_GetFermiMomentum(nucleus* nuc);
//  final nucleus fragmentation. Return List of particles
// which should be used for further tracking.
//extern G4ReactionProductVector* nucleus_Fragmentate(nucleus* nuc);
// momentum of absorbed Particles ..
extern void nucleus_AddMomentum(nucleus* nuc, const threeVector aMomentum);

/*
  // return particle to be absorbed.
     G4DynamicParticle* ReturnAbsorbingParticle(G4double weight);
*/

// ****************************** end ChV ******************************
