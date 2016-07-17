// source/particles/management/include/G3DynamicParticle.hh

// Class Description:
//  The dynamic particle is a class which contains the purely
//  dynamic aspects of a moving particle. It also has a
//  pointer to a G4ParticleDefinition object, which holds
//  all the static information.
//
#ifndef DYNAMICPARTICLE_H
#define DYNAMICPARTICLE_H 1

#include "debug.h"
#include "type.h"
#include "particleDefinition.h"
#include "threeVector.h"
#include "lorentzVector.h"
#include <stdbool.h>
//FIXME: G4ThreeVector, G4ParticleDefinition
typedef struct dynamicParticle {
	threeVector theMomentumDirection;
	//  The normalized momentum vector

	const particleDefinition *theParticleDefinition;
	//  Contains the static information of this particle.

	threeVector thePolarization;

	dREAL theKineticEnergy;
	dREAL theProperTime;
	dREAL theDynamicalMass;
	dREAL theDynamicalCharge;
	dREAL theDynamicalSpin;
	dREAL theDynamicalMagneticMoment;
//	G4ElectronOccupancy* theElectronOccupancy;
//	G4DecayProducts *thePreAssignedDecayProducts;
	dREAL thePreAssignedDecayTime;

//protected:
    int thePDGcode;

//	 protected:
//	   G4int verboseLevel;
//	 protected:
//	    G4PrimaryParticle* primaryParticle;
//	    // This void pointer is used by G4EventManager to maintain the
//	    // link between pre-assigned decay products and corresponding
//	    // primary particle.
} dynamicParticle;

//- constructors
     extern void constructDynamicParticle(dynamicParticle* dp);

     extern void constructDynamicParticle_KE(dynamicParticle* dp, const particleDefinition* aParticleDefinition,
                        const threeVector/*&*/ aMomentumDirection,
                        dREAL aKineticEnergy);

     extern void constructDynamicParticle_3VPM(dynamicParticle* dp, const particleDefinition * aParticleDefinition,
                        const threeVector/*&*/ aParticleMomentum);

     extern void constructDynamicParticle_LVPM(dynamicParticle* dp, const particleDefinition * aParticleDefinition,
                        const lorentzVector/*&*/ aParticleMomentum);

     extern void constructDynamicParticle_TE_3VPM(dynamicParticle* dp, const particleDefinition * aParticleDefinition,
    		 dREAL aTotalEnergy,
            const threeVector/*&*/ aParticleMomentum);

     extern void constructDynamicParticle_KE_DM(dynamicParticle* dp, const particleDefinition * aParticleDefinition,
		       const threeVector/*&*/ aMomentumDirection,
		       dREAL aKineticEnergy,
		       const dREAL dynamicalMass);

//     constructDynamicParticle(dynamicParticle dp, const G4DynamicParticle &right);

 //- Set/Get methods

     extern const threeVector/*&*/ dynamicParticle_GetMomentumDirection(dynamicParticle* dp);
	   //  Returns the normalized direction of the momentum
	  extern void dynamicParticle_SetMomentumDirection_3V(dynamicParticle* dp, const threeVector/*&*/ aDirection);
	   //  Sets the normalized direction of the momentum
	  extern void dynamicParticle_SetMomentumDirection(dynamicParticle* dp, dREAL px, dREAL py, dREAL pz);
	   //  Sets the normalized direction of the momentum by coordinates

	  extern threeVector dynamicParticle_GetMomentum(dynamicParticle* dp);
	   //  Returns the current particle momentum vector
	  extern void dynamicParticle_SetMomentum(dynamicParticle* dp, const threeVector/*&*/ momentum);
	   //  set the current particle momentum vector

	  extern lorentzVector dynamicParticle_Get4Momentum(dynamicParticle* dp);
	   //  Returns the current particle energy-momentum 4vector
	  extern void dynamicParticle_Set4Momentum(dynamicParticle* dp, const lorentzVector/*&*/ momentum);
	   //  Set the current particle energy-momentum 4vector


	  extern dREAL dynamicParticle_GetTotalMomentum(dynamicParticle* dp);
	   //  Returns the module of the momentum vector
	  extern dREAL dynamicParticle_GetTotalEnergy(dynamicParticle* dp);
	   //  Returns the total energy of the particle

	  extern dREAL dynamicParticle_GetKineticEnergy(dynamicParticle* dp);
	   //  Returns the kinetic energy of a particle
	  extern void dynamicParticle_SetKineticEnergy(dynamicParticle* dp, dREAL aEnergy);
	   //  Sets the kinetic energy of a particle


	  extern dREAL dynamicParticle_GetProperTime(dynamicParticle* dp);
	   //  Returns the current particle proper time
	  extern void dynamicParticle_SetProperTime(dynamicParticle* dp,  dREAL );
	   //  Set the current particle Proper Time


	  extern const threeVector/*&*/ dynamicParticle_GetPolarization(dynamicParticle* dp);
	  extern void dynamicParticle_SetPolarization(dynamicParticle* dp, dREAL polX, dREAL polY, dREAL polZ);
	   //   Set/Get polarization vector


	  extern dREAL dynamicParticle_GetMass(dynamicParticle* dp);
	  extern void     dynamicParticle_SetMass(dynamicParticle* dp, dREAL mass);
	  // set/get dynamical mass
	  // the dynamical mass is set to PDG mass in default


	  extern dREAL dynamicParticle_GetCharge(dynamicParticle* dp);
	  extern void     dynamicParticle_SetCharge(dynamicParticle* dp, dREAL charge);
	  extern void     dynamicParticle_SetCharge_eplus(dynamicParticle* dp, int    chargeInUnitOfEplus);
	  // set/get dynamical charge
	  // the dynamical mass is set to PDG charge in default

	  extern dREAL dynamicParticle_GetSpin(dynamicParticle* dp);
	  extern void dynamicParticle_SetSpin(dynamicParticle* dp, dREAL spin);
	  extern void dynamicParticle_SetSpin_int(dynamicParticle* dp, int spinInUnitOfHalfInteger);
	  // set/get dynamical spin
	  // the dynamical spin is set to PDG spin in default

	  extern dREAL dynamicParticle_GetMagneticMoment(dynamicParticle* dp);
	  extern void     dynamicParticle_SetMagneticMoment(dynamicParticle* dp, dREAL magneticMoment);
	  // set/get dynamical MagneticMoment
	  // the dynamical mass is set to PDG MagneticMoment in default


	  //extern const G4ElectronOccupancy* dynamicParticle_GetElectronOccupancy(dynamicParticle* dp);
	  // Get electron occupancy
	  // ElectronOccupancy is valid only if the particle is ion
	  extern int  dynamicParticle_GetTotalOccupancy(dynamicParticle* dp);
	  extern int  dynamicParticle_GetOccupancy(dynamicParticle* dp, int orbit);
	  extern void   dynamicParticle_AddElectron(dynamicParticle* dp, int orbit, int number/* = 1*/);
	  extern void   dynamicParticle_RemoveElectron(dynamicParticle* dp, int orbit, int number/* = 1*/);


	  extern const particleDefinition* dynamicParticle_GetParticleDefinition(dynamicParticle* dp);
	  extern void dynamicParticle_SetDefinition(dynamicParticle* dp, const particleDefinition * aParticleDefinition);
	  //   Set/Get particle definition
	  //  following method of GetDefinition remains
	  //  because of backward compatiblity. It will be removed in future
	  extern particleDefinition* dynamicParticle_GetDefinition(dynamicParticle* dp);


	  //extern const G4DecayProducts *dynamicParticle_GetPreAssignedDecayProducts(dynamicParticle* dp);
	 // extern void dynamicParticle_SetPreAssignedDecayProducts(dynamicParticle* dp, G4DecayProducts *aDecayProducts);
	   //   Set/Get pre-assigned decay channel

	  extern dREAL dynamicParticle_GetPreAssignedDecayProperTime(dynamicParticle* dp);
	  extern void dynamicParticle_SetPreAssignedDecayProperTime(dynamicParticle* dp, dREAL);
	   //   Set/Get pre-assigned proper time when the particle will decay

	  //extern void dynamicParticle_SetPrimaryParticle(dynamicParticle* dp, G4PrimaryParticle* p);
	  extern void dynamicParticle_SetPDGcode(dynamicParticle* dp, int c);

	  //extern G4PrimaryParticle* dynamicParticle_GetPrimaryParticle(dynamicParticle* dp);
		// Return the pointer to the corresponding G4PrimaryParticle object
		// if this particle is a primary particle OR is defined as a pre-assigned
		// decay product. Otherwise return null.

		extern int dynamicParticle_GetPDGcode(dynamicParticle* dp);
		// Return the PDG code of this particle. If the particle is known to Geant4
		// its PDG code defined in G4ParticleDefinition is returned. If it is unknown
		// (i.e. PDG code in G4ParticleDefinition is 0), PDG code defined in the
		// corresponding primary particle or pre-assigned decay product will be
		// returned if available. Otherwise (e.g. for geantino) returns 0.
extern bool OP_EQ_dynamicParticle(const dynamicParticle* a, const dynamicParticle* b);
extern void OP_ASSIGN_dynamicParticle(particleDefinition* a, const particleDefinition* b);
#endif
