#include <stdio.h>
#include <math.h>
#include "dynamicParticle.h"
#include <stdbool.h>
static const dREAL EnergyMomentumRelationAllowance = 1.e-3*1.;//KeV in systemUnit.h:123

//- constructors
/*
void constructDynamicParticle(dynamicParticle* dp){
			dp->theMomentumDirection(0.0,0.0,1.0);
			dp->theParticleDefinition = 0;
			dp->theKineticEnergy = 0.0;
			dp->theProperTime = 0.0;
			dp->theDynamicalMass = 0.0;
			dp->theDynamicalCharge = 0.0;
			dp->theDynamicalSpin = 0.0;
			dp->theDynamicalMagneticMoment = 0.0;
//			dp->theElectronOccupancy = 0;
//			dp->thePreAssignedDecayProducts = 0;
			dp->thePreAssignedDecayTime = -1.0;
//			dp->verboseLevel = 1;
//			dp->primaryParticle = 0;
			dp->thePDGcode = 0;
}
*/
/*
void constructDynamicParticle_KE(dynamicParticle* dp, const particleDefinition* aParticleDefinition,
				//const threeVector& aMomentumDirection,
				const threeVector* aMomentumDirection,
				dREAL aKineticEnergy){
			dp->theMomentumDirection = aMomentumDirection;
			dp->theParticleDefinition = aParticleDefinition;
			dp->theKineticEnergy = aKineticEnergy;
			dp->theProperTime = 0.0;
			dp->theDynamicalMass(particleDefinition_GetPDGMass(aParticleDefinition)),
			dp->theDynamicalCharge(particleDefinition_GetPDGCharge(aParticleDefinition)),
			dp->theDynamicalSpin(particleDefinition_GetPDGSpin(aParticleDefinition)),
			dp->theDynamicalMagneticMoment(particleDefinition_GetPDGMagneticMoment(aParticleDefinition)),
//			dp->theElectronOccupancy = 0;
//			dp->thePreAssignedDecayProducts = 0;
			dp->thePreAssignedDecayTime = -1.0;
//			dp->verboseLevel = 1;
//			dp->primaryParticle = 0;
			dp->thePDGcode = 0;
}
*/
/*
void constructDynamicParticle_3VPM(dynamicParticle* dp, const particleDefinition * aParticleDefinition,
				//const threeVector& aParticleMomentum){
				const threeVector* aParticleMomentum){
			dp->theMomentumDirection = aParticleDefinition;
			dp->theKineticEnergy = 0.0;
			dp->theProperTime = 0.0;
			dp->theDynamicalMass = particleDefinition_GetPDGMass(aParticleDefinition);
			dp->theDynamicalCharge = particleDefinition_GetPDGCharge(aParticleDefinition);
			dp->theDynamicalSpin = particleDefinition_GetPDGSpin(aParticleDefinition);
			dp->theDynamicalMagneticMoment = particleDefinition_GetPDGMagneticMoment(aParticleDefinition);
//			dp->theElectronOccupancy = 0;
//			dp->thePreAssignedDecayProducts = 0;
			dp->thePreAssignedDecayTime = -1.0;
//			dp->verboseLevel = 1;
//			dp->primaryParticle = 0;
			dp->thePDGcode = 0;
			dynamicParticle_SetMomentum(dp, aParticleMomentum);  // 3-dim momentum is given
}
*/
/*
void constructDynamicParticle_LVPM(dynamicParticle* dp, const particleDefinition * aParticleDefinition,
				//const lorentzVector& aParticleMomentum){
				const lorentzVector* aParticleMomentum){
			dp->theMomentumDirection = aParticleDefinition;
			dp->theKineticEnergy = 0.0;
			dp->theProperTime = 0.0;
			dp->theDynamicalMass = particleDefinition_GetPDGMass(aParticleDefinition);
			dp->theDynamicalCharge = particleDefinition_GetPDGCharge(aParticleDefinition);
			dp->theDynamicalSpin = particleDefinition_GetPDGSpin(aParticleDefinition);
			dp->theDynamicalMagneticMoment = particleDefinition_GetPDGMagneticMoment(aParticleDefinition);
//			dp->theElectronOccupancy = 0;
//			dp->thePreAssignedDecayProducts = 0;
		    dp->thePreAssignedDecayTime = -1.0;
//			dp->verboseLevel = 1;
//			dp->primaryParticle = 0;
	        dp->thePDGcode = 0;
			dynamicParticle_Set4Momentum(dp, aParticleMomentum);  // 4-momentum vector (Lorentz vector) is given
}
*/
//FIXME: threeVector
/*
void constructDynamicParticle_TE_3VPM(dynamicParticle* dp, const particleDefinition * aParticleDefinition,
	 //dREAL aTotalEnergy, const threeVector& aParticleMomentum){
	 dREAL aTotalEnergy, const threeVector* aParticleMomentum){
			dp->theMomentumDirection = aParticleDefinition;
			dp->theKineticEnergy = 0.0;
			dp->theProperTime = 0.0;
			dp->theDynamicalMass = particleDefinition_GetPDGMass(aParticleDefinition);
			dp->theDynamicalCharge = particleDefinition_GetPDGCharge(aParticleDefinition);
			dp->theDynamicalSpin = particleDefinition_GetPDGSpin(aParticleDefinition);
			dp->theDynamicalMagneticMoment = particleDefinition_GetPDGMagneticMoment(aParticleDefinition);
//			dp->theElectronOccupancy = 0;
//			dp->thePreAssignedDecayProducts = 0;
			dp->thePreAssignedDecayTime = -1.0;
//			dp->verboseLevel = 1;
//			dp->primaryParticle = 0;
			dp->thePDGcode = 0;

	   // total energy and 3-dim momentum are given
			dREAL pModule2 = threeVector_mag2(aParticleMomentum);
		 if (pModule2>0.0) {
			 dREAL mass2 = aTotalEnergy*aTotalEnergy - pModule2;
		   if(mass2 < EnergyMomentumRelationAllowance*EnergyMomentumRelationAllowance) {
			   dp->theDynamicalMass = 0.;
			 SetMomentumDirection(threeVector_unit(aParticleMomentum));
			 SetKineticEnergy(aTotalEnergy);
		   } else {
			   dp->theDynamicalMass = sqrt(mass2);
			   dynamicParticle_SetMomentum(dp, aParticleMomentum);
		   }
		 } else {
			 dynamicParticle_SetMomentumDirection(dp, 1.0,0.0,0.0);
			 dynamicParticle_SetKineticEnergy(dp, 0.0);
		 }
}
*/
/*
void constructDynamicParticle_KE_DM(dynamicParticle* dp, const particleDefinition * aParticleDefinition,
	   //const threeVector& aMomentumDirection,
	   const threeVector* aMomentumDirection,
	   dREAL aKineticEnergy,
	   const dREAL dynamicalMass){
			dp->theMomentumDirection = aMomentumDirection;
			dp->theParticleDefinition = aParticleDefinition;
			dp->theKineticEnergy = aKineticEnergy;
			dp->theProperTime = 0.0;
			dp->theDynamicalMass = dynamicalMass;
			dp->theDynamicalCharge = particleDefinition_GetPDGCharge(aParticleDefinition);
			dp->theDynamicalSpin = particleDefinition_GetPDGSpin(aParticleDefinition);
			dp->theDynamicalMagneticMoment = particleDefinition_GetPDGMagneticMoment(aParticleDefinition);
//			dp->theElectronOccupancy = 0;
//			dp->thePreAssignedDecayProducts = 0;
			dp->thePreAssignedDecayTime = -1.0;
//			dp->verboseLevel = 1;
//			dp->primaryParticle = 0;
			dp->thePDGcode = 0;
}
*/
//     constructDynamicParticle(dynamicParticle dp, const G4DynamicParticle &right);

//- Set/Get methods
/*
//const threeVector& dynamicParticle_GetMomentumDirection(dynamicParticle* dp){
const threeVector* dynamicParticle_GetMomentumDirection(dynamicParticle* dp){
//  Returns the normalized direction of the momentum
	return dp->theMomentumDirection;
}
*/
/*
//void dynamicParticle_SetMomentumDirection_3V(dynamicParticle* dp, const threeVector& aDirection){
void dynamicParticle_SetMomentumDirection_3V(dynamicParticle* dp, const threeVector* aDirection){
//  Sets the normalized direction of the momentum
	dp->theMomentumDirection = aDirection;
}
*/
//FIXME: threeVector
/*
void dynamicParticle_SetMomentumDirection(dynamicParticle* dp, dREAL px, dREAL py, dREAL pz){
//  Sets the normalized direction of the momentum by coordinates
	threeVector_setX(dp->theMomentumDirection, px);
	threeVector_setY(dp->theMomentumDirection, py);
	threeVector_setZ(dp->theMomentumDirection, pz);
}
*/

//FIXME: threeVector
/*
threeVector dynamicParticle_GetMomentum(dynamicParticle* dp){
//  Returns the current particle momentum vector
	dREAL pModule = sqrt(dp->theKineticEnergy*dp->theKineticEnergy +
	                        2*dp->theKineticEnergy*dp->theDynamicalMass);
	threeVector pMomentum(threeVector_x(dp->theMomentumDirection)*pModule,
							threeVector_y(dp->theMomentumDirection)*pModule,
							threeVector_z(dp->theMomentumDirection)*pModule);
	  return pMomentum;
}
*/
//FIXME: threeVector
/*
//void dynamicParticle_SetMomentum(dynamicParticle* dp, const threeVector& momentum){
void dynamicParticle_SetMomentum(dynamicParticle* dp, const threeVector* momentum){
//  set the current particle momentum vector
	dREAL pModule2 = momentum.mag2();
	  if (pModule2>0.0) {
		  dREAL mass = dp->theDynamicalMass;
		  dynamicParticle_SetMomentumDirection(dp, momentum.unit());
		  dynamicParticle_SetKineticEnergy(dp, sqrt(pModule2 + mass*mass)-mass);
	  } else {
		  dynamicParticle_SetMomentumDirection(dp, 1.0,0.0,0.0);
		  dynamicParticle_SetKineticEnergy(dp, 0.0);
	  }
}
*/
//FIXME: lorentzVector
/*
lorentzVector dynamicParticle_Get4Momentum(dynamicParticle* dp){
//  Returns the current particle energy-momentum 4vector
	dREAL mass      = dp->theDynamicalMass;
	dREAL energy    = dp->theKineticEnergy;
	dREAL momentum  = sqrt(energy*energy+2.0*mass*energy);
	lorentzVector    p4( threeVector_x(dp->theMomentumDirection)*momentum,
			  	  	  	  threeVector_y(dp->theMomentumDirection)*momentum,
			  	  	  	  threeVector_z(dp->theMomentumDirection)*momentum,
			  	  	  	  energy+mass);
	  return p4;
}
*/
//FIXME: lorentzVector
/*
//void dynamicParticle_Set4Momentum(dynamicParticle* dp, const lorentzVector& momentum){
void dynamicParticle_Set4Momentum(dynamicParticle* dp, const lorentzVector* momentum){
//  Set the current particle energy-momentum 4vector
	dREAL pModule2 = momentum.vect().mag2();
	  if (pModule2>0.0) {
		  dynamicParticle_SetMomentumDirection(dp, momentum.vect().unit());
	    dREAL totalenergy = momentum.t();
	    dREAL mass2 = totalenergy*totalenergy - pModule2;
	    if(mass2 < EnergyMomentumRelationAllowance*EnergyMomentumRelationAllowance) {
	    	dp->theDynamicalMass = 0.;
	    	dynamicParticle_SetKineticEnergy(dp, totalenergy);
	    } else {
	    	dp->theDynamicalMass = sqrt(mass2);
	    	dynamicParticle_SetKineticEnergy(dp, totalenergy-dp->theDynamicalMass);
	    }
	  } else {
		  dynamicParticle_SetMomentumDirection(dp, 1.0,0.0,0.0);
		  dynamicParticle_SetKineticEnergy(dp, 0.0);
	  }
}
*/
/*
dREAL dynamicParticle_GetTotalMomentum(dynamicParticle* dp){
//  Returns the module of the momentum vector
	// The momentum is returned in energy equivalent.
	  return sqrt((dp->theKineticEnergy + 2.*dp->theDynamicalMass)* dp->theKineticEnergy);
}
*/
/*
dREAL dynamicParticle_GetTotalEnergy(dynamicParticle* dp){
//  Returns the total energy of the particle
	return (dp->theKineticEnergy+dp->theDynamicalMass);
}
*/

dREAL dynamicParticle_GetKineticEnergy(dynamicParticle* dp){
//  Returns the kinetic energy of a particle
	return (*dp).theKineticEnergy;
}

/*
void dynamicParticle_SetKineticEnergy(dynamicParticle* dp, dREAL aEnergy){
//  Sets the kinetic energy of a particle
	dp->theKineticEnergy = aEnergy;
}
*/
/*
dREAL dynamicParticle_GetProperTime(dynamicParticle* dp){
//  Returns the current particle proper time
	return dp->theProperTime;
}
*/
/*
void dynamicParticle_SetProperTime(dynamicParticle* dp,  dREAL atime){
	//  Set the current particle Proper Time
	dp->theProperTime = atime;
}
*/
/*
//const threeVector& dynamicParticle_GetPolarization(dynamicParticle* dp){
const threeVector* dynamicParticle_GetPolarization(dynamicParticle* dp){
	return dp->thePolarization;
}
*/
//FIXME: threeVector
/*
void dynamicParticle_SetPolarization(dynamicParticle* dp, dREAL polX, dREAL polY, dREAL polZ){
//   Set/Get polarization vector
	threeVector_setX(dp->thePolarization, polX);
	threeVector_setY(dp->thePolarization, polY);
	threeVector_setZ(dp->thePolarization, polZ);
}
*/

// set/get dynamical mass
// the dynamical mass is set to PDG mass in default
/*
dREAL dynamicParticle_GetMass(dynamicParticle* dp){
	return dp->theDynamicalMass;
}
*/
/*
void dynamicParticle_SetMass(dynamicParticle* dp, dREAL mass){
	dp->theDynamicalMass = mass;
}
*/
/*
dREAL dynamicParticle_GetCharge(dynamicParticle* dp){
	return dp->theDynamicalCharge;
}
*/
// set/get dynamical charge
// the dynamical mass is set to PDG charge in default
/*
void     dynamicParticle_SetCharge(dynamicParticle* dp, dREAL charge){
	dp->theDynamicalCharge = charge;
}
*/
//FIXME: CLHEP::eplus
/*
void     dynamicParticle_SetCharge_eplus(dynamicParticle* dp, int    chargeInUnitOfEplus){
	dp->theDynamicalCharge = chargeInUnitOfEplus*CLHEP::eplus;
}
*/
// set/get dynamical spin
//uthe dynamical spin is set to PDG spin in default
/*
dREAL dynamicParticle_GetSpin(dynamicParticle* dp){
	return dp->theDynamicalSpin;
}
*/
/*
void     dynamicParticle_SetSpin(dynamicParticle* dp, dREAL spin){
	dp->theDynamicalSpin = spin;
}
*/
/*
void     dynamicParticle_SetSpin(dynamicParticle* dp, int    spinInUnitOfHalfInteger){
	dp->theDynamicalSpin =  spinInUnitOfHalfInteger * 0.5;
}
*/
// set/get dynamical MagneticMoment
// the dynamical mass is set to PDG MagneticMoment in default
/*
dREAL dynamicParticle_GetMagneticMoment(dynamicParticle* dp){
	return dp->theDynamicalMagneticMoment;
}
*/
/*
void     dynamicParticle_SetMagneticMoment(dynamicParticle* dp, dREAL magneticMoment){
	dp->theDynamicalMagneticMoment = magneticMoment;
}
*/

// Get electron occupancy
// ElectronOccupancy is valid only if the particle is ion

//const G4ElectronOccupancy* dynamicParticle_GetElectronOccupancy(dynamicParticle dp){
//	return dp.theElectronOccupancy;
//}

//int  dynamicParticle_GetTotalOccupancy(dynamicParticle dp){
//	int value = 0;
//	  if ( dp.theElectronOccupancy != 0) {
//	    value = dp.theElectronOccupancy->GetTotalOccupancy();
//	  }
//	  return value;
//}

//int  dynamicParticle_GetOccupancy(dynamicParticle dp, int orbit){
//	int value = 0;
//	  if ( dp.theElectronOccupancy != 0) {
//	    value = dp.theElectronOccupancy->GetOccupancy(orbit);
//	  }
//	  return value;
//}

//FIXME: CLHEP::eplus, G4ElectronOccupancy
void   dynamicParticle_AddElectron(dynamicParticle* dp, int orbit, int number/* = 1*/){
		BLURT;
		printf("G4ElectronOccupancy not implemented \n");
		exit(0);
//	if ( dp.theElectronOccupancy == 0) AllocateElectronOccupancy();
//	  if ( dp.theElectronOccupancy != 0) {
//	    G4int n = dp.theElectronOccupancy->AddElectron(orbit, number );
//	    dp.theDynamicalCharge -= CLHEP::eplus * n;
//	    dp.theDynamicalMass += GetElectronMass() * n;
//	  }
}

//FIXME: CLHEP::eplus, G4ElectronOccupancy
void   dynamicParticle_RemoveElectron(dynamicParticle* dp, int orbit, int number/* = 1*/){
	BLURT;
	printf("G4ElectronOccupancy not implemented \n");
	exit(0);
//	if ( dp.theElectronOccupancy == 0) AllocateElectronOccupancy();
//	 if ( dp.theElectronOccupancy != 0) {
//	    int n = theElectronOccupancy->RemoveElectron(orbit, number );
//	    dp.theDynamicalCharge += CLHEP::eplus * n;
//	    dp.theDynamicalMass -= GetElectronMass() * n;
//	  }
}


//   Set/Get particle definition
const particleDefinition* dynamicParticle_GetParticleDefinition(dynamicParticle* dp){
	return (*dp).theParticleDefinition;
}

void dynamicParticle_SetDefinition(dynamicParticle* dp, const particleDefinition * aParticleDefinition){
	BLURT;
	printf("G4ElectronOccupancy not implemented \n");
	exit(0);
//	// remove preassigned decay
//	  if (dp.thePreAssignedDecayProducts != 0) {
//	#ifdef G4VERBOSE
//	    if (verboseLevel>0) {
//	      G4cout << " G4DynamicParticle::SetDefinition()::"
//	             << "!!! Pre-assigned decay products is attached !!!! " << G4endl;
//	      G4cout << "!!! New Definition is " << aParticleDefinition->GetParticleName()
//		     << " !!! " << G4endl;
//	      G4cout << "!!! Pre-assigned decay products will be deleted !!!! " << G4endl;
//	    }
//	#endif
//	    delete thePreAssignedDecayProducts;
//	  }
//	  thePreAssignedDecayProducts = 0;
//
//	  theParticleDefinition = aParticleDefinition;
//
//	  // set Dynamic mass/chrge
//	  theDynamicalMass = theParticleDefinition->GetPDGMass();
//	  theDynamicalCharge = theParticleDefinition->GetPDGCharge();
//	  theDynamicalSpin = theParticleDefinition->GetPDGSpin();
//	  theDynamicalMagneticMoment = theParticleDefinition->GetPDGMagneticMoment();
//
//	  // Set electron orbits
//	  if (theElectronOccupancy != 0) delete theElectronOccupancy;
//	  theElectronOccupancy =0;
//	  //AllocateElectronOccupancy();
}

//  following method of GetDefinition remains
//  because of backward compatiblity. It will be removed in future
particleDefinition* dynamicParticle_GetDefinition(dynamicParticle* dp){
	return (*dp).theParticleDefinition;
}

/*
const G4DecayProducts *dynamicParticle_GetPreAssignedDecayProducts(dynamicParticle* dp){
	BLURT;
	printf("G4DecayProducts not implemented \n");
	exit(0);
}
*/
//   Set/Get pre-assigned decay channel
/*
void dynamicParticle_SetPreAssignedDecayProducts(dynamicParticle* dp, G4DecayProducts *aDecayProducts){
	BLURT;
	printf("G4DecayProducts not implemented \n");
	exit(0);
}
*/
dREAL dynamicParticle_GetPreAssignedDecayProperTime(dynamicParticle* dp){
	BLURT;
	printf("G4DecayProducts not implemented \n");
	exit(0);
}

//   Set/Get pre-assigned proper time when the particle will decay
void dynamicParticle_SetPreAssignedDecayProperTime(dynamicParticle* dp, dREAL aDecayProducts){
	BLURT;
	printf("G4DecayProducts not implemented \n");
	exit(0);
//	dp.thePreAssignedDecayProducts = aDecayProducts;
}

//void dynamicParticle_SetPrimaryParticle(dynamicParticle dp, G4PrimaryParticle* p){
//
//}

void dynamicParticle_SetPDGcode(dynamicParticle* dp, int c){
	dp->thePDGcode = c;
}

// Return the pointer to the corresponding G4PrimaryParticle object
// if this particle is a primary particle OR is defined as a pre-assigned
// decay product. Otherwise return null.

//G4PrimaryParticle* dynamicParticle_GetPrimaryParticle(dynamicParticle dp){
//	BLURT;
//	printf("G4PrimaryParticle not implemented \n");
//	exit(0);
////	return dp.primaryParticle;
//}

int dynamicParticle_GetPDGcode(dynamicParticle* dp){
	int code = particleDefinition_GetPDGEncoding(dp->theParticleDefinition);
	if(code==0) code = dp->thePDGcode;
	return code;
}

bool OP_EQ_dynamicParticle(const dynamicParticle* a, const dynamicParticle* b)
{
	return (*a).theKineticEnergy			== (*b).theKineticEnergy
		&&(*a).theProperTime				== (*b).theProperTime	
		&&(*a).theDynamicalMass			== (*b).theDynamicalMass
		&&(*a).theDynamicalCharge			== (*b).theDynamicalCharge
		&&(*a).theDynamicalSpin			== (*b).theDynamicalSpin
		&&(*a).theDynamicalMagneticMoment == (*b).theDynamicalMagneticMoment
		&&(*a).thePreAssignedDecayTime	== (*b).thePreAssignedDecayTime	
		&&(*a).thePDGcode					== (*b).thePDGcode;
}
void OP_ASSIGN_dynamicParticle(particleDefinition* a, const particleDefinition* b)
{
		(*a).theParticleName		= (*b).theParticleName;
		(*a).thePDGMass				= (*b).thePDGMass;
		(*a).thePDGWidth			= (*b).thePDGWidth;
		(*a).thePDGCharge			= (*b).thePDGCharge;
		(*a).thePDGiSpin			= (*b).thePDGiSpin;
}
