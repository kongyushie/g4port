
#include <math.h>
#include "g4/hadProjectile.h"


void constructHadProjectile(hadProjectile* had){
	had->theMat = 0;
	had->theDef = 0;
	had->theTime = 0.0;
	had->theBoundEnergy = 0.0;
}

void constructHadProjectile_track(hadProjectile* had, const G4Track /*&*/aT){
	hadProjectile_Initialise(had, aT);
}

void constructHadProjectile_dp(hadProjectile* had, const dynamicParticle /*&*/aT){
	had->theMat = 0;//NULL;
	had->theOrgMom = dynamicParticle_Get4Momentum(&aT);
	had->theDef = dynamicParticle_GetDefinition(&aT);
	lorentzRotation toZ;
	lorentzRotation_rotateZ(&toZ, -lorentzVector_phi(had->theOrgMom));
	lorentzRotation_rotateY(&toZ, -lorentzVector_theta(had->theOrgMom));

	//FIXME: operator overflow
	had->theMom = toZ*had->theOrgMom;
	had->toLabFrame = lorentzRotation_inverse(&toZ);
	had->theTime = 0.0;
	had->theBoundEnergy = 0.0;
}


//FIXME: G4Track
void hadProjectile_Initialise(hadProjectile* had, const G4Track /*&*/aT){
	had->theMat = aT.GetMaterial();
	had->theOrgMom = aT.GetDynamicParticle()->Get4Momentum();
	had->theDef = aT.GetDefinition();

	lorentzRotation toZ;
	lorentzRotation_rotateZ(&toZ, -lorentzVector_phi(had->theOrgMom));
	lorentzRotation_rotateY(&toZ, -lorentzVector_theta(had->theOrgMom));

	//FIXME: operator overflow
	had->theMom = toZ*had->theOrgMom;
	had->toLabFrame = lorentzRotation_inverse(&toZ);

	//VI time of interaction starts from zero
	//   not global time of a track
	had->theTime = 0.0;
	had->theBoundEnergy = 0.0;
}

const material * hadProjectile_GetMaterial(hadProjectile* had){
	return had->theMat;
}

const particleDefinition * hadProjectile_GetDefinition(hadProjectile* had){
	return had->theDef;
}

const lorentzVector /*&*/ hadProjectile_Get4Momentum(hadProjectile* had){
	  return had->theMom;
}

lorentzRotation /*&*/ hadProjectile_GetTrafoToLab(hadProjectile* had){
	return had->toLabFrame;
}

dREAL hadProjectile_GetKineticEnergy(hadProjectile* had){
	dREAL ekin = hadProjectile_GetTotalEnergy(had) - particleDefinition_GetPDGMass(hadProjectile_GetDefinition(had));
	if(ekin < 0.0) { ekin = 0.0; }
	return ekin;
}

dREAL hadProjectile_GetTotalEnergy(hadProjectile* had){
	return lorentzVector_e(hadProjectile_Get4Momentum(had));
}

dREAL hadProjectile_GetTotalMomentum(hadProjectile* had){
	return threeVector_mag(lorentzVector_vect(hadProjectile_Get4Momentum(had)));
}

dREAL hadProjectile_GetGlobalTime(hadProjectile* had){
	return had->theTime;
}

dREAL hadProjectile_GetBoundEnergy(hadProjectile* had){
	return had->theBoundEnergy;
}

void hadProjectile_SetGlobalTime(hadProjectile* had, dREAL t){
	had->theTime = t;
}

void hadProjectile_SetBoundEnergy(hadProjectile* had, dREAL e){
	had->theBoundEnergy = e;
}

