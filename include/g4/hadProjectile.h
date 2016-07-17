//source/processes/hadronic/util/include/G4HadProjectile.hh

#include <string.h>
#include "debug.h"
#include "type.h"
#include "material.h"
#include "lorentzVector.h"
#include "lorentzRotation.h"
#include "particleDefinition.h"

typedef struct hadProjectile{
	const material * theMat;
	  lorentzVector theOrgMom;
	  lorentzVector theMom;
	  const particleDefinition * theDef;
	  lorentzRotation toLabFrame;
	  dREAL theTime;
	  dREAL theBoundEnergy;
} hadProjectile;

  extern void constructHadProjectile(hadProjectile* had);
  extern void constructHadProjectile_track(hadProjectile* had, const G4Track /*&*/aT);
  extern void constructHadProjectile_dp(hadProjectile* had, const dynamicParticle /*&*/aT);
//  ~G4HadProjectile();

  extern void hadProjectile_Initialise(hadProjectile* had, const G4Track /*&*/aT);

  extern const material * hadProjectile_GetMaterial(hadProjectile* had);
  extern const particleDefinition * hadProjectile_GetDefinition(hadProjectile* had);
  extern const lorentzVector /*&*/ hadProjectile_Get4Momentum(hadProjectile* had);
  extern lorentzRotation /*&*/ hadProjectile_GetTrafoToLab(hadProjectile* had);
  extern dREAL hadProjectile_GetKineticEnergy(hadProjectile* had);
  extern dREAL hadProjectile_GetTotalEnergy(hadProjectile* had);
  extern dREAL hadProjectile_GetTotalMomentum(hadProjectile* had);
  extern dREAL hadProjectile_GetGlobalTime(hadProjectile* had);
  extern dREAL hadProjectile_GetBoundEnergy(hadProjectile* had);
  extern void hadProjectile_SetGlobalTime(hadProjectile* had, dREAL t);
  extern void hadProjectile_SetBoundEnergy(hadProjectile* had, dREAL e);
