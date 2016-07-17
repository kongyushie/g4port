//source/processes/hadronic/util/include/G4HadFinalState.h

#include <string.h>
#include "debug.h"
#include "type.h"
#include "g4/element.h"
#include "g4/dynamicParticle.h"
#include "g4/threeVector.h"
#include "g4/lorentzRotation.h"

typedef enum G4HadFinalStateStatus{isAlive, stopAndKill, suspend} G4HadFinalStateStatus;

typedef struct hadFinalState{
	threeVector theDirection;
	dREAL theEnergy;
//	std::vector<G4HadSecondary> theSecs;
	G4HadFinalStateStatus theStat;
	lorentzRotation theT;
	dREAL theW;
	dREAL theEDep;
} hadFinalState;


extern void constructHadFinalState(hadFinalState* had);

//extern int hadFinalState_GetNumberOfSecondaries(hadFinalState* had);
extern void hadFinalState_SetEnergyChange(hadFinalState* had, dREAL anEnergy);
extern dREAL hadFinalState_GetEnergyChange(hadFinalState* had);
extern void hadFinalState_SetMomentumChange_3V(hadFinalState* had, const threeVector/*&*/ aV);
extern void hadFinalState_SetMomentumChange(hadFinalState* had, dREAL x, dREAL y, dREAL z);
extern const threeVector/*&*/ hadFinalState_GetMomentumChange(hadFinalState* had);
//extern void hadFinalState_AddSecondary(hadFinalState* had, dynamicParticle* aP, int mod=-1);
//extern void hadFinalState_AddSecondary(hadFinalState* had, const G4HadSecondary& aHS);
extern void hadFinalState_SetStatusChange(hadFinalState* had, G4HadFinalStateStatus aS);
extern G4HadFinalStateStatus hadFinalState_GetStatusChange(hadFinalState* had);
extern void hadFinalState_Clear(hadFinalState* had);
extern const lorentzRotation/*&*/ hadFinalState_GetTrafoToLab(hadFinalState* had);
extern void hadFinalState_SetTrafoToLab(hadFinalState* had, const lorentzRotation /*&*/ aT);
extern void hadFinalState_SetWeightChange(hadFinalState* had, dREAL aW);
extern dREAL hadFinalState_GetWeightChange(hadFinalState* had);
//extern G4HadSecondary* hadFinalState_GetSecondary(hadFinalState* had, size_t i);
//extern const G4HadSecondary* hadFinalState_GetSecondary(hadFinalState* had, size_t i);
extern void hadFinalState_SetLocalEnergyDeposit(hadFinalState* had, dREAL aE);
extern dREAL hadFinalState_GetLocalEnergyDeposit(hadFinalState* had);
//  void SecondariesAreStale();    // Deprecated; not needed for values
//extern void hadFinalState_ClearSecondaries(hadFinalState* had);

//// Concatenate lists efficiently
//extern void hadFinalState_AddSecondaries(hadFinalState* had, const std::vector<G4HadSecondary>& addSecs);
//extern void hadFinalState_AddSecondaries(hadFinalState* had, const hadFinalState/*&*/ addHFS);
//extern void hadFinalState_AddSecondaries(hadFinalState* had, const hadFinalState* addHFS);
