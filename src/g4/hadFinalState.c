
#include <math.h>
#include "g4/hadFinalState.h"


void constructHadFinalState(hadFinalState* had){
	had->theDirection(0,0,1);
	had->theEnergy = -1;
	had->theStat = isAlive;
	had->theW = 1.;
	had->theEDep = 0.;
}

//int hadFinalState_GetNumberOfSecondaries(hadFinalState* had) { return had->theSecs.size(); }
void hadFinalState_SetEnergyChange(hadFinalState* had, dREAL anEnergy){
	had->theEnergy=anEnergy;
	  if(had->theEnergy<0){
		  BLURT;
		  printf("Final state energy was: E =  %f\n", had->theEnergy);
		  exit(0);
	    }
}

dREAL hadFinalState_GetEnergyChange(hadFinalState* had) { return had->theEnergy; }

void hadFinalState_SetMomentumChange_3V(hadFinalState* had, const threeVector/*&*/ aV){ had->theDirection=aV; }

void hadFinalState_SetMomentumChange(hadFinalState* had, dREAL x, dREAL y, dREAL z){
	threeVector_set(had->theDirection, x,y,z);
	if(fabs(x*x + y*y + z*z - 1.0)>0.001) {
		BLURT;
		printf("We have negative threeVector_mag(had->theDirection) = %f\n",
			  threeVector_mag(had->theDirection));
		exit(0);
	}
}

const threeVector/*&*/ hadFinalState_GetMomentumChange(hadFinalState* had)  { return had->theDirection; }
//void hadFinalState_AddSecondary(hadFinalState* had, dynamicParticle* aP, int mod=-1) {
//theSecs.push_back(G4HadSecondary(aP, theW, mod));
//};
//void hadFinalState_AddSecondary(hadFinalState* had, const G4HadSecondary& aHS)   { theSecs.push_back(aHS); }
void hadFinalState_SetStatusChange(hadFinalState* had, G4HadFinalStateStatus aS) { had->theStat=aS; }
G4HadFinalStateStatus hadFinalState_GetStatusChange(hadFinalState* had) { return had->theStat; }

void hadFinalState_Clear(hadFinalState* had){
	threeVector_set(had->theDirection,0,0,1);
	had->theEnergy = -1;
	had->theStat = isAlive;
	had->theW = 1.;
	had->theEDep = 0.;
//	hadFinalState_ClearSecondaries();
}

const lorentzRotation/*&*/ hadFinalState_GetTrafoToLab(hadFinalState* had) { return had->theT; }
void hadFinalState_SetTrafoToLab(hadFinalState* had, const lorentzRotation /*&*/ aT) { had->theT = aT; }
void hadFinalState_SetWeightChange(hadFinalState* had, dREAL aW) { had->theW=aW; }
dREAL hadFinalState_GetWeightChange(hadFinalState* had) { return had->theW; }
//G4HadSecondary* hadFinalState_GetSecondary(hadFinalState* had, size_t i);
//const G4HadSecondary* hadFinalState_GetSecondary(hadFinalState* had, size_t i) ;
void hadFinalState_SetLocalEnergyDeposit(hadFinalState* had, dREAL aE) { had->theEDep=aE; }
dREAL hadFinalState_GetLocalEnergyDeposit(hadFinalState* had) { return had->theEDep;}
//  void SecondariesAreStale();    // Deprecated; not needed for values
//void hadFinalState_ClearSecondaries(hadFinalState* had)  { theSecs.clear(); }

//// Concatenate lists efficiently
//void hadFinalState_AddSecondaries(hadFinalState* had, const std::vector<G4HadSecondary>& addSecs);
//void hadFinalState_AddSecondaries(hadFinalState* had, const hadFinalState/*&*/ addHFS) {
// hadFinalState_AddSecondaries(addHFS.theSecs);
//}
//void hadFinalState_AddSecondaries(hadFinalState* had, const hadFinalState* addHFS) {
//if (addHFS) AddSecondaries(addHFS->theSecs);
//}


