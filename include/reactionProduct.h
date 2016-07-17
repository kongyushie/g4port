// source/processes/hadronic/util/include/G4ReactionProduct.hh
#pragma once
#include "particleDefinition.h"
#include "threeVector.h"
#include "debug.h"
#include "type.h"

typedef struct reactionProduct {
	const particleDefinition *theParticleDefinition;

	// for use with string models and cascade.
	threeVector positionInNucleus;
	dREAL formationTime;
	char /*G4bool*/ hasInitialStateParton;

	// mass is included here, since pseudo-particles are created with masses different
	// than the standard particle masses, and we are not allowed to create particles
	dREAL mass;

	threeVector momentum;

	dREAL totalEnergy;
	dREAL kineticEnergy;

	dREAL timeOfFlight;

	//  side refers to how the particles are distributed in the
	//  forward (+) and backward (-) hemispheres in the center of mass system
	int side;

	int theCreatorModel;

	// NewlyAdded refers to particles added by "nuclear excitation", or as
	//  "black track" particles, or as deuterons, tritons, and alphas
	char /*G4bool*/ NewlyAdded;
	char /*G4bool*/ MayBeKilled;
} reactionProduct;

extern void constructReactionProduct(reactionProduct* rp);

extern void constructReactionProduct_pd(reactionProduct* rp, const particleDefinition *aParticleDefinition );

//    ~G4ReactionProduct() {}

//    G4ReactionProduct( const G4ReactionProduct &right );

//    // Override new and delete for use with G4Allocator
//    inline void* operator new(size_t) {
//      if (!aRPAllocator) aRPAllocator = new G4Allocator<G4ReactionProduct>  ;
//      return (void *)aRPAllocator->MallocSingle();
//    }
//#ifdef __IBMCPP__
//    inline void* operator new(size_t, void *p) {
//      return p;
//    }
//#endif
//    inline void operator delete(void* aReactionProduct) {
//      aRPAllocator->FreeSingle((G4ReactionProduct*)aReactionProduct);
//    }
//
//    G4ReactionProduct &operator= ( const G4ReactionProduct &right );
//
//    G4ReactionProduct &operator= ( const G4DynamicParticle &right );
//
//    G4ReactionProduct &operator= ( const G4HadProjectile &right );
//
//    inline G4bool operator== ( const G4ReactionProduct &right ) const
//    { return ( this == (G4ReactionProduct*) &right ); }
//
//    inline G4bool operator!= ( const G4ReactionProduct &right ) const
//    { return ( this != (G4ReactionProduct*) &right ); }

extern const particleDefinition* reactionProduct_GetDefinition(reactionProduct* rp);

extern void reactionProduct_SetDefinition(reactionProduct* rp, const particleDefinition* aParticleDefinition );

extern void reactionProduct_SetDefinitionAndUpdateE(reactionProduct* rp, const particleDefinition* aParticleDefinition );

extern void reactionProduct_SetMomentum_3d(reactionProduct* rp, const dREAL x, const dREAL y, const dREAL z );

extern void reactionProduct_SetMomentum_2d(reactionProduct* rp, const dREAL x, const dREAL y );

extern void reactionProduct_SetMomentum_1d(reactionProduct* rp, const dREAL z );

extern void reactionProduct_SetMomentum_3V(reactionProduct* rp, const threeVector /*&*/mom );

extern threeVector reactionProduct_GetMomentum(reactionProduct* rp);

extern dREAL reactionProduct_GetTotalMomentum(reactionProduct* rp);

extern dREAL reactionProduct_GetTotalEnergy(reactionProduct* rp);

extern void reactionProduct_SetKineticEnergy(reactionProduct* rp, const dREAL en );

extern dREAL reactionProduct_GetKineticEnergy(reactionProduct* rp);

extern void reactionProduct_SetTotalEnergy(reactionProduct* rp, const dREAL en );

extern void reactionProduct_SetMass(reactionProduct* rp, const dREAL mas );

extern dREAL reactionProduct_GetMass(reactionProduct* rp);

extern void reactionProduct_SetTOF(reactionProduct* rp, const dREAL t );

extern dREAL reactionProduct_GetTOF(reactionProduct* rp);

extern void reactionProduct_SetSide(reactionProduct* rp, const int sid );

extern int reactionProduct_GetSide(reactionProduct* rp);

extern void reactionProduct_SetCreatorModel(reactionProduct* rp, const int mod );

extern int reactionProduct_GetCreatorModel(reactionProduct* rp);

extern void reactionProduct_SetNewlyAdded(reactionProduct* rp, const char /*G4bool*/ f );

extern char /*G4bool*/ reactionProduct_GetNewlyAdded(reactionProduct* rp);

extern void reactionProduct_SetMayBeKilled(reactionProduct* rp, const char /*G4bool*/ f );

extern char /*G4bool*/ reactionProduct_GetMayBeKilled(reactionProduct* rp);

extern void reactionProduct_SetZero(reactionProduct* rp);

extern void reactionProduct_Lorentz(reactionProduct* rp, const reactionProduct /*&*/p1, const reactionProduct /*&*/p2 );

extern dREAL reactionProduct_Angle(reactionProduct* rp, const reactionProduct /*&*/p );

extern void reactionProduct_SetPositionInNucleus(reactionProduct* rp, dREAL x, dREAL y, dREAL z);
extern void reactionProduct_SetPositionInNucleus_3V(reactionProduct* rp, threeVector /*&*/ aPosition );

extern threeVector reactionProduct_GetPositionInNucleus(reactionProduct* rp);
extern dREAL reactionProduct_GetXPositionInNucleus(reactionProduct* rp);
extern dREAL reactionProduct_GetYPositionInNucleus(reactionProduct* rp);
extern dREAL reactionProduct_GetZPositionInNucleus(reactionProduct* rp);
extern void reactionProduct_SetFormationTime(reactionProduct* rp, dREAL aTime);
extern dREAL reactionProduct_GetFormationTime(reactionProduct* rp);
extern void reactionProduct_HasInitialStateParton_char(reactionProduct* rp, char /*G4bool*/ aFlag);
extern char /*G4bool*/ reactionProduct_HasInitialStateParton(reactionProduct* rp);

