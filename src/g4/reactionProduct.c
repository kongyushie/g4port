#include "g4/reactionProduct.h"


void constructReactionProduct(reactionProduct* rp){
	rp->theParticleDefinition = 0; //(NULL)
	rp->formationTime = 0.0;
	rp->hasInitialStateParton = 0; //(false),
	rp->mass = 0.0;
	rp->totalEnergy = 0.0;
	rp->kineticEnergy = 0.0;
	rp->timeOfFlight = 0.0;
	rp->side = 0;
	rp->theCreatorModel = -1;
	rp->NewlyAdded = 0; //(false),
	rp->MayBeKilled = 1;//(true)
	reactionProduct_SetMomentum(rp, 0.0, 0.0, 0.0 );
	reactionProduct_SetPositionInNucleus(rp, 0.0, 0.0, 0.0 );
}

void constructReactionProduct_pd(reactionProduct* rp, const particleDefinition *aParticleDefinition ){
	reactionProduct_SetMomentum(rp, 0.0, 0.0, 0.0 );
	reactionProduct_SetPositionInNucleus(rp, 0.0, 0.0, 0.0 );
	    rp->formationTime = 0.0;
	    rp->hasInitialStateParton = 0;//false;
	    rp->theParticleDefinition = aParticleDefinition;
	    rp->mass = particleDefinition_GetPDGMass(aParticleDefinition);//aParticleDefinition->GetPDGMass();
	    rp->totalEnergy = rp->mass;
	    rp->kineticEnergy = 0.0;
	    (particleDefinition_GetPDGEncoding(aParticleDefinition)<0) ? rp->timeOfFlight=-1.0 : rp->timeOfFlight=1.0;
	    rp->side = 0;
	    rp->theCreatorModel = -1;
	    rp->NewlyAdded = 0;//false;
	    rp->MayBeKilled = 1;//true;
}

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

const particleDefinition* reactionProduct_GetDefinition(reactionProduct* rp)
{ return rp->theParticleDefinition; }

void reactionProduct_SetDefinition(reactionProduct* rp, const particleDefinition* aParticleDefinition ){
	rp->theParticleDefinition = aParticleDefinition;
	rp->mass = particleDefinition_GetPDGMass(aParticleDefinition);
	rp->totalEnergy = rp->mass;
	rp->kineticEnergy = 0.0;
	    (particleDefinition_GetPDGEncoding(aParticleDefinition)<0) ?
	    		rp->timeOfFlight=-1.0 : rp->timeOfFlight=1.0;
}

void reactionProduct_SetDefinitionAndUpdateE(reactionProduct* rp, const particleDefinition* aParticleDefinition ){
	dREAL aKineticEnergy = reactionProduct_GetKineticEnergy(rp);
	dREAL pp = threeVector_mag(reactionProduct_GetMomentum(rp));
	threeVector aMomentum = reactionProduct_GetMomentum(rp);
	reactionProduct_SetDefinition(rp, aParticleDefinition );
	reactionProduct_SetKineticEnergy(rp, aKineticEnergy );
	if( pp > DBL_MIN )
		reactionProduct_SetMomentum(rp, aMomentum * (sqrt(aKineticEnergy*aKineticEnergy +
									2*aKineticEnergy*GetMass())/pp) );
}

void reactionProduct_SetMomentum_3d(reactionProduct* rp, const dREAL x, const dREAL y, const dREAL z ){
	threeVector_setX( rp->momentum, x );
	threeVector_setY( rp->momentum, y );
	threeVector_setZ( rp->momentum, z );
}

void reactionProduct_SetMomentum_2d(reactionProduct* rp, const dREAL x, const dREAL y ){
	threeVector_setX( rp->momentum, x );
	threeVector_setY( rp->momentum, y );
}

void reactionProduct_SetMomentum_1d(reactionProduct* rp, const dREAL z ){
	threeVector_setZ( rp->momentum, z );
}

void reactionProduct_SetMomentum_3V(reactionProduct* rp, const threeVector /*&*/mom )
{ rp->momentum = mom; }

threeVector reactionProduct_GetMomentum(reactionProduct* rp)
{ return rp->momentum; }

dREAL reactionProduct_GetTotalMomentum(reactionProduct* rp)
{ return sqrt(abs(rp->kineticEnergy*(rp->totalEnergy+rp->mass))); }

dREAL reactionProduct_GetTotalEnergy(reactionProduct* rp)
{ return rp->totalEnergy; }

void reactionProduct_SetKineticEnergy(reactionProduct* rp, const dREAL en )
{
	rp->kineticEnergy = en;
	rp->totalEnergy = rp->kineticEnergy + rp->mass;
}

dREAL reactionProduct_GetKineticEnergy(reactionProduct* rp)
{ return rp->kineticEnergy; }

void reactionProduct_SetTotalEnergy(reactionProduct* rp, const dREAL en )
{
	rp->totalEnergy = en;
	rp->kineticEnergy = rp->totalEnergy - rp->mass;
}

void reactionProduct_SetMass(reactionProduct* rp, const dREAL mas )
{ rp->mass = mas; }

dREAL reactionProduct_GetMass(reactionProduct* rp)
{ return rp->ass; }

void reactionProduct_SetTOF(reactionProduct* rp, const dREAL t )
{ rp->timeOfFlight = t; }

dREAL reactionProduct_GetTOF(reactionProduct* rp)
{ return rp->timeOfFlight; }

void reactionProduct_SetSide(reactionProduct* rp, const int sid )
{ rp->side = sid; }

int reactionProduct_GetSide(reactionProduct* rp)
{ return rp->side; }

void reactionProduct_SetCreatorModel(reactionProduct* rp, const int mod )
{ rp->theCreatorModel = mod; }

int reactionProduct_GetCreatorModel(reactionProduct* rp)
{ return rp->theCreatorModel; }

void reactionProduct_SetNewlyAdded(reactionProduct* rp, const char /*G4bool*/ f )
{ rp->NewlyAdded = f; }

char /*G4bool*/ reactionProduct_GetNewlyAdded(reactionProduct* rp) { return rp->NewlyAdded; }

void reactionProduct_SetMayBeKilled(reactionProduct* rp, const char /*G4bool*/ f )
{ rp->MayBeKilled = f; }

char /*G4bool*/ reactionProduct_GetMayBeKilled(reactionProduct* rp)
{ return rp->MayBeKilled; }

void reactionProduct_SetZero(reactionProduct* rp){

}

void reactionProduct_Lorentz(reactionProduct* rp, const reactionProduct /*&*/p1, const reactionProduct /*&*/p2 ){

}

dREAL reactionProduct_Angle(reactionProduct* rp, const reactionProduct /*&*/p ){

}

void reactionProduct_SetPositionInNucleus(reactionProduct* rp, dREAL x, dREAL y, dREAL z)
 {
	threeVector_setX( rp->positionInNucleus, x );
	threeVector_setY( rp->positionInNucleus, y );
	threeVector_setZ( rp->positionInNucleus, z );
 }

void reactionProduct_SetPositionInNucleus_3V(reactionProduct* rp, threeVector /*&*/ aPosition )
 {
	rp->positionInNucleus = aPosition;
 }

threeVector reactionProduct_GetPositionInNucleus(reactionProduct* rp) { return rp->positionInNucleus; }
dREAL reactionProduct_GetXPositionInNucleus(reactionProduct* rp) { return rp->positionInNucleus.x(); }
dREAL reactionProduct_GetYPositionInNucleus(reactionProduct* rp) { return rp->positionInNucleus.y(); }
dREAL reactionProduct_GetZPositionInNucleus(reactionProduct* rp) { return rp->positionInNucleus.z(); }

void reactionProduct_SetFormationTime(reactionProduct* rp, dREAL aTime) { rp->formationTime = aTime; }

dREAL reactionProduct_GetFormationTime(reactionProduct* rp) { return rp->formationTime; }

void reactionProduct_HasInitialStateParton(reactionProduct* rp, char /*G4bool*/ aFlag) { rp->hasInitialStateParton = aFlag; }

char /*G4bool*/ reactionProduct_HasInitialStateParton(reactionProduct* rp) {
	return rp->hasInitialStateParton;
}

