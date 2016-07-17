#include <stdio.h>
#include <math.h>
#include "nucleus.h"
#include "physicalConstant.h"

/*
//FIXME: G4ThreeVector
void constructNucleus(nucleus* nuc){
	nuc->theA = 0;
	nuc->theZ = 0;
	nuc->aEff = 0.0;
	nuc->zEff = 0;
	nuc->pnBlackTrackEnergy = 0.0;
	nuc->dtaBlackTrackEnergy = 0.0;
	nuc->pnBlackTrackEnergyfromAnnihilation = 0.0;
	nuc->dtaBlackTrackEnergyfromAnnihilation = 0.0;
	nuc->excitationEnergy = 0.0;
	nuc->momentum = G4ThreeVector(0.,0.,0.);
	nuc->fermiMomentum = 1.52*hbarc/fermi;
	nuc->theTemp = 293.16*kelvin;
	nuc->fIsotope = 0;
}

//FIXME: G4ThreeVector
void constructNucleus_AZ(nucleus* nuc, const dREAL A, const dREAL Z){
	nucleusSetParameters( A, Z );
	nuc->pnBlackTrackEnergy = 0.0;
	nuc->dtaBlackTrackEnergy = 0.0;
	nuc->pnBlackTrackEnergyfromAnnihilation = 0.0;
	nuc->dtaBlackTrackEnergyfromAnnihilation = 0.0;
	nuc->excitationEnergy = 0.0;
	nuc->momentum = G4ThreeVector(0.,0.,0.);
	nuc->fermiMomentum = 1.52*hbarc/fermi;
	nuc->theTemp = 293.16*kelvin;
	nuc->fIsotope = 0;
}

//FIXME: G4Material, G4ThreeVector
void constructNucleus_Material(nucleus* nuc, const material aMaterial){
	nucleusChooseParameters( aMaterial );
	nuc->pnBlackTrackEnergy = 0.0;
	nuc->dtaBlackTrackEnergy = 0.0;
	nuc->pnBlackTrackEnergyfromAnnihilation = 0.0;
	nuc->dtaBlackTrackEnergyfromAnnihilation = 0.0;
	nuc->excitationEnergy = 0.0;
	nuc->momentum = G4ThreeVector(0.,0.,0.);
	nuc->fermiMomentum = 1.52*hbarc/fermi;
	nuc->theTemp = material_GetTemperature( &aMaterial);
	nuc->fIsotope = 0;
}

//FIXME: G4ElementVector, G4UniformRand
void nucleus_ChooseParameters(nucleus* nuc, const material aMaterial ){
	dREAL random = G4UniformRand();
	dREAL sum = material_GetTotNbOfAtomsPerVolume(&aMaterial);
	const G4ElementVector* theElementVector = material_GetElementVector(&aMaterial);
	dREAL running = 0;
	//  G4Element* element(0);
	element* element = (*theElementVector)[material_GetNumberOfElements(&aMaterial)-1];

  for (unsigned int i = 0; i < material_GetNumberOfElements(&aMaterial); ++i) {
	running += material_GetVecNbOfAtomsPerVolume(&aMaterial)[i];
	if (running > random*sum) {
	  element = (*theElementVector)[i];
	  break;
	}
  }

  if (element_GetNumberOfIsotopes(element) > 0) {
	  dREAL randomAbundance = G4UniformRand();
	  dREAL sumAbundance = element_GetRelativeAbundanceVector(element)[0];
	  unsigned int iso=0;
	while (iso < element_GetNumberOfIsotopes(element) &&
		   sumAbundance < randomAbundance) {
	  ++iso;
	  sumAbundance += element_GetRelativeAbundanceVector(element)[iso];
	}
	nuc.theA=isotope_GetN(element_GetIsotope(element, iso));
	nuc.theZ=isotope_GetZ(element_GetIsotope(element, iso));
	nuc.aEff=nuc.theA;
	nuc.zEff=nuc.theZ;
  } else {
	nuc.aEff = element_GetN(element);
	nuc.zEff = element_GetZ(element);
	nuc.theZ = (int) (nuc.zEff + 0.5);
	nuc.theA = (int) (nuc.aEff + 0.5);
  }
}

void nucleus_SetParameters(nucleus* nuc, const dREAL A, const dREAL Z ){
	nuc->theZ = (int) Z;
	nuc->theA = (int) A;
	if (nuc->theA<1 || nuc->theZ<0 || nuc->theZ>nuc->theA) {
		BLURT;
		printf("SetParameters called with non-physical parameters\n");
		exit(0);
//		throw G4HadronicException(__FILE__, __LINE__,
//			"G4Nucleus::SetParameters called with non-physical parameters");
	}
	nuc->aEff = A;  // atomic weight
	nuc->zEff = Z;  // atomic number
	nuc->fIsotope = 0;
}

int nucleus_GetA_asInt(nucleus* nuc){
	return nuc->theA;
}

int nucleus_GetN_asInt(nucleus* nuc){
	return nuc->theA-nuc->theZ;
}

int nucleus_GetZ_asInt(nucleus* nuc){
	return nuc->theZ;
}

const isotope* nucleus_GetIsotope(nucleus* nuc){
	return nuc->fIsotope;
}
*/
void nucleus_SetIsotope(nucleus* nuc, const isotope* iso){
	nuc->fIsotope = iso;
	if(iso) {
		(*nuc).theZ = isotope_GetZ(iso);
		(*nuc).theA = isotope_GetN(iso);
		(*nuc).aEff = nuc->theA;
		(*nuc).zEff = nuc->theZ;
	}
}
/*
dREAL nucleus_GetPNBlackTrackEnergy(nucleus* nuc){
	return nuc->pnBlackTrackEnergy;
}

dREAL nucleus_GetDTABlackTrackEnergy(nucleus* nuc){
	return nuc->dtaBlackTrackEnergy;
}

dREAL nucleus_GetAnnihilationPNBlackTrackEnergy(nucleus* nuc){
	return nuc->pnBlackTrackEnergyfromAnnihilation;
}

dREAL nucleus_GetAnnihilationDTABlackTrackEnergy(nucleus* nuc){
	return nuc->dtaBlackTrackEnergyfromAnnihilation;
}

dREAL nucleus_Cinema(nucleus* nuc, dREAL kineticEnergy ){
	// derived from original FORTRAN code CINEMA by H. Fesefeldt (14-Oct-1987)
	//
	// input: kineticEnergy (MeV)
	// returns modified kinetic energy (MeV)
	//
	static const dREAL expxu =  82.;           // upper bound for arg. of exp
	static const dREAL expxl = -expxu;         // lower bound for arg. of exp

	dREAL ek = kineticEnergy/GeV;
	dREAL ekLog = log( ek );
	dREAL aLog = log( nuc->aEff );
	dREAL em = fmin( 1.0, 0.2390 + 0.0408*aLog*aLog );
	dREAL temp1 = -ek * fmin( 0.15, 0.0019*aLog*aLog*aLog );
	dREAL temp2 = exp( fmax( expxl, fmin( expxu, -(ekLog-em)*(ekLog-em)*2.0 ) ) );
	dREAL result = 0.0;
	if( abs( temp1 ) < 1.0 )
	{
	  if( temp2 > 1.0e-10 )result = temp1*temp2;
	}
	else result = temp1*temp2;
	if( result < -ek )result = -ek;
	return result*GeV;
}

//FIXME: G4UniformRand
dREAL nucleus_EvaporationEffects(nucleus* nuc, dREAL kineticEnergy ){
	// derived from original FORTRAN code EXNU by H. Fesefeldt (10-Dec-1986)
	    //
	    // Nuclear evaporation as function of atomic number
	    // and kinetic energy (MeV) of primary particle
	    //
	    // returns kinetic energy (MeV)
	    //
	    if( nuc->aEff < 1.5 ){
	    	nuc->pnBlackTrackEnergy = nuc->dtaBlackTrackEnergy = 0.0;
	    	return 0.0;
	    }
	    dREAL ek = kineticEnergy/GeV;
	    fREAL ekin = fmin( 4.0, fmax( 0.1, ek ) );
	    const fREAL atno = fmin( 120., nuc->aEff );
	    const fREAL gfa = 2.0*((nuc->aEff-1.0)/70.)*exp(-(nuc->aEff-1.0)/70.);
	    //
	    // 0.35 value at 1 GeV
	    // 0.05 value at 0.1 GeV
	    //
	    fREAL cfa = fmax( 0.15, 0.35 + ((0.35-0.05)/2.3)*log(ekin) );
	    fREAL exnu = 7.716 * cfa * exp(-cfa)
	      * ((atno-1.0)/120.)*exp(-(atno-1.0)/120.);
	    fREAL fpdiv = fmax( 0.5, 1.0-0.25*ekin*ekin );
	    //
	    // pnBlackTrackEnergy  is the kinetic energy (in GeV) available for
	    //                     proton/neutron black track particles
	    // dtaBlackTrackEnergy is the kinetic energy (in GeV) available for
	    //                     deuteron/triton/alpha black track particles
	    //
	    nuc->pnBlackTrackEnergy = exnu*fpdiv;
	    nuc->dtaBlackTrackEnergy = exnu*(1.0-fpdiv);

	    if( (int) (nuc->zEff+0.1) != 82 )
	    {
	    	dREAL ran1 = -6.0;
	    	dREAL ran2 = -6.0;
	      for( int i=0; i<12; ++i )
	      {
	        ran1 += G4UniformRand();
	        ran2 += G4UniformRand();
	      }
	      nuc->pnBlackTrackEnergy *= 1.0 + ran1*gfa;
	      nuc->dtaBlackTrackEnergy *= 1.0 + ran2*gfa;
	    }
	    nuc->pnBlackTrackEnergy = fmax( 0.0, nuc->pnBlackTrackEnergy );
	    nuc->dtaBlackTrackEnergy = fmax( 0.0, nuc->dtaBlackTrackEnergy );
	    while( nuc->pnBlackTrackEnergy+nuc->dtaBlackTrackEnergy >= ek )
	    {
	    	nuc->pnBlackTrackEnergy *= 1.0 - 0.5*G4UniformRand();
	    	nuc->dtaBlackTrackEnergy *= 1.0 - 0.5*G4UniformRand();
	    }
	//    G4cout << "EvaporationEffects "<<kineticEnergy<<" "
	//           <<pnBlackTrackEnergy+dtaBlackTrackEnergy<<endl;
	    return (nuc->pnBlackTrackEnergy+nuc->dtaBlackTrackEnergy)*GeV;
}

//FIXME: G4UniformRand
dREAL nucleus_AnnihilationEvaporationEffects(nucleus* nuc, dREAL kineticEnergy, dREAL ekOrg){
		// Nuclear evaporation as a function of atomic number and kinetic
	    // energy (MeV) of primary particle.  Modified for annihilation effects.
	    //
	    if( nuc->aEff < 1.5 || ekOrg < 0.)
	    {
	    	nuc->pnBlackTrackEnergyfromAnnihilation = 0.0;
	    	nuc->dtaBlackTrackEnergyfromAnnihilation = 0.0;
	      return 0.0;
	    }
	    dREAL ek = kineticEnergy/GeV;
	    fREAL ekin = fmin( 4.0, fmax( 0.1, ek ) );
	    const fREAL atno = fmin( 120., nuc->aEff );
	    const fREAL gfa = 2.0*((nuc->aEff-1.0)/70.)* exp(-(nuc->aEff-1.0)/70.);

	    fREAL cfa = fmax( 0.15, 0.35 + ((0.35-0.05)/2.3)*log(ekin) );
	    fREAL exnu = 7.716 * cfa * exp(-cfa)
	      * ((atno-1.0)/120.)*exp(-(atno-1.0)/120.);
	    fREAL fpdiv = fmax( 0.5, 1.0-0.25*ekin*ekin );

	    nuc->pnBlackTrackEnergyfromAnnihilation = exnu*fpdiv;
	    nuc->dtaBlackTrackEnergyfromAnnihilation = exnu*(1.0-fpdiv);

	    dREAL ran1 = -6.0;
	    dREAL ran2 = -6.0;
	    for( int i=0; i<12; ++i ) {
	      ran1 += G4UniformRand();
	      ran2 += G4UniformRand();
	    }
	    nuc->pnBlackTrackEnergyfromAnnihilation *= 1.0 + ran1*gfa;
	    nuc->dtaBlackTrackEnergyfromAnnihilation *= 1.0 + ran2*gfa;

	    nuc->pnBlackTrackEnergyfromAnnihilation = fmax( 0.0, nuc->pnBlackTrackEnergyfromAnnihilation);
	    nuc->dtaBlackTrackEnergyfromAnnihilation = fmax( 0.0, nuc->dtaBlackTrackEnergyfromAnnihilation);
	    dREAL blackSum = nuc->pnBlackTrackEnergyfromAnnihilation+nuc->dtaBlackTrackEnergyfromAnnihilation;
	    if (blackSum >= ekOrg/GeV) {
	    	nuc->pnBlackTrackEnergyfromAnnihilation *= ekOrg/GeV/blackSum;
	    	nuc->dtaBlackTrackEnergyfromAnnihilation *= ekOrg/GeV/blackSum;
	    }

	    return (nuc->pnBlackTrackEnergyfromAnnihilation+nuc->dtaBlackTrackEnergyfromAnnihilation)*GeV;
}

//FIXME: G4NucleiProperties
dREAL nucleus_AtomicMass(nucleus* nuc, const dREAL A, const dREAL Z ){
	// Now returns (atomic mass - electron masses)
	    return G4NucleiProperties::GetNuclearMass(A, Z);
}

//FIXME: G4NucleiProperties
dREAL nucleus_AtomicMass(nucleus* nuc, const int A, const int Z ){
	// Now returns (atomic mass - electron masses)
	    return G4NucleiProperties::GetNuclearMass(A, Z);
}

//FIXME: G4RandGauss
dREAL nucleus_GetThermalPz(nucleus* nuc, const dREAL mass, const dREAL temp ){
	dREAL result = G4RandGauss::shoot();
	    result *= sqrt(k_Boltzmann*temp*mass); // Das ist impuls (Pz),
	                                           // nichtrelativistische rechnung
	                                           // Maxwell verteilung angenommen
	return result;
}

//FIXME: G4UniformRand
dynamicParticle *nucleus_ReturnTargetParticle(nucleus* nuc){
	// choose a proton or a neutron as the target particle

	//FIXME: targetParticle is a local variable !!
	//TODO:  sol: use malloc to allocate it at heap
	dynamicParticle targetParticle;
	if( G4UniformRand() < nuc->zEff/nuc->aEff )
	  dynamicParticle_SetDefinition(&targetParticle, G4Proton::Proton() );
	else
	  dynamicParticle_SetDefinition(&targetParticle, G4Neutron::Neutron() );
	return targetParticle;
}
*/
//FIXME: G4Neutron
/*
//reactionProduct nucleus_GetThermalNucleus(nucleus* nuc, dREAL aMass, dREAL temp=-1){
reactionProduct nucleus_GetThermalNucleus(nucleus* nuc, dREAL aMass, dREAL temp){
	dREAL currentTemp = temp;
	if(currentTemp < 0) currentTemp = nuc->theTemp;
	reactionProduct theTarget;
	reactionProduct_SetMass(&theTarget, targetMass*G4Neutron::Neutron()->GetPDGMass());
	dREAL px, py, pz;
	px = GetThermalPz(reactionProduct_GetMass(&theTarget), currentTemp);
	py = GetThermalPz(reactionProduct_GetMass(&theTarget), currentTemp);
	pz = GetThermalPz(reactionProduct_GetMass(&theTarget), currentTemp);
	reactionProduct_SetMomentum(&theTarget, px, py, pz);
	dREAL tMom = sqrt(px*px+py*py+pz*pz);
	dREAL tEtot = sqrt((tMom+reactionProduct_GetMass(&theTarget))*
						  (tMom+reactionProduct_GetMass(&theTarget))-
						  2.*tMom*reactionProduct_GetMass(&theTarget));
	if(1-tEtot/reactionProduct_GetMass(&theTarget)>0.001)
	{
		reactionProduct_SetTotalEnergy(&theTarget, tEtot);
	}
	else
	{
		reactionProduct_SetKineticEnergy(&theTarget, tMom*tMom/(2.*reactionProduct_GetMass(&theTarget)));
	}
	return theTarget;
}
*/
//FIXME: G4Neutron
/*
reactionProduct nucleus_GetBiasedThermalNucleus(nucleus* nuc,
		dREAL aMass, threeVector aVelocity, dREAL temp){
	dREAL velMag = aVelocity.mag();
	  reactionProduct result;
	  dREAL value = 0;
	  dREAL random = 1;
	  dREAL norm = 3. * sqrt(k_Boltzmann*temp*aMass*G4Neutron::Neutron()->GetPDGMass());
	  norm /= G4Neutron::Neutron()->GetPDGMass();
	  norm *= 5.;
	  norm += velMag;
	  norm /= velMag;
	  while(value/norm<random)
	  {
	     result = GetThermalNucleus(aMass, temp);
	     G4ThreeVector targetVelocity = 1./result.GetMass()*result.GetMomentum();
	     value = (targetVelocity+aVelocity).mag()/velMag;
	     random = G4UniformRand();
	  }
	  return result;
}

// ******************  methods introduced by ChV ***********************
// return excitation Energy
dREAL nucleus_GetEnergyDeposit(nucleus* nuc) {
	return nuc->excitationEnergy;
}
// excitation Energy...
void nucleus_AddExcitationEnergy(nucleus* nuc, dREAL anEnergy){
	nuc->excitationEnergy+=anEnergy;
}


//FIXME: G4ReactionProductVector, G4ThreeVector

// return fermi momentum
threeVector nucleus_GetFermiMomentum(nucleus* nuc){
		// chv: .. we assume zero temperature!

	    // momentum is equally distributed in each phasespace volume dpx, dpy, dpz.
	dREAL ranflat1=
	      G4RandFlat::shoot((G4double)0.,(G4double)fermiMomentum);
	dREAL ranflat2=
	      G4RandFlat::shoot((G4double)0.,(G4double)fermiMomentum);
	dREAL ranflat3=
	      G4RandFlat::shoot((G4double)0.,(G4double)fermiMomentum);
	dREAL ranmax = (ranflat1>ranflat2? ranflat1: ranflat2);
	    ranmax = (ranmax>ranflat3? ranmax : ranflat3);

	    // Isotropic momentum distribution
	    dREAL costheta = 2.*G4UniformRand() - 1.0;
	    dREAL sintheta = sqrt(1.0 - costheta*costheta);
	    dREAL phi = 2.0*pi*G4UniformRand();

	    dREAL pz=costheta*ranmax;
	    dREAL px=sintheta*cos(phi)*ranmax;
	    dREAL py=sintheta*sin(phi)*ranmax;

	    //FIXME: p is a local variable !!
	    threeVector p;
	    constructThreeVector(&p, px,py,pz);
	    return p;
}
//  final nucleus fragmentation. Return List of particles
// which should be used for further tracking.
G4ReactionProductVector* nucleus_Fragmentate(nucleus* nuc){
	// needs implementation!
	    return 0;//NULL;
}
// momentum of absorbed Particles ..
void nucleus_AddMomentum(nucleus* nuc, const threeVector aMomentum){
	nuc->momentum+=(aMomentum);
}
*/
