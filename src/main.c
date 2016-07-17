#include <stdio.h>
#include "crossSectionDataStore.h"
#include "dynamicParticle.h"
#include "particleDefinition.h"
#include "element.h"
#include "elementVector.h"
#include "nucleus.h"
int main()
{
	//material
	material mat;
	mat.fName = "G4_Al";
	mat.fChemicalFormula = "";
	mat.fDensity = 1.6845834537578004e+19;
	mat.fState = kStateSolid;
	mat.fTemp = 293.14999999999998;
	mat.fPressure = 632420964.9944762;
	mat.maxNbComponents = 1;
	mat.fArrayLength = 1;
	mat.fNumberOfElements = 1;
	mat.fNumberOfComponents = 1;
	element ele;
	ele.fName = "Al";
	ele.fSymbol = "Al";
	ele.fZeff = 13;
	ele.fNeff = 26.9815;
	elementVector eleVec;
	eleVec.v[0] = &ele;
	mat.theElementVector = &eleVec; 
	//dynamicParticle
	dynamicParticle dp;
	particleDefinition pd;
	pd.theParticleName = "photon";
	pd.thePDGMass = 938.27201300000002;
	pd.thePDGWidth = 0;
	pd.thePDGCharge = 1;
	pd.thePDGiSpin = 1;
	dp.theParticleDefinition = &pd;
	//crossSectionDataStore
	crossSectionDataStore XSDS;
	//nucleus
	nucleus target;
	crossSectionDataStore_SampleZandA(
			&XSDS,
			&dp,
			&mat,
			&target);

	printf("I love coding\n");



}
