#include <stdio.h>
#include "isotope.h"
#include "debug.h"


//source/materials/include/G4NistManager.hh, G4NistElementBuilder.hh
//static G4double G4NistElementBuilder::GetAtomicMass(G4int Z, G4int N) const
//{
//  G4double mass = 0.0;
//  G4int i = N - nFirstIsotope[Z];
//  if(i >= 0 && i <nIsotopes[Z]) {
//    mass = massIsotopes[i + idxIsotopes[Z]] +
//      Z*CLHEP::electron_mass_c2 - bindingEnergy[Z];
//  }
//  return mass;
//}

void constructIsotope(isotope* iso, const char* name, int z,	int n, dREAL a, int m){

	iso->fName = name;
	iso->fZ = z;
	iso->fN = n;
	iso->fA = a;
	iso->fm = m;
	//, fCountUse(0)

//  if (Z<1) {
//	G4ExceptionDescription ed;
//	ed << "Wrong Isotope " << Name << " Z= " << Z << G4endl;
//	G4Exception ("G4Isotope::G4Isotope()", "mat001", FatalException, ed);
//  }
//  if (N<Z) {
//	G4ExceptionDescription ed;
//	ed << "Wrong Isotope " << Name << " Z= " << Z << " > N= " << N << G4endl;
//	G4Exception ("G4Isotope::G4Isotope()", "mat002", FatalException, ed);
//  }

	if (a <=0.0) {
		BLURT;
		printf("unhandled situation: \"a <=0.0\"\n");
		exit(0);
//		iso.fA = (G4NistManager::Instance()->GetAtomicMass(Z,N))*g/(mole*amu_c2);
	}

//	theIsotopeTable.push_back(this);
//	fIndexInTable = theIsotopeTable.size() - 1;
}

const char* isotope_GetName(isotope* iso){return iso->fName;}
int isotope_GetZ(isotope* iso){return iso->fZ;}
int isotope_GetN(isotope* iso){return iso->fN;}
int isotope_GetA(isotope* iso){return iso->fA;}
int isotope_Getm(isotope* iso){return iso->fm;}
void isotope_SetName(isotope* iso, const char* name) {iso->fName=name;}

