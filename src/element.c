#include "element.h"
#include <math.h>

// Constructor to Build an element directly; no reference to isotopes
  //
/*
  void constructElement(element* ele,
						const char* name,		//its name
						const char* symbol,		//its symbol
						dREAL  zeff,		//atomic number
						dREAL  aeff){		//mass of mole
	ele->fName = name;
	ele->fSymbol = symbol;
	int iz = (int) zeff;
	if (iz < 1) {
	  BLURT;
	  printf("Fail to create G4Element %s Z= %f < 1 !\n",
			  name, zeff);
	  exit(0);
	}

	if (abs(zeff - iz) > perMillion) {
	  BLURT;
	  printf("G4Element Warning: %s Z= %f  A= %f\n",
			  name, zeff, aeff/(g/mole));
	  exit(0);
	}

	element_InitializePointers(ele);

	ele->fZeff   = zeff;
	ele->fAeff   = aeff;
	ele->fNeff   = ele->fAeff/(g/mole);

	if(ele->fNeff < 1.0) ele->fNeff = 1.0;

	if (ele->fNeff < zeff) {
	  BLURT;
	  printf("Fail to create G4Element %s with Z= %f  N= %f  N < Z is not allowed\n",
			  name, zeff, ele->fNeff);
	  exit(0);
	}

	//		  ele->fNbOfAtomicShells      = G4AtomicShells::GetNumberOfShells(iz);
	//		  ele->fAtomicShells          = new G4double[fNbOfAtomicShells];
	//		  ele->fNbOfShellElectrons    = new G4int[fNbOfAtomicShells];

	element_AddNaturalIsotopes(ele);

	for (int i=0;i<ele->fNbOfAtomicShells;i++)
	{
	  BLURT;
	  printf("Not implemented\n");
	  exit(0);
	//			  ele->fAtomicShells[i] = G4AtomicShells::GetBindingEnergy(iz, i);
	//			  ele->fNbOfShellElectrons[i] = G4AtomicShells::GetNumberOfElectrons(iz, i);
	}
	element_ComputeDerivedQuantities(ele);
  }

  //
  // Constructor to Build an element from isotopes via AddIsotope
  //
//FIXME: G4IsotopeVector
void constructElement_iso(element* ele,
							const char* name,		//its name
							const char* symbol,		//its symbol
							int nIsotopes){			//nb of isotopes
	ele->fName = name;
	ele->fSymbol = symbol;
	element_InitializePointers(ele);

	size_t n = sizeof(nIsotopes);

	if(0 >= nIsotopes) {
		BLURT;
		printf("Fail to create G4Element %s <%s> with %d isotopes\n",
					  name, symbol, nIsotopes);
		exit(0);
	} else {
		ele->theIsotopeVector         = new G4IsotopeVector(n,0);
		ele->fRelativeAbundanceVector = new G4double[nIsotopes];
	}
}

  //
  // Add an isotope to the element
  //
//FIXME: G4AtomicShells
void element_AddIsotope(element* ele,
							  isotope* isotope,			//isotope
							  dREAL   abundance){	//fraction of nb of
								//atomes per volume
	if ( ele->theIsotopeVector == 0) {
		BLURT;
		printf("Fail to add Isotope to G4Element %s with Z= %f N= %f\n",
				ele->fName, ele->fZeff, ele->fNeff);
		exit(0);
	    //return;
	  }
	  int iz = isotope_GetZ(isotope);

	  // filling ...
	  if ( ele->fNumberOfIsotopes < (int) ele->theIsotopeVector->size() ) {
	    // check same Z
	    if (ele->fNumberOfIsotopes==0) { ele->fZeff = G4double(iz); }
	    else if ((dREAL)iz != ele->fZeff) {
	    	BLURT;
			printf("Fail to add Isotope Z= %f to G4Element %s with different Z= %f %f\n",
					iz, ele->fName, ele->fZeff, ele->fNeff);
			exit(0);
//	      return;
	    }
	    //Z ok
	    ele->fRelativeAbundanceVector[ele->fNumberOfIsotopes] = abundance;
	    ele->(*theIsotopeVector)[ele->fNumberOfIsotopes] = isotope;
	    ++ele->fNumberOfIsotopes;

	  } else {
			BLURT;
			printf("Fail to add Isotope Z= %f to G4Element %s  - more isotopes than declaired \n",
				iz, ele->fName);
			exit(0);
//	    return;
	  }

	  // filled.
	  if ( ele->fNumberOfIsotopes == (int) ele->theIsotopeVector->size() ) {
	    G4double wtSum=0.0;
	    ele->fAeff = 0.0;
	    for (int i=0; i<ele->fNumberOfIsotopes; i++) {
	    	ele->fAeff +=  ele->fRelativeAbundanceVector[i]*ele->(*theIsotopeVector)[i]->GetA();
	      wtSum +=  ele->fRelativeAbundanceVector[i];
	    }
	    if(wtSum > 0.0) { ele->fAeff  /= wtSum; }
	    ele->fNeff   = ele->fAeff/(g/mole);

	    if(wtSum != 1.0) {
	      for(int i=0; i<ele->fNumberOfIsotopes; ++i) {
	    	  ele->fRelativeAbundanceVector[i] /= wtSum;
	      }
	    }

	    ele->fNbOfAtomicShells = G4AtomicShells::GetNumberOfShells(iz);
	    ele->fAtomicShells     = new dREAL[ele->fNbOfAtomicShells];
	    ele->fNbOfShellElectrons = new int[ele->fNbOfAtomicShells];

	    for ( int j = 0; j < ele->fNbOfAtomicShells; j++ )
	    {
	    	ele->fAtomicShells[j]       = G4AtomicShells::GetBindingEnergy(iz, j);
	    	ele->fNbOfShellElectrons[j] = G4AtomicShells::GetNumberOfElectrons(iz, j);
	    }
	    element_ComputeDerivedQuantities(ele);

	  }
}
*/
char* element_GetName(element* ele){return (*ele).fName;}
char* element_GetSymbol(element* ele){return (*ele).fSymbol;}

// Atomic number
dREAL element_GetZ(element* ele){return (*ele).fZeff;}
/*
// Atomic weight in atomic units
dREAL element_GetN(element* ele) {return ele->fNeff;}
dREAL element_GetAtomicMassAmu(element* ele) {return ele->fNeff;}

// Mass of a mole in Geant4 units for atoms with atomic shell
dREAL element_GetA(element* ele) {return ele->fAeff;}
*/
//int/*G4bool*/ element_GetNaturalAbundanceFlag(element* ele){
//	return ele->fNaturalAbundance;
//}

//void element_SetNaturalAbundanceFlag(element* ele, int/*G4bool*/ flag){
//	ele->fNaturalAbundance = flag;
//}
/*
//the number of atomic shells in this element:
//
int element_GetNbOfAtomicShells(element* ele) {return ele->fNbOfAtomicShells;}

//the binding energy of the shell, ground shell index=0
//
dREAL element_GetAtomicShell(element* ele, int i){
	if (i<0 || i>=ele->fNbOfAtomicShells) {
		BLURT;
		printf("Invalid argument %d in for G4Element %s with Z= %f and Nshells= %d \n",
			i, ele->fName, ele->fZeff, ele->fNbOfAtomicShells );
		exit(0);
//	    return 0.0;
	  }
	  return ele->fAtomicShells[i];
}

//the number of electrons at the shell, ground shell index=0
//
int element_GetNbOfShellElectrons(element* ele, int i){
	if (i<0 || i>=ele->fNbOfAtomicShells) {
		BLURT;
		printf("Invalid argument %d in for G4Element %s with Z= %f and Nshells= %d \n",
			i, ele->fName, ele->fZeff, ele->fNbOfAtomicShells );
		exit(0);
//	    return 0;
	  }
	  return ele->fNbOfShellElectrons[i];
}
*/
//number of isotopes constituing this element:
//
size_t element_GetNumberOfIsotopes(element* ele) {return (*ele).fNumberOfIsotopes;}

//vector of pointers to isotopes constituing this element:
//
isotopeVector* element_GetIsotopeVector(element* ele) {return (*ele).theIsotopeVector;}

//vector of relative abundance of each isotope:
//
dREAL* element_GetRelativeAbundanceVector(element* ele)
				{return (*ele).fRelativeAbundanceVector;}
/*
const isotope* element_GetIsotope(element* ele, int iso)
				{return (*theIsotopeVector)[iso];}

//the (static) Table of Elements:
//
//  static
G4ElementTable* element_GetElementTable(element* ele){
	return ele->theElementTable;
}

//  static
size_t element_GetNumberOfElements(element* ele){
	return theElementTable.size();
}

//the index of this element in the Table:
//
size_t element_GetIndex(element* ele) {return ele->fIndexInTable;}
*/
//return pointer to an element, given its name:
//

//  static
//element* element_GetElement(element* ele, char* name, int /*G4bool*/ warning/*=true*/){
//	// search the element by its name
//	  for (size_t J=0 ; J<ele->theElementTable.size() ; J++)
//	   {
//	     if (ele->theElementTable[J]->GetName() == name)
//	       return ele->theElementTable[J];
//	   }
//
//	  // the element does not exist in the table
//	  if (warning) {
//		    BLURT;
//			printf("\n---> warning from G4Element::GetElement(). The element: %s",	name);
//			printf(" does not exist in the table. Return NULL pointer.\n");
//			exit(0);
//	  }
//	  return 0;
//}
/*
//Coulomb correction factor:
//
dREAL element_GetfCoulomb(element* ele) {return ele->fCoulomb;}

//Tsai formula for the radiation length:
//
dREAL element_GetfRadTsai(element* ele) {return ele->fRadTsai;}

void element_SetName(element* ele, const char* name)  {ele->fName=name;}
*/
//int/*G4bool*/ element_GetNaturalAbundanceFlag(element* ele){
//  return ele->fNaturalAbundance;
//}

//void element_SetNaturalAbundanceFlag(element* ele, int/*G4bool*/ val){
//	ele->fNaturalAbundance = val;
//}
/*
void element_InitializePointers(element* ele){
//	ele->theIsotopeVector = 0;
	ele->fRelativeAbundanceVector = 0;
	ele->fAtomicShells = 0;
	ele->fNbOfShellElectrons = 0;
//	ele->fIonisation = 0;
	ele->fNumberOfIsotopes = 0;
	ele->fNaturalAbundance = 0;//false;

	  // add initialisation of all remaining members
	ele->fZeff = 0;
	ele->fNeff = 0;
	ele->fAeff = 0;
	ele->fNbOfAtomicShells = 0;
	ele->fIndexInTable = 0;
	ele->fCoulomb = 0.0;
	ele->fRadTsai = 0.0;
}

void element_ComputeDerivedQuantities(element* ele){
	 // some basic functions of the atomic number
	BLURT;
	printf("Not implemented\n");
	exit(0);

	  // Store in table
//	  theElementTable.push_back(this);
//	  fIndexInTable = theElementTable.size() - 1;

	  // Radiation Length
	  element_ComputeCoulombFactor(ele);
	  element_ComputeLradTsaiFactor(ele);

	  // parameters for energy loss by ionisation
//	  if (ele->fIonisation) { delete fIonisation; }
//	  ele->fIonisation = new G4IonisParamElm(fZeff);
}

void element_ComputeCoulombFactor(element* ele){
	//
	//  Compute Coulomb correction factor (Phys Rev. D50 3-1 (1994) page 1254)

	static const dREAL k1 = 0.0083 , k2 = 0.20206 ,k3 = 0.0020 , k4 = 0.0369 ;

	dREAL az2 = (fine_structure_const*ele->fZeff)*(fine_structure_const*ele->fZeff);
	dREAL az4 = az2 * az2;

	ele->fCoulomb = (k1*az4 + k2 + 1./(1.+az2))*az2 - (k3*az4 + k4)*az4;
}

void element_ComputeLradTsaiFactor(element* ele){
	//
	  //  Compute Tsai's Expression for the Radiation Length
	  //  (Phys Rev. D50 3-1 (1994) page 1254)

	  static const dREAL Lrad_light[]  = {5.31  , 4.79  , 4.74 ,  4.71} ;
	  static const dREAL Lprad_light[] = {6.144 , 5.621 , 5.805 , 5.924} ;

	  const dREAL logZ3 = log(ele->fZeff)/3.;

	  dREAL Lrad, Lprad;
	  int iz = (int)(ele->fZeff+0.5) - 1 ;
	  if (iz <= 3) { Lrad = Lrad_light[iz] ;  Lprad = Lprad_light[iz] ; }
	    else { Lrad = log(184.15) - logZ3 ; Lprad = log(1194.) - 2*logZ3;}

	  ele->fRadTsai = 4*alpha_rcl2*ele->fZeff*(ele->fZeff*(Lrad-ele->fCoulomb) + Lprad);
}

//FIXME: G4NistManager
void element_AddNaturalIsotopes(element* ele){
	int Z = (int)ele->fZeff;
	G4NistManager* nist = G4NistManager::Instance();
	int n = nist->GetNumberOfNistIsotopes(Z);
	int N0 = nist->GetNistFirstIsotopeN(Z);

	if("" == ele->fSymbol) {
	const std::vector<G4String> elmnames =
	  G4NistManager::Instance()->GetNistElementNames();
	if(Z < (int)elmnames.size()) { ele->fSymbol = elmnames[Z]; }
	else { ele->fSymbol = ele->fName; }
	}

	ele->fNumberOfIsotopes = 0;
	for(int i=0; i<n; ++i) {
	if(nist->GetIsotopeAbundance(Z, N0+i) > 0.0) { ++ele->fNumberOfIsotopes; }
	}
	ele->theIsotopeVector = new G4IsotopeVector((unsigned int)ele->fNumberOfIsotopes,0);
	ele->fRelativeAbundanceVector = new G4double[ele->fNumberOfIsotopes];
	int idx = 0;
	dREAL xsum = 0.0;
	for(int i=0; i<n; ++i) {
	int N = N0 + i;
	dREAL x = nist->GetIsotopeAbundance(Z, N);
	if(x > 0.0) {
	  std::ostringstream strm;
	  strm << ele->fSymbol << N;
	  (*theIsotopeVector)[idx] = new G4Isotope(strm.str(), Z, N, 0.0, 0);
	  ele->fRelativeAbundanceVector[idx] = x;
	  xsum += x;
	  ++idx;
	}
	}
	if(xsum != 0.0 && xsum != 1.0) {
	for(int i=0; i<idx; ++i) { ele->fRelativeAbundanceVector[i] /= xsum; }
	}
	ele->fNaturalAbundance = 1; //true;
}
*/
