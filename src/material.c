#include <stdio.h>
#include <math.h>
#include "material.h"
#include "elementVector.h"
#include <stdbool.h>

// Constructor to create a material from single element
//
//void constructMaterial_SE(material* ma, const char* name,				//its name
//							dREAL  z, 				//atomic number
//							dREAL  a,					//mass of mole
//							dREAL  density, 				//density
//							MaterialState   state    /*= kStateUndefined*/,	//solid,gas
//							dREAL  temp     /*= NTP_Temperature*/,	//temperature
//							dREAL  pressure /*= CLHEP::STP_Pressure*/){	//pressure
//	ma->fName = name;
//	InitializePointers();

//	  if (density < universe_mean_density){
//		  printf("G4Material WARNING: define a material with density=0 is not allowed. \n");
//		  printf("The material %s will be constructed with the default minimal density: %f g/cm3\n",
//				  name, universe_mean_density/(g/cm3));
//		  density = universe_mean_density;
//		}

//	  ma->fDensity  = density;
//	  ma->fState    = state;
//	  ma->fTemp     = temp;
//	  ma->fPressure = pressure;

	  // Initialize theElementVector allocating one
	  // element corresponding to this material
//	  ma->maxNbComponents        = ma->fNumberOfComponents = ma->fNumberOfElements = 1;
//	  ma->fArrayLength           = ma->maxNbComponents;
//	  ma->theElementVector       = new G4ElementVector();

//	  const char* elmnames = G4NistManager::Instance()->GetNistElementNames();
//	  char* enam, snam;
//	  int iz = G4lrint(z);
//	  if(iz < strlen(elmnames)/*(int)elmnames.size()*/) {
//		snam = elmnames[iz];
//		enam = snam;
//	  } else {
//		enam = "ELM_" + name;
//		snam = name;
//	  }
//	  theElementVector->push_back(new G4Element(enam, snam, z, a));

//	  ma->fMassFractionVector    = new G4double[1];
//	  ma->fMassFractionVector[0] = 1. ;
//	  ma->fMassOfMolecule        = a/Avogadro;

//	  if (ma->fState == kStateUndefined)
//		{
//		  if (ma->fDensity > kGasThreshold) { ma->fState = kStateSolid; }
//		  else                          	{ ma->fState = kStateGas; }
//		}

//	  ComputeDerivedQuantities();
//}

//
// Constructor to create a material from a combination of elements
// and/or materials subsequently added via AddElement and/or AddMaterial
//
//void constructMaterial_ME(material* ma, const char* name,				//its name
//								dREAL  density, 				//density
//								int     nComponents,			//nbOfComponents
//								MaterialState   state    /*= kStateUndefined*/,	//solid,gas
//								dREAL  temp     /*= NTP_Temperature*/,	//temperature
//								dREAL  pressure /*= CLHEP::STP_Pressure*/){	//pressure
//	ma->fName = name;
//	InitializePointers();

//	  if (density < universe_mean_density){
//		  printf("G4Material WARNING: define a material with density=0 is not allowed. \n");
//		  printf("The material %s will be constructed with the default minimal density: %f g/cm3\n",
//				  name, universe_mean_density/(g/cm3));
//	      density = universe_mean_density;
//	    }

//	  ma->fDensity  = density;
//	  ma->fState    = state;
//	  ma->fTemp     = temp;
//	  ma->fPressure = pressure;

//	  ma->maxNbComponents     = nComponents;
//	  ma->fArrayLength        = ma->maxNbComponents;
//	  ma->fNumberOfComponents = ma->fNumberOfElements = 0;
//	  ma->theElementVector    = new G4ElementVector();
//	  ma->theElementVector->reserve(maxNbComponents);

//	  if (ma->fState == kStateUndefined)
//	    {
//	      if (ma->Density > kGasThreshold) { ma->fState = kStateSolid; }
//	      else                          { ma->fState = kStateGas; }
//	    }

//}

//
// Constructor to create a material from the base material
//
//void constructMaterial_BM(material* ma, const char* name,				//its name
//							dREAL  density, 				//density
//							const material* baseMaterial,			//base material
//							MaterialState   state    /*= kStateUndefined*/,	//solid,gas
//							dREAL  temp     /*= NTP_Temperature*/,	//temperature
//							dREAL  pressure /*= CLHEP::STP_Pressure*/){	//pressure
//	ma->fName = name;
//	InitializePointers();

//	  if (density < universe_mean_density){
//		  printf("G4Material WARNING: define a material with density=0 is not allowed. \n");
//		  printf("The material %s will be constructed with the default minimal density: %f g/cm3\n",
//		  name, universe_mean_density/(g/cm3));
//	      density = universe_mean_density;
//	    }

//	  ma->fDensity  = density;
//	  ma->fState    = state;
//	  ma->fTemp     = temp;
//	  ma->fPressure = pressure;

//	  ma->fBaseMaterial = baseMaterial;
//	  ma->fChemicalFormula = ma->fBaseMaterial->GetChemicalFormula();
//	  ma->fMassOfMolecule  = ma->fBaseMaterial->GetMassOfMolecule();

//	  ma->fNumberOfElements = ma->fBaseMaterial->GetNumberOfElements();
//	  ma->maxNbComponents = ma->fNumberOfElements;
//	  ma->fNumberOfComponents = ma->fNumberOfElements;

//	  ma->fMaterialPropertiesTable = ma->fBaseMaterial->GetMaterialPropertiesTable();

//	  CopyPointersOfBaseMaterial();
//}

//
// Add an element, giving number of atoms
//
//void material_AddElement_i (material* ma, element* element,				//the element
//					  int      nAtoms){				//nb of atoms in
								// a molecule
	// initialization
//	  if ( ma->fNumberOfElements == 0 ) {
//		  ma->fAtomsVector        = new int   	[ma->fArrayLength];
//		  ma->fMassFractionVector = new dREAL	[ma->fArrayLength];
//	  }

	  // filling ...
//	  if ( ma->fNumberOfElements < ma->maxNbComponents ) {
//		  ma->theElementVector->push_back(element);
//		  ma->fAtomsVector[ma->fNumberOfElements] = nAtoms;
//		  ma->fNumberOfComponents = ++ma->fNumberOfElements;
//	  } else {
//		  BLURT;
//		  printf("G4Material::AddElement ERROR for %s  nElement= %d\n",
//				  ma->fName, ma->fNumberOfElements);
//		  exit(0);
//	  }
	  // filled.
//	  if ( ma->fNumberOfElements == ma->maxNbComponents ) {
	    // compute proportion by mass
//	    int i=0;
//	    dREAL Amol = 0.;
//	    for (i=0; i<ma->fNumberOfElements; ++i) {
//	    	dREAL w = ma->fAtomsVector[i]*(*theElementVector)[i]->GetA();
//	    	Amol += w;
//	    	ma->fMassFractionVector[i] = w;
//	    }
//	    for (i=0; i<ma->fNumberOfElements; ++i) {
//	    	ma->fMassFractionVector[i] /= Amol;
//	    }

//	    ma->fMassOfMolecule = Amol/Avogadro;
//	    ComputeDerivedQuantities();
//	  }
//}
//
// Add an element or material, giving fraction of mass
//
void material_AddElement_f (material* ma, element* element ,				//the element
					  dREAL   fraction){			//fractionOfMass
	if(fraction < 0.0 || fraction > 1.0) {
		BLURT;
		printf("G4Material::AddElement ERROR for %s and %s mass fraction = %d is wrong \n",
				ma->fName, element_GetName(element), fraction);
		exit(0);
	  }
	  // initialization
/*	  if (ma->fNumberOfComponents == 0) {
		  ma->fMassFractionVector = new dREAL [ma->fArrayLength];
		  ma->fAtomsVector        = new int   [ma->fArrayLength];
	  }
	  // filling ...
	  if (ma->fNumberOfComponents < ma->maxNbComponents) {
	    int el = 0;
	    // Loop checking, 07-Aug-2015, Vladimir Ivanchenko
	    while ((el<ma->fNumberOfElements)&&(element!=(*theElementVector)[el])) { ++el; }
	    if (el<ma->fNumberOfElements) ma->fMassFractionVector[el] += fraction;
	    else {
	    	ma->theElementVector->push_back(element);
	    	ma->fMassFractionVector[el] = fraction;
	      ++ma->fNumberOfElements;
	    }
	    ++ma->fNumberOfComponents;
	  } else {
		  BLURT;
		  printf("G4Material::AddElement ERROR for %s nElement = %d\n",
				  ma->fName, ma->fNumberOfElements);
		  exit(0);
	  }

	  // filled.
	  if (ma->fNumberOfComponents == ma->maxNbComponents) {

	    int i=0;
	    dREAL Zmol = 0., Amol = 0.;
	    // check sum of weights -- OK?
	    dREAL wtSum = 0.0;
	    for (i=0; i<ma->fNumberOfElements; ++i) {
	      wtSum += ma->fMassFractionVector[i];
	      Zmol +=  ma->fMassFractionVector[i]*(*theElementVector)[i]->GetZ();
	      Amol +=  ma->fMassFractionVector[i]*(*theElementVector)[i]->GetA();
	    }
	    if (fabs(1.-wtSum) > perThousand) {
	    	BLURT;
	    	printf("WARNING !! for %s sum of fractional masses %f is not 1 - results may be wrong\n",
	    			ma->fName, wtSum);
	    	exit(0);
	    }
	    for (i=0; i<ma->fNumberOfElements; ++i) {
	    	ma->fAtomsVector[i] =
		G4lrint(ma->fMassFractionVector[i]*Amol/(*theElementVector)[i]->GetA());
	    }

	    ComputeDerivedQuantities();
	  }
*/
}

/*void material_AddMaterial(material* ma, material* material,			//the material
					  dREAL   fraction){			//fractionOfMass
	if(fraction < 0.0 || fraction > 1.0) {
		BLURT;
		printf("G4Material::AddMaterial ERROR for %s and %s mass fraction = %d is wrong \n",
						ma->fName, material_GetName(material), fraction);
		exit(0);
	  }
	  // initialization
	  if (ma->fNumberOfComponents == 0) {
		  ma->fMassFractionVector = new dREAL	[ma->fArrayLength];
		  ma->fAtomsVector        = new int   	[ma->fArrayLength];
	  }

	  int nelm = material_GetNumberOfElements(material);

	  // arrays should be extended
	  if(nelm > 1) {
	    int nold    = ma->fArrayLength;
	    ma->fArrayLength += nelm - 1;
	    dREAL* v1 = new dREAL[ma->fArrayLength];
	    int* i1    = new int[ma->fArrayLength];
	    for(int i=0; i<nold; ++i) {
	      v1[i] = ma->fMassFractionVector[i];
	      i1[i] = ma->fAtomsVector[i];
	    }
//	    delete [] fAtomsVector;
//	    delete [] fMassFractionVector;
	    ma->fMassFractionVector = v1;
	    ma->fAtomsVector = i1;
	  }

	  // filling ...
	  if (ma->fNumberOfComponents < ma->maxNbComponents) {
	    for (int elm=0; elm<nelm; ++elm)
	      {
	        element* element = (*(material->GetElementVector()))[elm];
	        int el = 0;
		// Loop checking, 07-Aug-2015, Vladimir Ivanchenko
	        while ((el<ma->fNumberOfElements)&&(element!=(*theElementVector)[el])) el++;
	        if (el < ma->fNumberOfElements) ma->fMassFractionVector[el] += fraction
	                                          *(material_GetFractionVector(material))[elm];
	        else {
	        	ma->theElementVector->push_back(element);
	        	ma->fMassFractionVector[el] = fraction
		                                  *(material_GetFractionVector(material))[elm];
	          ++ma->fNumberOfElements;
	        }
	      }
	    ++ma->fNumberOfComponents;
	    ///store massFraction of material component
	    ma->fMatComponents[material] = fraction;

	  } else {
		  BLURT;
		  printf("G4Material::AddMaterial ERROR for %s nElement= %d\n",
				  ma->fName, ma->fNumberOfElements);
		  exit(0);
	  }

	  // filled.
	  if (ma->fNumberOfComponents == ma->maxNbComponents) {
	    int i=0;
	    dREAL Zmol = 0., Amol = 0.;
	    // check sum of weights -- OK?
	    dREAL wtSum = 0.0;
	    for (i=0; i<ma->fNumberOfElements; ++i) {
	      wtSum += ma->fMassFractionVector[i];
	      Zmol +=  ma->fMassFractionVector[i]*(*theElementVector)[i]->GetZ();
	      Amol +=  ma->fMassFractionVector[i]*(*theElementVector)[i]->GetA();
	    }
	    if (fabs(1.-wtSum) > perThousand) {
	    	printf("G4Material::AddMaterial WARNING !! for %s "
	    			"sum of fractional masses %f is not 1 - results may be wrong\n",
	    					  ma->fName, wtSum);
	    }
	    for (i=0; i<ma->fNumberOfElements; ++i) {
	    	ma->fAtomsVector[i] =
		G4lrint(ma->fMassFractionVector[i]*Amol/(*theElementVector)[i]->GetA());
	    }

	    ComputeDerivedQuantities();
	  }
}


void material_SetChemicalFormula (material* ma, const char* chF) {
	ma->fChemicalFormula=chF;
}

//
// retrieval methods
//
const char* material_GetName(material* ma){
	return ma->fName;
}

const char* material_GetChemicalFormula(material* ma) {
	return ma->fChemicalFormula;
}

dREAL material_GetDensity(material* ma){
	return ma->fDensity;
}
MaterialState  material_GetState(material* ma){
	return ma->fState;
}

dREAL material_GetTemperature(material* ma){
	return ma->fTemp;
}

dREAL material_GetPressure(material* ma){
	return ma->fPressure;
}
*/
//number of elements constituing this material:
int material_GetNumberOfElements(material* ma){
	return (*ma).fNumberOfElements;
}

////vector of pointers to elements constituing this material:
const elementVector* material_GetElementVector(material* ma) {return (*ma).theElementVector;}
/*
//vector of fractional mass of each element:
const dREAL* material_GetFractionVector(material* ma) {return ma->fMassFractionVector;}

//vector of atom count of each element:
const int*    material_GetAtomsVector(material* ma){return ma->fAtomsVector;}

////return a pointer to an element, given its index in the material:
//const element* material_GetElement(material* ma, int iel) {return (*theElementVector)[iel];}
*/
//vector of nb of atoms per volume of each element in this material:
const dREAL* material_GetVecNbOfAtomsPerVolume(material* ma) {return (*ma).VecNbOfAtomsPerVolume;}
/*
//total number of atoms per volume:
dREAL  material_GetTotNbOfAtomsPerVolume(material* ma) {return ma->TotNbOfAtomsPerVolume;}
//total number of electrons per volume:
dREAL  material_GetTotNbOfElectPerVolume(material* ma) {return ma->TotNbOfElectPerVolume;}

//obsolete names (5-10-98) see the 2 functions above
const dREAL* material_GetAtomicNumDensityVector(material* ma) {return ma->VecNbOfAtomsPerVolume;}
dREAL  material_GetElectronDensity(material* ma) {return ma->TotNbOfElectPerVolume;}

// Radiation length:
dREAL  material_GetRadlen(material* ma){return ma->fRadlen;}

// Nuclear interaction length
dREAL material_GetNuclearInterLength(material* ma) {return ma->fNuclInterLen;}

//// ionisation parameters:
//G4IonisParamMat* material_GetIonisation(material* ma) {return ma->fIonisation;}
//
//// Sandia table:
//G4SandiaTable*  material_GetSandiaTable(material* ma) {return ma->fSandiaTable;}
//
//// Base material:
//const material* material_GetBaseMaterial(material* ma) {
//	return ma->fBaseMaterial;
//}

// material components:
//const std::map<G4Material*,G4double>& material_GetMatComponents(material* ma){
//	return ma->fMatComponents;
//}

//// for chemical compound
//dREAL material_GetMassOfMolecule(material* ma) {
//	return ma->fMassOfMolecule;
//}

//FIXME: theElementVector, element
dREAL material_GetZ(material* ma){
// meaningful only for single material:
	if (ma->fNumberOfElements > 1) {
		BLURT;
		printf("G4Material ERROR in GetZ. The material: %s is a mixture.\n", ma->fName);
		exit(0);
	  }
	  return element_GetZ((*theElementVector)[0]);
}

//FIXME: theElementVector, element
dREAL material_GetA(material* ma){
	if (ma->fNumberOfElements > 1) {
		BLURT;
		printf("G4Material ERROR in GetA. The material: %s is a mixture.\n", ma->fName);
		exit(0);
	  }
	  return  element_GetA((*theElementVector)[0]);
}

////the MaterialPropertiesTable (if any) attached to this material:
//void material_SetMaterialPropertiesTable(material* ma, G4MaterialPropertiesTable* anMPT){
//	ma->fMaterialPropertiesTable = anMPT;
//}
//
//G4MaterialPropertiesTable* material_GetMaterialPropertiesTable(material* ma) {
//	return ma->fMaterialPropertiesTable;
//}
//
////the index of this material in the Table:
//size_t material_GetIndex(material* ma) {
//	return ma->fIndexInTable;
//}


void material_SetName (material* ma, const char* name){
	ma->fName=name;
}

//void material_InitializePointers(material* ma){
//	theElementVector         = 0;
//	  fMassFractionVector      = 0;
//	  fAtomsVector             = 0;
//	  fMaterialPropertiesTable = 0;
//
//	  VecNbOfAtomsPerVolume    = 0;
//	  fBaseMaterial            = 0;
//
//	  fChemicalFormula         = "";
//
//	  // initilized data members
//	  fDensity  = 0.0;
//	  fState    = kStateUndefined;
//	  fTemp     = 0.0;
//	  fPressure = 0.0;
//	  maxNbComponents     = 0;
//	  fArrayLength        = 0;
//	  TotNbOfAtomsPerVolume = 0;
//	  TotNbOfElectPerVolume = 0;
//	  fRadlen = 0.0;
//	  fNuclInterLen = 0.0;
//	  fMassOfMolecule = 0.0;
//
//	  fIonisation = 0;
//	  fSandiaTable = 0;
//
//	  // Store in the static Table of Materials
//	  fIndexInTable = theMaterialTable.size();
//	  for(size_t i=0; i<fIndexInTable; ++i) {
//	    if(theMaterialTable[i]->GetName() == fName) {
//	      G4cout << "G4Material WARNING: duplicate name of material "
//		     << fName << G4endl;
//	      break;
//	    }
//	  }
//	  theMaterialTable.push_back(this);
//}

// Header routine for all derived quantities
void material_ComputeDerivedQuantities(material* ma){
	BLURT;
	printf("Not Implemented\n");
	exit(0);

//	// Header routine to compute various properties of material.
//	//
//
//	// Number of atoms per volume (per element), total nb of electrons per volume
//	dREAL Zi, Ai;
//	ma->TotNbOfAtomsPerVolume = 0.;
//	if (ma->VecNbOfAtomsPerVolume) { delete [] ma->VecNbOfAtomsPerVolume; }
//	ma->VecNbOfAtomsPerVolume = new G4double[ma->fNumberOfElements];
//	ma->TotNbOfElectPerVolume = 0.;
//	for (int i=0; i<ma->fNumberOfElements; ++i) {
//	 Zi = (*theElementVector)[i]->GetZ();
//	 Ai = (*theElementVector)[i]->GetA();
//	 ma->VecNbOfAtomsPerVolume[i] = Avogadro*ma->fDensity*ma->fMassFractionVector[i]/Ai;
//	 ma->TotNbOfAtomsPerVolume += ma->VecNbOfAtomsPerVolume[i];
//	 ma->TotNbOfElectPerVolume += ma->VecNbOfAtomsPerVolume[i]*Zi;
//	}
//
//	ComputeRadiationLength();
//	ComputeNuclearInterLength();
//
//	if (ma->fIonisation) { delete fIonisation; }
//	ma->fIonisation  = new G4IonisParamMat(this);
//	if (ma->fSandiaTable) { delete fSandiaTable; }
//	ma->fSandiaTable = new G4SandiaTable(this);
}

//FIXME: DBL_MAX, theElementVector
void material_ComputeRadiationLength(material* ma){
// Compute Radiation length
	dREAL radinv = 0.0 ;
	  for (int i=0;i<ma->fNumberOfElements;++i) {
	     radinv += ma->VecNbOfAtomsPerVolume[i]*((*theElementVector)[i]->GetfRadTsai());
	  }
	  ma->fRadlen = (radinv <= 0.0 ? DBL_MAX : 1./radinv);
}

//FIXME: DBL_MAX, theElementVector
void material_ComputeNuclearInterLength(material* ma){
// Compute Nuclear interaction length
	static const dREAL lambda0  = 35*g/cm2;
	  static const dREAL twothird = 2.0/3.0;
	  dREAL NILinv = 0.0;
	  for (int i=0; i<ma->fNumberOfElements; ++i) {
	    int Z = G4lrint( (*theElementVector)[i]->GetZ());
	    dREAL A = (*theElementVector)[i]->GetN();
	    if(1 == Z) {
	      NILinv += ma->VecNbOfAtomsPerVolume[i]*A;
	    } else {
	      NILinv += ma->VecNbOfAtomsPerVolume[i]*G4Exp(twothird*G4Log(A));
	    }
	  }
	  NILinv *= amu/lambda0;
	  ma->fNuclInterLen = (NILinv <= 0.0 ? DBL_MAX : 1./NILinv);
}

//// Copy pointers of base material
//void material_CopyPointersOfBaseMaterial(material* ma){
//	G4double factor = fDensity/fBaseMaterial->GetDensity();
//	  TotNbOfAtomsPerVolume = factor*fBaseMaterial->GetTotNbOfAtomsPerVolume();
//	  TotNbOfElectPerVolume = factor*fBaseMaterial->GetTotNbOfElectPerVolume();
//
//	  theElementVector =
//	    const_cast<G4ElementVector*>(fBaseMaterial->GetElementVector());
//	  fMassFractionVector =
//	    const_cast<G4double*>(fBaseMaterial->GetFractionVector());
//	  fAtomsVector = const_cast<G4int*>(fBaseMaterial->GetAtomsVector());
//
//	  const G4double* v = fBaseMaterial->GetVecNbOfAtomsPerVolume();
//	  if (VecNbOfAtomsPerVolume)  { delete [] VecNbOfAtomsPerVolume; }
//	  VecNbOfAtomsPerVolume = new G4double[fNumberOfElements];
//	  for (G4int i=0; i<fNumberOfElements; ++i) {
//	    VecNbOfAtomsPerVolume[i] = factor*v[i];
//	  }
//	  fRadlen = fBaseMaterial->GetRadlen()/factor;
//	  fNuclInterLen = fBaseMaterial->GetNuclearInterLength()/factor;
//
//	  if (fIonisation) { delete fIonisation; }
//	  fIonisation = new G4IonisParamMat(this);
//
//	  fSandiaTable = fBaseMaterial->GetSandiaTable();
//	  fIonisation->SetMeanExcitationEnergy(fBaseMaterial->GetIonisation()->GetMeanExcitationEnergy());
//}
*/
bool OP_EQ_material(const material *a,const material *b)
{
	return ((*a).fTemp == (*b).fTemp) 
			&& ((*a).fPressure == (*b).fPressure) 
			&&((*a).maxNbComponents == (*b).maxNbComponents)	
			&&((*a).fArrayLength == (*b).fArrayLength) 
			&&((*a).fNumberOfComponents == (*b).fNumberOfComponents) 
			&&((*a).fNumberOfElements == (*b).fNumberOfElements);
}
bool OP_ASSIGN_material(material *a,const material *b)
{
	(*a).fTemp				 = (*b).fTemp; 
	(*a).fPressure			 = (*b).fPressure;
	(*a).maxNbComponents	 = (*b).maxNbComponents; 
	(*a).fArrayLength		 = (*b).fArrayLength;
	(*a).fNumberOfComponents = (*b).fNumberOfComponents;
	(*a).fNumberOfElements	 = (*b).fNumberOfElements;
		
}
