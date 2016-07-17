// source/materials/include/G4Material.hh
#ifndef MATERIAL
#define MATERIAL 1
#include "debug.h"
#include "type.h"
#include <string.h>
#include "element.h"
#include "elementVector.h"
#include "physicalConstant.h"
#include <stdbool.h>
typedef enum MaterialState {
	kStateUndefined = 0,
	kStateSolid,
	kStateLiquid,
	kStateGas
} MaterialState;

static const dREAL NTP_Temperature = 293.15;

typedef struct material {
	//
	  // Basic data members ( To define a material)
	  //
	  char*         fName;                 // Material name
	  char*         fChemicalFormula;      // Material chemical formula
	  dREAL         fDensity;              // Material density

	  MaterialState   fState;                // Material state (determined
	                                          // internally based on density)
	  dREAL         fTemp;                 // Temperature (defaults: STP)
	  dREAL         fPressure;             // Pressure    (defaults: STP)

	  int            maxNbComponents;       // totalNbOfComponentsInTheMaterial
	  int            fArrayLength;          // the length of fAtomsVector
	  int            fNumberOfComponents;   // Nb of components declared so far

	  int            fNumberOfElements;     // Nb of Elements in the material
	  elementVector* theElementVector;      // vector of constituent Elements
//	  typedef std::vector<G4Element*> G4ElementVector;
	  dREAL*        fMassFractionVector;   // composition by fractional mass
	  int*           fAtomsVector;          // composition by atom count

//	  G4MaterialPropertiesTable* fMaterialPropertiesTable;

//	  static G4MaterialTable theMaterialTable;       // the material table
//	  size_t fIndexInTable;                   // the position in the table

	  //
	  // Derived data members (computed from the basic data members)
	  //
	  // some general atomic properties

	  dREAL* VecNbOfAtomsPerVolume;        // vector of nb of atoms per volume
	  dREAL  TotNbOfAtomsPerVolume;        // total nb of atoms per volume
	  dREAL  TotNbOfElectPerVolume;        // total nb of electrons per volume
	  dREAL  fRadlen;                      // Radiation length
	  dREAL  fNuclInterLen;                // Nuclear interaction length

//	  G4IonisParamMat* fIonisation;           // ionisation parameters
//	  G4SandiaTable*   fSandiaTable;          // Sandia table

	  // utilities
	  //
	  const struct material* fBaseMaterial;        // Pointer to the base material
	  dREAL fMassOfMolecule; 		  // for materials built by atoms count
//	  std::map<G4Material*,G4double> fMatComponents; // for composites built via
	                                                 // AddMaterial()
} material;



//
  // Constructor to create a material from single element
  //
  extern void constructMaterial_SE(material* ma, const char* name,				//its name
								dREAL  z, 				//atomic number
								dREAL  a,					//mass of mole
								dREAL  density, 				//density
								MaterialState   state    /*= kStateUndefined*/,	//solid,gas
								dREAL  temp     /*= NTP_Temperature*/,	//temperature
								dREAL  pressure /*= CLHEP::STP_Pressure*/);	//pressure

  //
  // Constructor to create a material from a combination of elements
  // and/or materials subsequently added via AddElement and/or AddMaterial
  //
  extern void constructMaterial_ME(material* ma, const char* name,				//its name
									dREAL  density, 				//density
									int     nComponents,			//nbOfComponents
									MaterialState   state    /*= kStateUndefined*/,	//solid,gas
									dREAL  temp     /*= NTP_Temperature*/,	//temperature
									dREAL  pressure /*= CLHEP::STP_Pressure*/);	//pressure

  //
  // Constructor to create a material from the base material
  //
  extern void constructMaterial_BM(material* ma, const char* name,				//its name
								dREAL  density, 				//density
								const material* baseMaterial,			//base material
								MaterialState   state    /*= kStateUndefined*/,	//solid,gas
								dREAL  temp     /*= NTP_Temperature*/,	//temperature
								dREAL  pressure /*= CLHEP::STP_Pressure*/);	//pressure

  //
  // Add an element, giving number of atoms
  //
  extern void material_AddElement_i (material* ma, element* element,				//the element
						  int      nAtoms);				//nb of atoms in
		    						// a molecule
  //
  // Add an element or material, giving fraction of mass
  //
  extern void material_AddElement_f (material* ma, element* element ,				//the element
		  	  	  	  	  dREAL   fraction);			//fractionOfMass

  extern void material_AddMaterial(material* ma, material* material,			//the material
		  	  	  	  	  dREAL   fraction);			//fractionOfMass


  extern void material_SetChemicalFormula (material* ma, const char* chF);

  //
  // retrieval methods
  //
  extern const char* material_GetName(material* ma);
  extern const char* material_GetChemicalFormula(material* ma);
  extern dREAL material_GetDensity(material* ma);
  extern MaterialState  material_GetState(material* ma);
  extern dREAL material_GetTemperature(material* ma);
  extern dREAL material_GetPressure(material* ma);

  //number of elements constituing this material:
  extern int material_GetNumberOfElements(material* ma);

  //vector of pointers to elements constituing this material:
//  extern const G4ElementVector* material_GetElementVector(material* ma);

  //vector of fractional mass of each element:
  extern const
  dREAL* material_GetFractionVector(material* ma);

  //vector of atom count of each element:
  extern const
  int*    material_GetAtomsVector(material* ma);

  //return a pointer to an element, given its index in the material:
  extern const element* material_GetElement(material* ma, int iel);

  //vector of nb of atoms per volume of each element in this material:
  extern const dREAL* material_GetVecNbOfAtomsPerVolume(material* ma);
  //total number of atoms per volume:
  extern dREAL  material_GetTotNbOfAtomsPerVolume(material* ma);
  //total number of electrons per volume:
  extern dREAL  material_GetTotNbOfElectPerVolume(material* ma);

  //obsolete names (5-10-98) see the 2 functions above
  extern const dREAL* material_GetAtomicNumDensityVector(material* ma);
  extern dREAL material_GetElectronDensity(material* ma);

  // Radiation length:
  extern dREAL  material_GetRadlen(material* ma);

  // Nuclear interaction length
  extern dREAL material_GetNuclearInterLength(material* ma);

  // ionisation parameters:
//  extern G4IonisParamMat* material_GetIonisation(material* ma);

  // Sandia table:
//  extern G4SandiaTable*  material_GetSandiaTable(material* ma);

  // Base material:
  extern
  const material* material_GetBaseMaterial(material* ma);

  // material components:
//  extern const std::map<G4Material*,G4double>& material_GetMatComponents(material* ma);

  // for chemical compound
  extern
  dREAL material_GetMassOfMolecule(material* ma);

  // meaningful only for single material:
  extern dREAL material_GetZ(material* ma);
  extern dREAL material_GetA(material* ma);

  //the MaterialPropertiesTable (if any) attached to this material:
//  extern void material_SetMaterialPropertiesTable(material* ma, G4MaterialPropertiesTable* anMPT);

//  extern G4MaterialPropertiesTable* material_GetMaterialPropertiesTable(material* ma);

  // the static Table of Materials:
  //
//  static G4MaterialTable* material_GetMaterialTable(material* ma);

  static size_t material_GetNumberOfMaterials(material* ma);

  //the index of this material in the Table:
  extern size_t material_GetIndex(material* ma);

  //return  pointer to a material, given its name:
  static material* material_GetMaterial(material* ma, const char* name, char/*G4bool*/ warning/*=true*/);

  extern void material_SetName (material* ma, const char* name);

  extern void material_InitializePointers(material* ma);

  // Header routine for all derived quantities
  extern void material_ComputeDerivedQuantities(material* ma);

  // Compute Radiation length
  extern void material_ComputeRadiationLength(material* ma);

  // Compute Nuclear interaction length
  extern void material_ComputeNuclearInterLength(material* ma);

  // Copy pointers of base material
  extern void material_CopyPointersOfBaseMaterial(material* ma);

  extern bool OP_EQ_material(const material* a,const material* b);
  extern bool OP_ASSIGN_material(material* a,const material* b);
#endif /* material */
