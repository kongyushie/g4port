// source/materials/include/G4Element.hh

// Class description:
//
// An element is a chemical element either directly defined in terms of
// its characteristics: its name, symbol,
//                      Z (effective atomic number)
//                      N (effective number of nucleons)
//                      A (effective mass of a mole)
// or in terms of a collection of constituent isotopes with specified
// relative abundance (i.e. fraction of nb of atoms per volume).
//
// Quantities, with physical meaning or not, which are constant in a given
// element are computed and stored here as Derived data members.
//
// The class contains as a private static member the table of defined
// elements (an ordered vector of elements).
//
// Elements can be assembled singly or in mixtures into materials used
// in volume definitions via the G4Material class.
//
// It is strongly recommended do not delete G4Element instance in the
// user code. All G4Elements will be automatically deleted at the end
// of Geant4 session
#pragma once
#include "debug.h"
#include "type.h"
#include "isotope.h"
#include "isotopeVector.h"
#include <string.h>

typedef struct element {
	//
	  // Basic data members (which define an Element)
	  //
	  char* fName;              // name
	  char* fSymbol;            // symbol
	  dREAL fZeff;              // Effective atomic number
	  dREAL fNeff;              // Effective number of nucleons
	  dREAL fAeff;              // Effective mass of a mole

	  int fNbOfAtomicShells;     // number  of atomic shells
	  dREAL* fAtomicShells ;    // Pointer to atomic shell binding energies
	  int* fNbOfShellElectrons;  // Pointer to the number of subshell electrons

	  // Isotope vector contains constituent isotopes of the element
	  int fNumberOfIsotopes;     // Number of isotopes added to the element
	  isotopeVector* theIsotopeVector;
	  dREAL* fRelativeAbundanceVector;     // Fraction nb of atomes per volume
	                                          // for each constituent

	  // Set up the static Table of Elements
//	  static G4ElementTable theElementTable;
	  size_t fIndexInTable;
	  int/*G4bool*/ fNaturalAbundance;

	  //
	  // Derived data members (computed from the basic data members)
	  //
	  dREAL fCoulomb;             // Coulomb correction factor
	  dREAL fRadTsai;             // Tsai formula for the radiation length
//	  G4IonisParamElm* fIonisation;  // Pointer to ionisation parameters
} element;


//
  // Constructor to Build an element directly; no reference to isotopes
  //
  extern void constructElement(element* ele,
								const char* name,		//its name
								const char* symbol,		//its symbol
								dREAL  Zeff,		//atomic number
								dREAL  Aeff);		//mass of mole

  //
  // Constructor to Build an element from isotopes via AddIsotope
  //
  extern void constructElement_iso(element* ele,
								const char* name,		//its name
								const char* symbol,		//its symbol
								int nbIsotopes);			//nb of isotopes

  //
  // Add an isotope to the element
  //
  extern void element_AddIsotope(element* ele,
		  	  	  	  	  	  	  isotope* isotope,			//isotope
		  	  	  	  	  	  	  dREAL   RelativeAbundance);	//fraction of nb of
                  					//atomes per volume
//  virtual ~G4Element();

  //
  // Retrieval methods
  //
  extern char* element_GetName(element* ele);
  extern char* element_GetSymbol(element* ele);

  // Atomic number
  extern dREAL element_GetZ(element* ele);

  // Atomic weight in atomic units
  extern dREAL element_GetN(element* ele);
  extern dREAL element_GetAtomicMassAmu(element* ele);

  // Mass of a mole in Geant4 units for atoms with atomic shell
  extern dREAL element_GetA(element* ele);

  extern int/*G4bool*/   element_GetNaturalAbundanceFlag(element* ele);

  extern void     element_SetNaturalAbundanceFlag(element* ele, int/*G4bool*/);

  //the number of atomic shells in this element:
  //
  extern int element_GetNbOfAtomicShells(element* ele);

  //the binding energy of the shell, ground shell index=0
  //
  extern dREAL element_GetAtomicShell(element* ele, int index);

  //the number of electrons at the shell, ground shell index=0
  //
  extern int element_GetNbOfShellElectrons(element* ele, int index);

  //number of isotopes constituing this element:
  //
  extern size_t element_GetNumberOfIsotopes(element* ele);

  //vector of pointers to isotopes constituing this element:
  //
//  extern G4IsotopeVector* element_GetIsotopeVector(element* ele);

  //vector of relative abundance of each isotope:
  //
  extern dREAL* element_GetRelativeAbundanceVector(element* ele);

  extern const isotope* element_GetIsotope(element* ele, int iso);

  //the (static) Table of Elements:
  //
//  static
//  extern G4ElementTable* element_GetElementTable(element* ele);

//  static
  extern size_t element_GetNumberOfElements(element* ele);

  //the index of this element in the Table:
  //
  extern size_t element_GetIndex(element* ele);

  //return pointer to an element, given its name:
  //

//  static
  extern element* element_GetElement(element* ele, char* name, int /*G4bool*/ warning/*=true*/);

  //Coulomb correction factor:
  //
  extern dREAL element_GetfCoulomb(element* ele);

  //Tsai formula for the radiation length:
  //
  extern dREAL element_GetfRadTsai(element* ele);

  //pointer to ionisation parameters:
  //
//  extern G4IonisParamElm* element_GetIonisation() const {return fIonisation;}

//  // printing methods
//  //
//  friend std::ostream& operator<<(std::ostream&, const G4Element*);
//  friend std::ostream& operator<<(std::ostream&, const G4Element&);
//  friend std::ostream& operator<<(std::ostream&, G4ElementTable);

//public:  // without description
//
//  G4int operator==(const G4Element&) const;
//  G4int operator!=(const G4Element&) const;
//
//  G4Element(__void__&);
//    // Fake default constructor for usage restricted to direct object
//    // persistency for clients requiring preallocation of memory for
//    // persistifiable objects.

  extern void element_SetName(element* ele, const char* name);

//private:
//
//  G4Element(G4Element&);
//  const G4Element & operator=(const G4Element&);

  extern void element_InitializePointers(element* ele);
  extern void element_ComputeDerivedQuantities(element* ele);
  extern void element_ComputeCoulombFactor(element* ele);
  extern void element_ComputeLradTsaiFactor(element* ele);
  extern void element_AddNaturalIsotopes(element* ele);

extern int/*G4bool*/ element_GetNaturalAbundanceFlag(element* ele);

extern void element_SetNaturalAbundanceFlag(element* ele, int/*G4bool*/ val);

