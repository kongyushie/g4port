// source/materials/include/G4Isotope.hh
#ifndef ISOTOPE
#define ISOTOPE 1
#include "type.h"
#include "isotopeVector.h"
// class description
//
// An isotope is a chemical isotope defined by its name,
//                                                 Z: atomic number
//                                                 N: number of nucleons
//                                                 A: mass of a mole (optional)
//                                                 m: isomer state (optional)
// If A is not defined it is taken from Geant4 database
//
// The class contains as a private static member the table of defined
// isotopes (an ordered vector of isotopes).
//
// Isotopes can be assembled into elements via the G4Element class.

typedef struct isotope{
	const char* fName;         // name of the Isotope
	int    fZ;                 // atomic number
	int    fN;                 // number of nucleons
	dREAL  fA;                 // atomic mass of a mole
	int    fm;                 // isomer level

//  int    fCountUse;          // nb of elements which use this isotope

//	Lawrence Cheng; 20160424
//	static G4IsotopeTable theIsotopeTable;
//	size_t fIndexInTable;        // index in the Isotope table
} isotope;

/* function in src/g4/isotope.c */
// Make an isotope
extern void constructIsotope(isotope* iso,
							 const char* name,	//its name
							 int     z,			//atomic number
							 int     n,			//number of nucleons
							 dREAL 	 a,			//mass of mole
							 int     m);		//isomer level


// Retrieval methods
extern const char* isotope_GetName(isotope* iso);

// Atomic number
extern int isotope_GetZ(isotope* iso);

// Number of nucleous
extern int isotope_GetN(isotope* iso);

// Atomic mass of mole in Geant4 units with electron shell
extern int isotope_GetA(isotope* iso);

// Isomer level
extern int isotope_Getm(isotope* iso);

// G4int GetCountUse() const {return fCountUse;}

//    static G4Isotope* GetIsotope(const G4String& name, G4bool warning=false);
//    static const G4IsotopeTable* GetIsotopeTable();

//    static size_t GetNumberOfIsotopes();
//    size_t GetIndex() const {return fIndexInTable;}

//    friend std::ostream& operator<<(std::ostream&, const G4Isotope*);
//    friend std::ostream& operator<<(std::ostream&, const G4Isotope&);
//    friend std::ostream& operator<<(std::ostream&, G4IsotopeTable);

//    G4int operator==(const G4Isotope&) const;
//    G4int operator!=(const G4Isotope&) const;

//    G4Isotope(__void__&);
      // Fake default constructor for usage restricted to direct object
      // persistency for clients requiring preallocation of memory for
      // persistifiable objects.

extern void isotope_SetName(isotope* iso, const char* name);
  //    void increaseCountUse()  {fCountUse++;}
  //  void decreaseCountUse()  {fCountUse--;}

#endif /* isotope */
