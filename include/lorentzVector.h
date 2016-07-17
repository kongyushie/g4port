//source/global/HEPGeometry/include/G4LorentzVector.hh
//source/externals/clhep/include/CLHEP/Vector/LorentzVector.h
#pragma once
 
#include "debug.h"
#include "type.h"
#include "threeVector.h"

enum { lorentzVector_X=0, 
	   lorentzVector_Y=1, 
	   lorentzVector_Z=2, 
	   lorentzVector_T=3, 
	   lorentzVector_NUM_COORDINATES=4, 
	   lorentzVector_SIZE=lorentzVector_NUM_COORDINATES };
  // Safe indexing of the coordinates when using with matrices, arrays, etc.
  // (BaBar)

typedef struct lorentzVector{
	threeVector pp;
	dREAL  ee;

	/*DLL_API static*/ dREAL tolerance;
	/*DLL_API static*/ dREAL metric;
} lorentzVector;
