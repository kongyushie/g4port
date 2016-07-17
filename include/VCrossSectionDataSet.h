
#ifndef G4VCrossSectionDataSet_h
#define G4VCrossSectionDataSet_h 1

//#include "globals.hh"
#include "type.h"
#include "element.h"
#include "material.h"
#include "particleDefinition.h"
#include "element.h"
#include "isotope.h"
#include "dynamicParticle.h"
//#include "test.h"
#include <stdbool.h>
//#include "G4HadTmpUtil.hh"

//class G4DynamicParticle;
//class G4Isotope;
//class G4Material;
//class G4CrossSectionDataSetRegistry;

typedef struct VCrossSectionDataSet
{
	int verboseLevel;

	//CrossSectionDataSetRegistry* registry;

	double minKinEnergy;
	double maxKinEnergy;

	//G4string name;
	char* name;
}VCrossSectionDataSet;
//public: //with description

extern void constructVCrossSectionDataSet(VCrossSectionDataSet* self, const char* nam);

//G4VCrossSectionDataSet(const G4VCrossSectionDataSet&);

//G4VCrossSectionDataSet & operator=(const G4VCrossSectionDataSet &right);

extern void destructVCrossSectionDataSet(VCrossSectionDataSet* self);

//virtual func :  need run time test
extern
bool VCrossSectionDataSet_IsElementApplicable(
		dREAL self,
		const dynamicParticle* dp, 
		int Z, 
		const material* mat);

extern
bool G4VCrossSectionDataSet_IsIsoApplicable(const dynamicParticle* dp, int Z, int A,    
		const element* elm,
		const material* mat);


//G4double ComputeCrossSection(const G4DynamicParticle*, 
//		const G4Element*,
//		const G4Material* mat = 0);

//virtual
dREAL GetElementCrossSection(const dynamicParticle*, int Z,
		const material* mat);

//virtual
//G4double GetIsoCrossSection(const G4DynamicParticle*, G4int Z, G4int A,  
dREAL VCrossSectionDataSet_GetIsoCrossSection(
		dREAL self,
		const dynamicParticle*, int Z, int A,  
		//const isotope* iso = 0,
		const isotope* iso,
		//const element* elm = 0,
		const element* elm,
		//const material* mat = 0);
		const material* mat);


extern isotope* VCrossSectionDataSet_SelectIsotope(VCrossSectionDataSet* self, 
			const element* anElement, 
			dREAL kinEnergy);

//virtual
//void BuildPhysicsTable(const G4ParticleDefinition&);

//virtual
//void DumpPhysicsTable(const G4ParticleDefinition&);

//virtual void CrossSectionDescription(std::ostream&) const;

/*
virtual G4int GetVerboseLevel() const;

virtual void SetVerboseLevel(G4int value);


	inline G4double 
G4VCrossSectionDataSet::GetCrossSection(const G4DynamicParticle* dp, 
		const G4Element* elm,
		const G4Material* mat)
{
	return ComputeCrossSection(dp, elm, mat);
}


inline G4int G4VCrossSectionDataSet::GetVerboseLevel() const
{
	return verboseLevel;
}

inline void G4VCrossSectionDataSet::SetVerboseLevel(G4int value)
{
	verboseLevel = value;
}

inline void G4VCrossSectionDataSet::SetMinKinEnergy(G4double value)
{
	minKinEnergy = value;
}

inline G4double G4VCrossSectionDataSet::GetMinKinEnergy() const
{
	return minKinEnergy;
}

inline void G4VCrossSectionDataSet::SetMaxKinEnergy(G4double value)
{
	maxKinEnergy = value;
}

inline G4double G4VCrossSectionDataSet::GetMaxKinEnergy() const
{
	return maxKinEnergy;
}

inline const G4String& G4VCrossSectionDataSet::GetName() const
{
	return name;
}

inline void G4VCrossSectionDataSet::SetName(const G4String& nam)
{
	name = nam;
}
*/
#endif
