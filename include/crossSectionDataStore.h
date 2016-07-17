#pragma once
//#include "globals.hh"
//#include "G4VCrossSectionDataSet.hh"
//#include "FastPathHadronicCrossSection.hh"
#include "dynamicParticle.h"
//#include "G4PhysicsVector.hh"
#include "vector.h"
//#include "test.h"
#include "material.h"
#include "particleDefinition.h"
#include "nucleus.h" 
//class G4Nucleus;
//class G4ParticleDefinition;
//class G4Isotope;
//class G4Element;
//class G4Material;
//class G4NistManager;

//class G4CrossSectionDataStore
typedef struct crossSectionDataStore
{
	//G4NistManager* nist;

	//vector<G4VCrossSectionDataSet*> dataSetList;
	vector_VCrossSectionDataSet dataSetList;
	vector_dREAL xsecelm;
	vector_dREAL xseciso;

	const material* currentMaterial;
	const particleDefinition* matParticle;
	dREAL matKinEnergy;
	dREAL matCrossSection;

	const material* elmMaterial;
	const element* currentElement;
	const particleDefinition* elmParticle;
	dREAL elmKinEnergy;
	dREAL elmCrossSection;

	int nDataSetList;
	int verboseLevel;
	//private:
	//friend struct G4FastPathHadronicCrossSection::fastPathEntry;
	//The following method is called by the public one GetCrossSection(const G4DynamicParticle*, const G4Material*)
	//The third parameter is used to force the calculation of cross-sections skipping the fast-path mechanism
	
	//G4FastPathHadronicCrossSection_controlFlag fastPathFlags;
	//G4FastPathHadronicCrossSection_fastPathParameters fastPathParams;
	//Counters
	//G4FastPathHadronicCrossSection_getCrossSectionCount counters;
	//TODO: share this among threads
	//G4FastPathHadronicCrossSection_G4CrossSectionDataStore_Cache fastPathCache;
	//G4FastPathHadronicCrossSection_timing timing;
	//G4FastPathHadronicCrossSection_G4CrossSectionDataStore_Requests requests;
	
}crossSectionDataStore;
//public:

extern void constructcrossSectionDataStore(crossSectionDataStore* self);

extern void destructcrossSectionDataStore(crossSectionDataStore* self);


// Cross section per element is computed
// getXS(Dp, M ,E)
/*extern G4double G4CrossSectionDataStore_GetCrossSection_E(
		G4CrossSectionDataStore* self, 
		const G4DynamicParticle*, 
		const G4Element*, 
		const G4Material*);

// Cross section per isotope is computed
// getXS(Z,A)
extern G4double G4CrossSectionDataStore_GetCrossSection_ZA(
		G4CrossSectionDataStore* self,
		const G4DynamicParticle*, 
		G4int Z, 
		G4int A,
		const G4Isotope*,
		const G4Element*, 
		const G4Material*);
*/
// getXS(DP, M)
extern dREAL crossSectionDataStore_GetCrossSection(
		crossSectionDataStore* self,
		const dynamicParticle*, 
		const material*);
// Sample Z and A of a target nucleus and upload into G4Nucleus

extern element* crossSectionDataStore_SampleZandA(
		crossSectionDataStore* self, 
		const dynamicParticle* dp, 
		const material*,
		nucleus* target);

/*
// Initialisation before run
void BuildPhysicsTable(const G4ParticleDefinition&);

// Dump store to G4cout
//void DumpPhysicsTable(const G4ParticleDefinition&);

// Dump store as html
//void DumpHtml(const G4ParticleDefinition&, std::ofstream&) const;
//void PrintCrossSectionHtml(const G4VCrossSectionDataSet *cs) const;
//private:
*/
extern dREAL crossSectionDataStore_GetIsoCrossSection(
		crossSectionDataStore* self,
		const dynamicParticle*, int Z, int A,
		const isotope*,
		const element*, const material* aMaterial,
		int index);
// no define operator
//G4CrossSectionDataStore & operator=(const G4CrossSectionDataStore &right);
/*
void cnostructG4CrossSectionDataStore_G4CrossSectionDataStore(crossSectionDataStore* self,const crossSectionDataStore* ano);

//public:

inline const G4FastPathHadronicCrossSection_fastPathParameters*
   G4CrossSectionDataStore_GetFastPathParameters(G4CrossSectionDataStore* self) const { return &self.fastPathParams; }
inline const G4FastPathHadronicCrossSection_controlFlag*
   G4CrossSectionDataStore_GetFastPathControlFlags(G4CrossSectionDataStore* self) const { return &self.fastPathFlags; }

void DumpFastPath( const G4ParticleDefinition* , const G4Material* , std::ostream& os);
 
inline G4double G4CrossSectionDataStore_GetCrossSection_inline(
		G4CrossSectionDataStore* self,
		const G4DynamicParticle* particle , 
		const G4Material* material ) {
	//By default tries to use the fast-path mechanism
	return G4CrossSectionDataStore_GetCrossSection_withbool(self, particle , material , false);
}

inline void G4CrossSectionDataStore::AddDataSet(G4VCrossSectionDataSet* p)
{
	dataSetList.push_back(p);
	++nDataSetList;
}

inline void G4CrossSectionDataStore::SetVerboseLevel(G4int value)
{
	verboseLevel = value;
}
*/
