#ifndef FastPathHadronicCrossSection_hh
#define FastPathHadronicCrossSection_hh
#include "testglobal.h"
//#include "G4PhysicsFreeVector.hh"
//#include "G4ParticleDefinition.hh"
//#include "G4Material.hh"
//#include <functional>
//#include <utility>
//#include <unordered_map>
//#include <iostream>
//#include <set>
//#include <stdint.h>

//class G4DynamicParticle;
//class G4Material;
//class G4CrossSectionDataStore;

//namespace G4FastPathHadronicCrossSection {
#define XSParam G4PhysicsFreeVector
//The key used to search in the cache.
#define G4CrossSectionDataStore_Key pair<const G4ParticleDefinition*,const G4Material*>
//This represents the fast XS implementation.
typedef struct G4FastPathHadronicCrossSection_fastPathEntry{
	const G4ParticleDefinition * const particle;
	const G4Material * const material;
	const G4double min_cutoff;

	XSParam *physicsVector;
#       ifdef FPDEBUG
	//stats for debug
	G4int count;
	G4double slowpath_sum; //sum of all slowpath xs
	G4double max_delta;
	G4double min_delta;
	G4double sum_delta;
	G4double sum_delta_square;
#		endif
}G4FastPathHadronicCrossSection_fastPathEntry;
extern void constructG4FastPathHadronicCrossSection_fastPathEntry(
		G4FastPathHadronicCrossSection_fastPathEntry* self, 
		const G4ParticleDefinition *par,
		const G4Material* mat,
		G4double min_cutoff);
extern void deconstructG4FastPathHadronicCrossSection_fastPathEntry(G4FastPathHadronicCrossSection_fastPathEntry* self);
inline G4double G4FastPathHadronicCrossSection_fastPathEntry_GetCrossSection(
		G4FastPathHadronicCrossSection_fastPathEntry* self,
		G4double ene) const { return self.physicsVector->Value(ene); }
extern void G4FastPathHadronicCrossSection_fastPathEntry_Initialize(
		G4FastPathHadronicCrossSection_fastPathEntry* self,
		G4CrossSectionDataStore* );

//A cache entry.
typedef struct G4FastPathHadronicCrossSection_cycleCountEntry{
	//const G4String& particle;
	const G4String* particle;
	const G4Material * const material;

	//optional fastPathEntry
	G4FastPathHadronicCrossSection_fastPathEntry* fastPath;

	//cache per element of material test
	G4double energy;
	G4double crossSection;
#	  ifdef FPDEBUG
	uint64_t cacheHitCount;//
	uint64_t initCyclesFastPath;
	uint64_t invocationCountSlowPath;
	uint64_t totalCyclesSlowPath;
	uint64_t invocationCountFastPath;
	uint64_t totalCyclesFastPath;
	uint64_t invocationCountTriedOneLineCache;//
	uint64_t invocationCountOneLineCache;//
#	  endif
}G4FastPathHadronicCrossSection_cycleCountEntry;
extern void constructG4FastPathHadronicCrossSection_cycleCountEntry(
		G4FastPathHadronicCrossSection_cycleCountEntry* self,
		const G4String& pname , 
		const G4Material* mat);

extern void destructcycleCountEntry(G4FastPathHadronicCrossSection_cycleCountEntry* self);

typedef struct G4FastPathHadronicCrossSection_timing {
	unsigned long long rdtsc_start;
	unsigned long long rdtsc_stop;
}G4FastPathHadronicCrossSection_timing;

typedef struct G4FastPathHadronicCrossSection_getCrossSectionCount {
#ifdef FPDEBUG
	uint64_t methodCalled;
	uint64_t hitOneLineCache;
	uint64_t fastPath;
	uint64_t slowPath;
	uint64_t sampleZandA;
#endif
}G4FastPathHadronicCrossSection_getCrossSectionCount;
constructG4FastPathHadronicCrossSection_getCrossSectionCount(G4FastPathHadronicCrossSection_getCrossSectionCount* self);
inline void G4FastPathHadronicCrossSection_getCrossSectionCount_MethodCalled(
		G4FastPathHadronicCrossSection_getCrossSectionCount* self);
inline void G4FastPathHadronicCrossSection_getCrossSectionCount_HitOneLine(
		G4FastPathHadronicCrossSection_getCrossSectionCount* self);
inline void G4FastPathHadronicCrossSection_getCrossSectionCount_FastPath(
		G4FastPathHadronicCrossSection_getCrossSectionCount* self);
inline void G4FastPathHadronicCrossSection_getCrossSectionCount_SlowPath(
		G4FastPathHadronicCrossSection_getCrossSectionCount* self);
inline void G4FastPathHadronicCrossSection_getCrossSectionCount_SampleZandA(
		G4FastPathHadronicCrossSection_getCrossSectionCount* self);

//Hashing the key
typedef struct G4FastPathHadronicCrossSection_G4CrossSectionDataStore_Key_Hash {
	hash<uint64_t> hash_uint64_t;
}G4FastPathHadronicCrossSection_G4CrossSectionDataStore_Key_Hash;
//inline size_t operator()(const G4CrossSectionDataStore_Key& x) const throw() {
//	return hash_uint64_t(hash_uint64_t( ((uint64_t)(x.first)) ) +  hash_uint64_t(((uint64_t)(x.second))));
//	}
inline size_t OP_G4FastPathHadronicCrossSection_G4CrossSectionDataStore_Key_Hash(
		G4FastPathHadronicCrossSection_G4CrossSectionDataStore_Key_Hash* self,
		const G4CrossSectionDataStore_Key* x){
	return hash_uint64_t(self.hash_uint64_t( ((uint64_t)(x.first)) )+self.hash_uint64_t(((uint64_t)(x.second))));
}
//Equality for two key elements
typedef struct G4FastPathHadronicCrossSection_G4CrossSectionDataStore_Key_EqualTo {
}G4FastPathHadronicCrossSection_G4CrossSectionDataStore_Key_EqualTo;
inline bool OP_G4FastPathHadronicCrossSection_G4CrossSectionDataStore_Key_EqualTo(
//inline bool operator()(const G4CrossSectionDataStore_Key& lhs, const G4CrossSectionDataStore_Key& rhs ) const {
		const G4CrossSectionDataStore_Key* lhs, 
		const G4CrossSectionDataStore_Key* rhs)
{
		//TODO: Verify this: particles are singletons, materials use operator==
		//TODO: in ref-10, G4Material::operator== becomes deleted, investigating why
	return (lhs.first==rhs.first)&&(lhs.second == rhs.second);
}
//	The cache itself
using G4CrossSectionDataStore_Cache=std::unordered_map<G4CrossSectionDataStore_Key,cycleCountEntry*,G4CrossSectionDataStore_Key_Hash,G4CrossSectionDataStore_Key_EqualTo>;

typedef struct G4FastPathHadronicCrossSection_fastPathRequestConfig_t {
	G4CrossSectionDataStore_Key part_mat;
	G4double min_cutoff;
}G4FastPathHadronicCrossSection_fastPathRequestConfig_t;

//Two of the elements are identical if the part_mat part is
typedef struct G4FastPathHadronicCrossSection_fastPathRequestConfig_Less {
	//less<G4CrossSectionDataStore_Key> less;
}G4FastPathHadronicCrossSection_fastPathRequestConfig_Less;
/*inline bool operator()(const fastPathRequestConfig_t& lhs,const fastPathRequestConfig_t& rhs ) const {
	return less(lhs.part_mat,rhs.part_mat);
}*/

using G4CrossSectionDataStore_Requests=std::set<fastPathRequestConfig_t,fastPathRequestConfig_Less>;

//Configure the caching mechanism
typedef struct G4CrossSectionDataStore_controlFlag{
	G4bool prevCalcUsedFastPath;
	G4bool useFastPathIfAvailable;
	G4bool initializationPhase;
}G4CrossSectionDataStore_controlFlag;
inline void constructG4CrossSectionDataStore_controlFlag(G4CrossSectionDataStore_controlFlag* self)
{
	self.prevCalcUsedFastPath = false;
	self.useFastPathIfAvailable = false;
	self.initializationPhase = false;
}
//Parameters to control sampling
struct G4CrossSectionDataStore_fastPathParameters {
	//PRUTH vars for sampling and surragate model
	G4double queryMax;
	G4double sampleMin;
	G4double sampleMax;
	G4int sampleCount;
	G4double dpTol;
}G4CrossSectionDataStore_fastPathParameters;
inline void constructG4CrossSectionDataStore_fastPathParameters(G4CrossSectionDataStore_fastPathParameters* self) {
		//default
		//TODO: are these ok?
		self.queryMax = 10000;
		self.sampleMin = 0.0001;
		self.sampleMax = 10000;
		self.sampleCount = 200000;
		self.dpTol = 0.01;
	}

//Logging functionalities, disabled if not in FPDEBUG mode
/*static inline void logInvocationTriedOneLine( cycleCountEntry* );
static inline void logInvocationOneLine( cycleCountEntry* );
static inline void logHit(cycleCountEntry*);
static inline void logInvocationCountFastPath( cycleCountEntry* );
static inline void logInvocationCountSlowPAth( cycleCountEntry* );
#ifdef FPDEBUG
void logStartCountCycles( timing& );
void logStopCountCycles( timing& );
#else
inline void logStartCountCycles(timing&) {}
inline void logStopCountCycles(timing&) {}
#endif
static inline void logInitCyclesFastPath( cycleCountEntry* , timing& );
static inline void logTotalCyclesFastPath( cycleCountEntry* , timing& );
static inline void logTotalCyclesSlowPath( cycleCountEntry* , timing& );
static inline void logTiming( cycleCountEntry* , fastPathEntry* , timing& );
//}

inline std::ostream& operator<<(std::ostream& os, const G4FastPathHadronicCrossSection::fastPathEntry& fp);

//Implementation of inline functions. Note the ifdef

namespace G4FastPathHadronicCrossSection {

#ifdef FPDEBUG
	inline void logInvocationTriedOneLine(cycleCountEntry* cl ) {
		if ( cl != nullptr ) ++(cl->invocationCountTriedOneLineCache);
	}
	inline void logInvocationOneLine( cycleCountEntry* cl ) {
		if ( cl != nullptr ) ++(cl->invocationCountOneLineCache);
	}
	inline void logHit(cycleCountEntry* cl) {
		if ( cl != nullptr ) ++(cl->cacheHitCount);
	}
	inline void logInvocationCountFastPath( cycleCountEntry* cl )
	{
		if ( cl != nullptr ) ++(cl->invocationCountFastPath);
	}
	inline void logInvocationCountSlowPAth( cycleCountEntry* cl)
	{
		if ( cl != nullptr ) ++(cl->invocationCountSlowPath);
	}

	inline void logInitCyclesFastPath(cycleCountEntry* cl,timing& tm)
	{
		if ( cl != nullptr ) cl->initCyclesFastPath = tm.rdtsc_stop - tm.rdtsc_start;
	}
	inline void logTotalCyclesFastPath( cycleCountEntry* cl,timing& tm)
	{
		if ( cl!=nullptr ) cl->totalCyclesFastPath = tm.rdtsc_stop - tm.rdtsc_start;
	}
	inline void logTotalCyclesSlowPath( cycleCountEntry* cl,timing& tm)
	{
		if ( cl!=nullptr ) cl->totalCyclesSlowPath = tm.rdtsc_stop - tm.rdtsc_start;
	}
	inline void logTiming( cycleCountEntry* entry , fastPathEntry* fast_entry, timing& timing)
	{
		if (fast_entry != nullptr ) {
			if ( entry->invocationCountFastPath == 0 ) {
				//PRUTH style initialization
				G4FastPathHadronicCrossSection::logInitCyclesFastPath(entry,timing);
				G4FastPathHadronicCrossSection::logInvocationCountFastPath(entry);
			} else {
				//PRUTH comment to understand:
				//the first one includes the initialization... don't count it for now
				G4FastPathHadronicCrossSection::logTotalCyclesFastPath(entry,timing);
				G4FastPathHadronicCrossSection::logInvocationCountFastPath(entry);
			}
		} else {
			G4FastPathHadronicCrossSection::logInvocationCountSlowPAth(entry);
			G4FastPathHadronicCrossSection::logTotalCyclesSlowPath(entry,timing);
		}
	}
#else
	inline void logInvocationTriedOneLine(cycleCountEntry*){}
	inline void logInvocationOneLine( cycleCountEntry*){}
	inline void logHit(cycleCountEntry*){}
	inline void logInvocationCountFastPath( cycleCountEntry*){}
	inline void logInvocationCountSlowPAth( cycleCountEntry*){}
	inline void logInitCyclesFastPath( cycleCountEntry* , timing& ){}
	inline void logTotalCyclesFastPath( cycleCountEntry* , timing& ){}
	inline void logTotalCyclesSlowPath( cycleCountEntry* , timing& ){}
	inline void logTiming( cycleCountEntry* , fastPathEntry* , timing& ) {}
#endif

	inline void G4FastPathHadronicCrossSection_getCrossSectionCount_MethodCalled(
			G4FastPathHadronicCrossSection_getCrossSectionCount* self) {
#ifdef FPDEBUG
		++(self.methodCalled);
#endif
	}

	inline void G4FastPathHadronicCrossSection_getCrossSectionCount_HitOneLine(
			G4FastPathHadronicCrossSection_getCrossSectionCount* self) {
#ifdef FPDEBUG
		++(self.hitOneLineCache);
#endif
	}

	inline void G4FastPathHadronicCrossSection_getCrossSectionCount_FastPath(G4FastPathHadronicCrossSection_getCrossSectionCount* self) {
#ifdef FPDEBUG
		++(self.fastPath);
#endif
	}

	inline void G4FastPathHadronicCrossSection_getCrossSectionCount_SlowPath(G4FastPathHadronicCrossSe    ction_getCrossSectionCount* self) {
#ifdef FPDEBUG
		++(self.slowPath);
#endif
	}

	inline void G4FastPathHadronicCrossSection_getCrossSectionCount_SampleZandA(
			G4FastPathHadronicCrossSection_getCrossSectionCount* self) {
#ifdef FPDEBUG
		++(self.sampleZandA);
#endif
	}
}//namespace

inline std::ostream& operator<<(std::ostream& os, const G4FastPathHadronicCrossSection::fastPathEntry& fp) {
	using CLHEP::MeV;
	os<<"#Particle: "<<(fp.particle!=nullptr?fp.particle->GetParticleName():"UNDEFINED")<<"\n";
	os<<"#Material: "<<(fp.material!=nullptr?fp.material->GetName():"UNDEFINED")<<"\n";
	os<<"#min_cutoff(MeV): "<<fp.min_cutoff/MeV<<"\n";
#ifdef FPDEBUG
	os<<"#DEBUG COUNTERS: count="<<fp.count<<" slowpath_sum="<<fp.slowpath_sum<<" max_delta="<<fp.max_delta;
	os<<" min_delta="<<fp.min_delta<<" sum_delta="<<fp.sum_delta<<" sum_delta_square="<<fp.sum_delta_square<<"\n";
#endif
	os<<*(fp.physicsVector)<<"\n";
	return os;
}*/

#endif //G4FastPathHadronicCrossSection_hh
