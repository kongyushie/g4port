
#include "VCrossSectionDataSet.h"
//#include "G4SystemOfUnits.hh"
//#include "G4CrossSectionDataSetRegistry.hh"
#include "dynamicParticle.h"
#include "material.h"
#include "element.h"
#include "elementVector.h"
#include "isotope.h"
#include "isotopeVector.h"
#include <stdbool.h>
//#include "G4NistManager.hh"
//#include "G4HadronicException.hh"
//#include "G4HadTmpUtil.hh"
//#include "Randomize.hh"

//#include <string.h>
//G4VCrossSectionDataSet::G4VCrossSectionDataSet(const G4String& nam) :
void constructVCrossSectionDataSet(VCrossSectionDataSet* self, const char* nam)
  //verboseLevel(0),minKinEnergy(0.0),maxKinEnergy(100*TeV),name(nam) 
{
  BLURT;
  printf("with porting CrossSectionDataSetRegistry\n");
  exit(0);
  /*
  self.verboseLevel = 0;
  self.minKinEnergy = 0.0;
  self.maxKinEnergy = 100*TeV;
  //self.name = nam;
  strcpy(self.name,nam);
  self.registry = CrossSectionDataSetRegistry_Instance();
  CrossSectionDataSetRegistry_Register(self.registry, self);
  */
}

//G4VCrossSectionDataSet::~G4VCrossSectionDataSet()
void destructVCrossSectionDataSet(VCrossSectionDataSet* self)
{
  BLURT;
  printf("with porting CrossSectionDataSetRegistry\n");
  exit(0);
  //CrossSectionDataSetRegistry_DeRegister(self.registry, self);
}

bool VCrossSectionDataSet_IsElementApplicable(
		dREAL self,
		const dynamicParticle* dp, 
		int Z,
        const material* mat) 
{
  return false;
}

bool 
VCrossSectionDataSet_IsIsoApplicable(const dynamicParticle* dp, 
                                        int Z, int A,
                                        const element* elm,  
                                        const material* mat)
{
  return false;
}

/*
G4double 
G4VCrossSectionDataSet::ComputeCrossSection(const G4DynamicParticle* part, 
					    const G4Element* elm,
					    const G4Material* mat)
{
  G4int Z = G4lrint(elm->GetZ());

  if (IsElementApplicable(part, Z, mat)) { 
    return GetElementCrossSection(part, Z, mat);
  }

  // isotope-wise cross section making sum over available
  // isotope cross sections, which may be incomplete, so
  // the result is corrected 
  G4int nIso = elm->GetNumberOfIsotopes();    
  G4double fact = 0.0;
  G4double xsec = 0.0;
  G4Isotope* iso = 0;

  if (0 < nIso) { 

    // user-defined isotope abundances        
    G4IsotopeVector* isoVector = elm->GetIsotopeVector();
    G4double* abundVector = elm->GetRelativeAbundanceVector();

    for (G4int j = 0; j<nIso; ++j) {
      iso = (*isoVector)[j];
      G4int A = iso->GetN();
      if(abundVector[j] > 0.0 && IsIsoApplicable(part, Z, A, elm, mat)) {
        fact += abundVector[j];
	xsec += abundVector[j]*GetIsoCrossSection(part, Z, A, iso, elm, mat);
      }
    }

  } else {

    // natural isotope abundances
    G4NistManager* nist = G4NistManager::Instance();
    G4int n0 = nist->GetNistFirstIsotopeN(Z);
    G4int nn = nist->GetNumberOfNistIsotopes(Z);
    for (G4int A = n0; A < n0+nn; ++A) {
      G4double abund = nist->GetIsotopeAbundance(Z, A);
      if(abund > 0.0 && IsIsoApplicable(part, Z, A, elm, mat)) {
        fact += abund;
	xsec += abund*GetIsoCrossSection(part, Z, A, iso, elm, mat);
      }
    }
  }
  if(fact > 0.0) { xsec /= fact; }
  return xsec;
}
*/
dREAL 
VCrossSectionDataSet_GetElementCrossSection(const dynamicParticle* dynPart,
					       int Z,
					       const material* mat)
{
//  G4cout << "G4VCrossSectionDataSet::GetCrossSection per element ERROR: "
//	 << " there is no cross section for "
//	 << dynPart->GetDefinition()->GetParticleName()
//	 << "  E(MeV)= "  << dynPart->GetKineticEnergy()/MeV;
 // if(mat) { G4cout << "  inside " << mat->GetName(); }
 // G4cout << " for Z= " << Z << G4endl;
 // throw G4HadronicException(__FILE__, __LINE__,
 //       "G4VCrossSectionDataSet::GetElementCrossSection is absent");
  return 0.0;
}

dREAL 
VCrossSectionDataSet_GetIsoCrossSection(
					   //VCrossSectionDataStore* self,
					   dREAL self,
					   const dynamicParticle* dynPart,
					   int Z, int A,
					   const isotope* iso,
					   const element* elm,
					   const material* mat)
{
  /*G4cout << "G4VCrossSectionDataSet::GetCrossSection per isotope ERROR: "
	 << " there is no cross section for "
	 << dynPart->GetDefinition()->GetParticleName()
	 << "  E(MeV)= "  << dynPart->GetKineticEnergy()/MeV;
  if(mat) { G4cout << "  inside " << mat->GetName(); }
  if(elm) { G4cout << " for " << elm->GetName(); }
  G4cout << "  Z= " << Z << " A= " << A << G4endl;
  throw G4HadronicException(__FILE__, __LINE__,
        "G4VCrossSectionDataSet::GetIsoCrossSection is absent");
  */
  return 0.0;
}

isotope* 
VCrossSectionDataSet_SelectIsotope(VCrossSectionDataSet* self,
		const element* anElement, 
		dREAL kinEnergy)
{
  int nIso = element_GetNumberOfIsotopes(anElement);
  isotopeVector* isoVector = element_GetIsotopeVector(anElement);
  isotope* iso = &((*isoVector).v[0]);

  // more than 1 isotope
  if(1 < nIso) {
    dREAL* abundVector = element_GetRelativeAbundanceVector(anElement);
    dREAL sum = 0.0;
    //dREAL q = Randomize_UniformRand();//may use Rand()
	dREAL q = (dREAL)(rand()%100)/100;
    for (int j = 0; j<nIso; ++j) {
      sum += abundVector[j];
      if(q <= sum) {
	iso = &((*isoVector).v[j]);
	break;
      }
    }
  }
  return iso;
}
/*
void G4VCrossSectionDataSet::BuildPhysicsTable(const G4ParticleDefinition&)
{}

void G4VCrossSectionDataSet::DumpPhysicsTable(const G4ParticleDefinition&)
{}

void G4VCrossSectionDataSet::CrossSectionDescription(std::ostream& outFile) const
{
  outFile << "The description for this cross section data set has not been written yet.\n";
}
*/
