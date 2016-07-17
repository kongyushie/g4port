//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....
#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
//#include "G4SystemOfUnits.hh"
//#include "G4UnitsTable.hh"
//#include "G4HadronicException.hh"
//#include "G4HadTmpUtil.hh"
//#include "Randomize.hh" in source/global/HEPRandom/include
//#include "nucleus.h"
//#include "material.h"
#include "debug.h"
#include "crossSectionDataStore.h"
#include "elementVector.h"
#include "isotopeVector.h"
#include "VCrossSectionDataSet.h"
#include "dynamicParticle.h"
//#include "G4Isotope.hh"
#include "element.h"
#include "material.h"
//#include "G4NistManager.hh"
//#include <algorithm>


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

/*G4CrossSectionDataStore::G4CrossSectionDataStore() :
  nDataSetList(0), verboseLevel(0),fastPathFlags(),fastPathParams(),
  counters(),fastPathCache()
  {
  nist = G4NistManager::Instance();
  currentMaterial = elmMaterial = 0;
  currentElement = 0;  //ALB 14-Aug-2012 Coverity fix.
  matParticle = elmParticle = 0;
  matKinEnergy = elmKinEnergy = matCrossSection = elmCrossSection = 0.0;
  }*/
void constructG4CrossSectionDataStore(crossSectionDataStore* self){
	(*self).nDataSetList = 0;
	(*self).verboseLevel = 0;
	//self.nist = G4NistManager::Instance();
	(*self).currentMaterial = NULL;
	(*self).elmMaterial = NULL;
	(*self).currentElement = NULL;  //ALB 14-Aug-2012 Coverity fix.
	(*self).matParticle = NULL;
	(*self).elmParticle = NULL;
	(*self).matKinEnergy = 0.0;
	(*self).elmKinEnergy = 0.0;
	(*self).matCrossSection = 0.0;
	(*self).elmCrossSection = 0.0;
}
void cnostructcrossSectionDataStore_crossSectionDataStore(
		crossSectionDataStore* self,
		const crossSectionDataStore* ano)
{
	BLURT;
	printf("without porting\n");
	exit(0);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

//G4CrossSectionDataStore::~G4CrossSectionDataStore(){}
void destructcrossSectionDataStore(crossSectionDataStore* self)
{return;}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....
	dREAL
crossSectionDataStore_GetCrossSection(
		crossSectionDataStore* self,
		const dynamicParticle* part,
		const material* mat)
{
	//if(mat == currentMaterial && part->GetDefinition() == matParticle
	//&& part->GetKineticEnergy() == matKinEnergy)
	if(   OP_EQ_material(mat, (*self).currentMaterial) 
	   && OP_EQ_dynamicParticle(dynamicParticle_GetDefinition(part), (*self).matParticle)
	   && (dynamicParticle_GetKineticEnergy(part) == (*self).matKinEnergy) )
	{ return (*self).matCrossSection; }

	//currentMaterial = mat;
	OP_ASSIGN_material((*self).currentMaterial, mat);
	//matParticle = part->GetDefinition();
	OP_ASSIGN_dynamicParticle((*self).matParticle, dynamicParticle_GetDefinition(part));
	(*self).matKinEnergy = dynamicParticle_GetKineticEnergy(part);
	(*self).matCrossSection = 0;

	//G4int nElements = mat->GetNumberOfElements();
	int nElements;
	nElements = material_GetNumberOfElements(mat);
	printf("d\n");
	printf("%d\n",nElements);
	//const G4double* nAtomsPerVolume = mat->GetVecNbOfAtomsPerVolume();
	dREAL* nAtomsPerVolume;
	nAtomsPerVolume = material_GetVecNbOfAtomsPerVolume(mat);
	printf("d\n");

	//if(G4int(xsecelm.size()) < nElements) { xsecelm.resize(nElements); }
	if((int)(*self).xsecelm.size < nElements){vector_dREAL_resize(&((*self).xsecelm),nElements);}
	for(int i=0; i < nElements; ++i) {
		//matCrossSection += nAtomsPerVolume[i] *
		//	GetCrossSection(part, (*mat->GetElementVector())[i], mat);
		(*self).matCrossSection += nAtomsPerVolume[i] * 
			crossSectionDataStore_GetCrossSection(self, material_GetElementVector(mat, i), mat);
		(*self).xsecelm.v[i] = (*self).matCrossSection;
	}
	return (*self).matCrossSection;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

/*	void
	G4CrossSectionDataStore::DumpFastPath(const G4ParticleDefinition* pd, const G4Material* mat,std::ostream& os)
	{
	const G4FastPathHadronicCrossSection::cycleCountEntry* entry = fastPathCache[{pd,mat}];
	if ( entry != nullptr ) {
	if ( entry->fastPath != nullptr ) {
	os<<*entry->fastPath;
	} else {
	os<<"#Cache entry for {"<<(pd!=nullptr?pd->GetParticleName():"UNDEFINED")<<",";
	os<<(mat!=nullptr?mat->GetName():"UNDEFINED")<<"} found, but no fast path defined";
	}
	} else {
	os<<"#Cache entry for {"<<(pd!=nullptr?pd->GetParticleName():"UNDEFINED")<<",";
	os<<(mat!=nullptr?mat->GetName():"UNDEFINED")<<"} not found.";
	}
	}
 */
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....
/*
   G4double
   G4CrossSectionDataStore_GetCrossSection(
   G4CrossSectionDataStore* self,
   const G4DynamicParticle* part,
   const G4Element* elm,
   const G4Material* mat)
   {
//if(mat == elmMaterial && elm == currentElement &&
//   part->GetDefinition() == elmParticle &&
//   part->GetKineticEnergy() == elmKinEnergy)
if(OP_EQ_G4Material(mat, self.elmMaterial) &&
OP_EQ_G4Element(elm, self.currentElement) &&
G4DynamicParticle_GetKineticEnergy(part) == self.elmKinEnergy) 
{ return self.elmCrossSection; }

//elmMaterial = mat;
OP_ASSIGN_G4Material(self.elmMaterial, mat);
//currentElement = elm;
OP_ASSIGN_G4Element(self.currentElement, elm);
//elmParticle = part->GetDefinition();
OP_ASSIGN_G4ParticleDefinition(self.elmParticle, G4ParticleDefinition_GetDefinition(part));
elmKinEnergy = part->GetKineticEnergy();
elmCrossSection = 0.0;

G4int i = nDataSetList-1;  
G4int Z = G4lrint(G4Element_GetZ(elm));//?
if (G4Element_GetNaturalAbundanceFlag(elm) &&
G4VCrossSectionDataSet_IsElementApplicable(self.dataSetList[i], part, Z, mat)) 
{
// element wise cross section
self.elmCrossSection = G4VCrossSectionDataSet_GetElementCrossSection(self.dataSetList[i], part, Z, mat);
} else {
// isotope wise cross section
G4int nIso = G4Element_GetNumberOfIsotopes(elm);    
G4Isotope* iso = 0;

// user-defined isotope abundances        
G4IsotopeVector* isoVector = G4Element_GetIsotopeVector(elm);
G4double* abundVector = G4Element_GetRelativeAbundanceVector(elm);

for (G4int j = 0; j<nIso; ++j) {
if(abundVector[j] > 0.0) {
iso = (*isoVector)[j];
self.elmCrossSection += abundVector[j]*
G4CrossSectionDataStore_GetIsoCrossSection(self, part, Z, G4Isotope_GetN(iso), iso, elm, mat, i);
}
}
}
return self.elmCrossSection;
}
 */
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

dREAL
crossSectionDataStore_GetIsoCrossSection(
		crossSectionDataStore* self,
		const dynamicParticle* part,
		int Z, int A, 
		const isotope* iso,
		const element* elm,
		const material* mat, 
		int idx)
{
	// this methods is called after the check that dataSetList[idx] 
	// depend on isotopes, so for this DataSet only isotopes are checked

	// isotope-wise cross section does exist
	if(VCrossSectionDataSet_IsIsoApplicable((*self).dataSetList.v[idx],part, Z, A, elm, mat) ) {
		return VCrossSectionDataSet_GetIsoCrossSection((*self).dataSetList.v[idx],part, Z, A, iso, elm, mat);

	} else {
		// seach for other dataSet
		for (int j = (*self).nDataSetList-1; j >= 0; --j) { 
			if (VCrossSectionDataSet_IsElementApplicable((*self).dataSetList.v[j],part, Z, mat)) {
				return VCrossSectionDataSet_GetElementCrossSection((*self).dataSetList.v[j],part, Z, mat);
			} else if (VCrossSectionDataSet_IsIsoApplicable((*self).dataSetList.v[j],part, Z, A, elm, mat)) {
				return VCrossSectionDataSet_GetIsoCrossSection((*self).dataSetList.v[j],part, Z, A, iso, elm, mat);
			}
		}
	}
	//G4cout << "G4CrossSectionDataStore::GetCrossSection ERROR: "
	//<< " no isotope cross section found"
	//<< G4endl;
	//G4cout << "  for " << part->GetDefinition()->GetParticleName() 
	//<< " off Element " << elm->GetName()
	//<< "  in " << mat->GetName() 
	//<< " Z= " << Z << " A= " << A
	//<< " E(MeV)= " << part->GetKineticEnergy()/MeV << G4endl; 
	//throw G4HadronicException(__FILE__, __LINE__, 
	//" no applicable data set found for the isotope");
	return 0.0;
	//return dataSetList[idx]->ComputeCrossSection(part, elm, mat);
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

/*	G4double 
	G4CrossSectionDataStore_GetCrossSection_ZA(
	G4CrossSectionDataStore* self,
	const G4DynamicParticle* part,
	G4int Z, 
	G4int A,
	const G4Isotope* iso,
	const G4Element* elm,
	const G4Material* mat)
	{
	for (G4int i = self.nDataSetList-1; i >= 0; --i) {
	if (G4VCrossSectionDataSet_IsIsoApplicable(self.dataSetList[i], part, Z, A, elm, mat) ) {
	return GetIsoCrossSection(self.dataSetList[i], part, Z, A, iso, elm, mat);
	}
	}
	G4cout << "G4CrossSectionDataStore::GetCrossSection ERROR: "
	<< " no isotope cross section found"
	<< G4endl;
	G4cout << "  for " << part->GetDefinition()->GetParticleName() 
	<< " off Element " << //elm->GetName()// *(G4Element_GetName(elm))
	<< "  in " << //mat->GetName()// *(G4Material_GetName(mat))
	<< " Z= " << Z << " A= " << A
	<< " E(MeV)= " << part->GetKineticEnergy()/MeV << G4endl; 
	throw G4HadronicException(__FILE__, __LINE__, 
	" no applicable data set found for the isotope");
	return 0.0;
	}
 */
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

	element*
crossSectionDataStore_SampleZandA(
		crossSectionDataStore* self,
		const dynamicParticle* part, 
		const material* mat,
		nucleus* target)
{
	//G4FastPathHadronicCrossSection_SampleZandA(self.counters);

	int nElements = material_GetNumberOfElements(mat);
	const elementVector* theElementVector = material_GetElementVector(mat);
	element* anElement = &((*theElementVector).v[0]);
	printf("a\n");

	dREAL cross = crossSectionDataStore_GetCrossSection(self, part, mat);
	printf("a\n");
	// select element from a compound 
	if(1 < nElements) {
		//cross *= Randomize_G4UniformRand();
		dREAL rand_temp = (dREAL)(rand()%100)/100;
		cross *= rand_temp;
		for(int i=0; i<nElements; ++i) {
			if(cross <= (*self).xsecelm.v[i]) {
				anElement = (*theElementVector).v[i];
				break;
			}
		}
	}

	//int Z = G4lrint(element_GetZ(anElement));
	int Z = (element_GetZ(anElement) > 0)? (int)(element_GetZ(anElement)+0.5):(int)(element_GetZ(anElement)-0.5);
	isotope* iso = NULL;

	int i = (*self).nDataSetList-1; 
	if (VCrossSectionDataSet_IsElementApplicable((*self).dataSetList.v[i],part, Z, mat)) {

		//----------------------------------------------------------------
		// element-wise cross section
		// isotope cross section is not computed
		//----------------------------------------------------------------
		int nIso = element_GetNumberOfIsotopes(anElement);
		if (0 >= nIso) { 
			//G4cout << " Element " << anElement->GetName() << " Z= " << Z 
			//<< " has no isotopes " << G4endl;
			printf("Element %s Z= %d has no isotopes\n",element_GetName(anElement),Z); 
			BLURT;
			exit(0);
			//throw G4HadronicException(__FILE__, __LINE__, 
			//" Isotope vector is not defined");
			return anElement;
		}
		// isotope abundances        
		isotopeVector* isoVector = element_GetIsotopeVector(anElement);
		iso = &((*isoVector).v[0]);

		// more than 1 isotope
		if(1 < nIso) { 
			iso = VCrossSectionDataSet_SelectIsotope(
					&((*self).dataSetList.v[i]), 
					anElement, 
					dynamicParticle_GetKineticEnergy(part));
		}	

	} else {

		//----------------------------------------------------------------
		// isotope-wise cross section
		// isotope cross section is computed
		//----------------------------------------------------------------
		int nIso = element_GetNumberOfIsotopes(anElement);
		cross = 0.0;

		if (0 >= nIso) { 
			//G4cout << " Element " << anElement->GetName() << " Z= " << Z 
			//<< " has no isotopes " << G4endl; 
			printf("Element %s Z= %d has no isotopes\n",element_GetName(anElement),Z); 
			BLURT;
			exit(0);
			//throw G4HadronicException(__FILE__, __LINE__, 
			//" Isotope vector is not defined");
			return anElement;
		}		

		// user-defined isotope abundances        
		isotopeVector* isoVector = element_GetIsotopeVector(anElement);
		iso = &((*isoVector).v[0]);

		// more than 1 isotope
		if(1 < nIso) {
			double* abundVector = element_GetRelativeAbundanceVector(anElement);
			if((*self).xseciso.size < nIso) { vector_dREAL_resize(&((*self).xseciso),nIso); }

			for (int j = 0; j<nIso; ++j) {
				double xsec = 0.0;
				if(abundVector[j] > 0.0) {
					iso = &((*isoVector).v[j]);
					xsec = abundVector[j]*
						crossSectionDataStore_GetIsoCrossSection(
								self,
								part, 
								Z, 
								isotope_GetN(iso), 
								iso, 
								anElement, 
								mat, 
								i);
				}	
				cross += xsec;
				(*self).xseciso.v[j] = cross;
			}
			//cross *= G4UniformRand();
			dREAL rand_temp = (dREAL)(rand()%100)/100;
			cross *= rand_temp;
			for (int j = 0; j<nIso; ++j) {
				if(cross <= (*self).xseciso.v[j]) {
					iso = &((*isoVector).v[j]);
					break;
				}
			}
		}
	}
	nucleus_SetIsotope(target,iso);
	return anElement;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

/*	void
	G4CrossSectionDataStore::BuildPhysicsTable(const G4ParticleDefinition& aParticleType)
	{
	if (nDataSetList == 0) 
	{
	throw G4HadronicException(__FILE__, __LINE__, 
	"G4CrossSectionDataStore: no data sets registered");
	return;
	}
	for (G4int i=0; i<nDataSetList; ++i) {
	dataSetList[i]->BuildPhysicsTable(aParticleType);
	} 
//A.Dotti: if fast-path has been requested we can now create the surrogate
//         model for fast path.
if ( fastPathFlags.useFastPathIfAvailable ) {
fastPathFlags.initializationPhase = true;
using my_value_type=G4FastPathHadronicCrossSection::G4CrossSectionDataStore_Requests::value_type;
//Loop on all requests, if particle matches create the corresponding fsat-path
std::for_each( requests.begin() , requests.end() ,
[&aParticleType,this](const my_value_type& req) {
if ( aParticleType == *req.part_mat.first ) {
G4FastPathHadronicCrossSection::cycleCountEntry* entry =
new G4FastPathHadronicCrossSection::cycleCountEntry(aParticleType.GetParticleName(),req.part_mat.second);
entry->fastPath =
new G4FastPathHadronicCrossSection::fastPathEntry(&aParticleType,req.part_mat.second,req.min_cutoff);
entry->fastPath->Initialize(this);
fastPathCache[req.part_mat] = entry;
}
}
);
fastPathFlags.initializationPhase = false;
}
}
 */
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

/*void G4CrossSectionDataStore::ActivateFastPath( const G4ParticleDefinition* pdef, const G4Material* mat, G4double min_cutoff)
  {
  assert(pdef!=nullptr&&mat!=nullptr);
  G4FastPathHadronicCrossSection::G4CrossSectionDataStore_Key key={pdef,mat};
  if ( requests.insert( { key , min_cutoff } ).second ) {
  std::ostringstream msg;
  msg<<"Attempting to request FastPath for couple: "<<pdef->GetParticleName()<<","<<mat->GetName();
  msg<<" but combination already exists";
  G4HadronicException(__FILE__,__LINE__,msg.str());
  }
  }
 */
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

/*void 
  G4CrossSectionDataStore::DumpPhysicsTable(const G4ParticleDefinition& aParticleType)
  {
// Print out all cross section data sets used and the energies at
// which they apply

if (nDataSetList == 0) {
G4cout << "WARNING - G4CrossSectionDataStore::DumpPhysicsTable: "
<< " no data sets registered" << G4endl;
return;
}

for (G4int i = nDataSetList-1; i >= 0; --i) {
G4double e1 = dataSetList[i]->GetMinKinEnergy();
G4double e2 = dataSetList[i]->GetMaxKinEnergy();
G4cout 
<< "     Cr_sctns: " << std::setw(25) << dataSetList[i]->GetName() << ": "
<<  G4BestUnit(e1, "Energy")
<< " ---> "
<<  G4BestUnit(e2, "Energy") << "\n";
if (dataSetList[i]->GetName() == "G4CrossSectionPairGG") {
dataSetList[i]->DumpPhysicsTable(aParticleType);
}
}
}
 */
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....
//#include <typeinfo>
//void G4CrossSectionDataStore::DumpHtml(const G4ParticleDefinition& /* pD */,
//                                       std::ofstream& outFile) const
/*{
// Write cross section data set info to html physics list
// documentation page

G4double ehi = 0;
G4double elo = 0;
G4String physListName(getenv("G4PhysListName"));
for (G4int i = nDataSetList-1; i > 0; i--) {
elo = dataSetList[i]->GetMinKinEnergy()/GeV;
ehi = dataSetList[i]->GetMaxKinEnergy()/GeV;
outFile << "      <li><b><a href=\"" << physListName << "_"
<< dataSetList[i]->GetName() << ".html\"> "
<< dataSetList[i]->GetName() << "</a> from "
<< elo << " GeV to " << ehi << " GeV </b></li>\n";
//G4cerr << i << ": XS for " << pD.GetParticleName() << " : " << dataSetList[i]->GetName() 
//       << " typeid : " << typeid(dataSetList[i]).name()<< G4endl;			
PrintCrossSectionHtml(dataSetList[i]);			
}

G4double defaultHi = dataSetList[0]->GetMaxKinEnergy()/GeV;
if (ehi < defaultHi) {
outFile << "      <li><b><a href=\"" << dataSetList[0]->GetName() << ".html\"> "
<< dataSetList[0]->GetName() << "</a> from "
<< ehi << " GeV to " << defaultHi << " GeV </b></li>\n";
PrintCrossSectionHtml(dataSetList[0]);			
}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void G4CrossSectionDataStore::PrintCrossSectionHtml(const G4VCrossSectionDataSet *cs) const
{
G4String dirName(getenv("G4PhysListDocDir"));
G4String physListName(getenv("G4PhysListName"));

G4String pathName = dirName + "/" + physListName + "_" + HtmlFileName(cs->GetName());
std::ofstream outCS;
outCS.open(pathName);
outCS << "<html>\n";
outCS << "<head>\n";
outCS << "<title>Description of " << cs->GetName() 
<< "</title>\n";
outCS << "</head>\n";
outCS << "<body>\n";

cs->CrossSectionDescription(outCS);

outCS << "</body>\n";
outCS << "</html>\n";

}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

G4String G4CrossSectionDataStore::HtmlFileName(const G4String & in) const
{
G4String str(in);
// replace blanks by _  C++11 version:
#ifdef G4USE_STD11
std::transform(str.begin(), str.end(), str.begin(), [](char ch) {
return ch == ' ' ? '_' : ch;
});
#else	
// and now in ancient language
for(std::string::iterator it = str.begin(); it != str.end(); ++it) {
if(*it == ' ') *it = '_';
}
#endif
str=str + ".html";		
return str;
}

*/
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....
