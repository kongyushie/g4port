#include "particleDefinition.h"

void constructParticleDefinition(particleDefinition* def,
								const char*    	 aName,
								dREAL            mass,
								dREAL            width,
								dREAL            charge,
								int               iSpin,
								int               iParity,
								int               iConjugation,
								int               iIsospin,
								int               iIsospin3,
								int               gParity,
								const char*     pType,
								int               lepton,
								int               baryon,
								int               encoding,
								char              stable,
								dREAL            lifetime,
//								G4DecayTable        *decaytable,
								char              shortlived,
								const char*     subType,
								int               anti_encoding,
								dREAL            magneticMoment)
{
	def->theParticleName = aName;
	def->thePDGMass = mass;
	def->thePDGWidth = width;
	def->thePDGCharge = charge;
	def->thePDGiSpin = iSpin;
	def->thePDGSpin = iSpin*0.5;
	def->thePDGiParity = iParity;
	def->thePDGiConjugation = iConjugation;
	def->thePDGiGParity = gParity;
	def->thePDGiIsospin = iIsospin;
	def->thePDGiIsospin3 = iIsospin3;
	def->thePDGIsospin = iIsospin*0.5;
	def->thePDGIsospin3 = iIsospin3*0.5;
	def->thePDGMagneticMoment = magneticMoment;
	def->theLeptonNumber = lepton;
	def->theBaryonNumber = baryon;
	def->theParticleType = pType;
	def->theParticleSubType = subType;
	def->thePDGEncoding = encoding;
	def->theAntiPDGEncoding = -1*encoding;
	def->fShortLivedFlag = shortlived;
	def->thePDGStable = stable;
	def->thePDGLifeTime = lifetime;
	//def->theDecayTable(decaytable),
	def->theAtomicNumber = 0;
	def->theAtomicMass = 0;
	def->verboseLevel = 1;
	def->fApplyCutsFlag = 0; //false;
	def->isGeneralIon = 0; //false;

   static char* nucleus = "nucleus";

   def->g4particleDefinitionInstanceID = -1;
//   def.theProcessManagerShadow = 0;

//   def.theParticleTable = G4ParticleTable::GetParticleTable();

   //set verboseLevel equal to ParticleTable
//   def.verboseLevel = theParticleTable->GetVerboseLevel();

   if (anti_encoding !=0) def->theAntiPDGEncoding = anti_encoding;

   // check quark contents
   if (particleDefinition_FillQuarkContents(def) != def->thePDGEncoding) {
#ifdef G4VERBOSE
     if (verboseLevel>0) {
       // Using G4cout expecting that it is available in construction of static objects
//       G4cout << "Particle " << aName << " has a strange PDGEncoding " <<G4endl;
    	 printf("Particle %s has a strange PDGEncoding\n", aName);
     }
#endif
	BLURT;
	printf("Strange PDGEncoding \n");
	exit(0);
//     G4Exception( "G4ParticleDefintion::G4ParticleDefintion",
//		  "PART102", JustWarning,
//		  "Strange PDGEncoding ");
   }

   // check initialization is in Pre_Init state except for ions
//   G4ApplicationState currentState = G4StateManager::GetStateManager()->GetCurrentState();

	BLURT;
	printf("not implemented: (G4State_PreInit) \n");
	exit(0);

//   if ( !def.fShortLivedFlag && (def.theParticleType!=nucleus) && (def.currentState!=G4State_PreInit)){
//#ifdef G4VERBOSE
//     if (GetVerboseLevel()>0) {
//       G4cout << "G4ParticleDefintion (other than ions and shortlived) should be created in Pre_Init state  "
//              << aName << G4endl;
//     }
//#endif
//	BLURT;
//	printf("G4ParticleDefinition should be created in PreInit state \n");
//	exit(0);
////     G4Exception( "G4ParticleDefintion::G4ParticleDefintion",
////		  "PART101", JustWarning,
////		  "G4ParticleDefinition should be created in PreInit state");
//   }


//   if (theParticleTable->GetIonTable()->IsIon(this)) {
//     SetAtomicNumber( G4int(GetPDGCharge()/eplus) );
//     SetAtomicMass( GetBaryonNumber() );
//   }

//   if (theParticleTable->GetIonTable()->IsAntiIon(this)) {
//     SetAtomicNumber( std::abs(G4int(GetPDGCharge()/eplus)) );
//     SetAtomicMass( std::abs(GetBaryonNumber()) );
//   }

   // check name and register this particle into ParticleTable
//   theParticleTable->Insert(this);

}

//G4ParticleDefinition::G4ParticleDefinition(const G4ParticleDefinition &)
//{
//  G4Exception("G4ParticleDefinition::G4ParticleDefinition()",
//	      "PART001", FatalException,
//	      "Illegal call of copy Constructor for G4ParticleDefinition ");
//}
//
//G4ParticleDefinition::G4ParticleDefinition()
//{
//  G4Exception("G4ParticleDefinition::G4ParticleDefinition()",
//	      "PART001", FatalException,
//	      "Illegal call of default Constructor for G4ParticleDefinition ");
//}

const char* particleDefinition_GetParticleName(particleDefinition* def) {
	return def->theParticleName;
}

dREAL particleDefinition_GetPDGMass(particleDefinition* def){
	return def->thePDGMass;
}

dREAL particleDefinition_GetPDGWidth(particleDefinition* def) {
	return def->thePDGWidth;
}

dREAL particleDefinition_GetPDGCharge(particleDefinition* def) {
	return def->thePDGCharge;
}

dREAL particleDefinition_GetPDGSpin(particleDefinition* def) {
	return def->thePDGSpin;
}

int    particleDefinition_GetPDGiSpin(particleDefinition* def) {
	return def->thePDGiSpin;
}

int    particleDefinition_GetPDGiParity(particleDefinition* def) {
	return def->thePDGiParity;
}

int    particleDefinition_GetPDGiConjugation(particleDefinition* def) {
	return def->thePDGiConjugation;
}

dREAL particleDefinition_GetPDGIsospin(particleDefinition* def) {
	return def->thePDGIsospin;
}

dREAL particleDefinition_GetPDGIsospin3(particleDefinition* def) {
	return def->thePDGIsospin3;
}

int    particleDefinition_GetPDGiIsospin(particleDefinition* def) {
	return def->thePDGiIsospin;
}

int    particleDefinition_GetPDGiIsospin3(particleDefinition* def) {
	return def->thePDGiIsospin3;
}

int    particleDefinition_GetPDGiGParity(particleDefinition* def) {
	return def->thePDGiGParity;
}

dREAL particleDefinition_GetPDGMagneticMoment(particleDefinition* def) {
	return def->thePDGMagneticMoment;
}

//void  particleDefinition_SetPDGMagneticMoment(particleDefinition def, dREAL mageticMoment){
//
//}

dREAL particleDefinition_CalculateAnomaly(particleDefinition* def){
        // Gives the anomaly of magnetic moment for spin 1/2 particles
	BLURT;
    printf("CalculateAnomaly() method will be removed in next release \n");
	exit(0);

//	G4Exception( "G4ParticleDefintion::G4ParticleDefintion",
//	               "PART114", JustWarning,
//	               "CalculateAnomaly() method will be removed in next release");
//	  // gives the anomaly of magnetic moment for spin 1/2 particles
//	  if (thePDGiSpin==1) {
//	    G4double muB = 0.5*CLHEP::eplus*CLHEP::hbar_Planck/(thePDGMass/CLHEP::c_squared);
//	    return 0.5*std::fabs(thePDGMagneticMoment/muB - 2.*thePDGCharge/CLHEP::eplus);
//	  } else {
//	    return 0.0;
//	  }
}

const char* particleDefinition_GetParticleType(particleDefinition* def) {
	return def->theParticleType;
}

const char* particleDefinition_GetParticleSubType(particleDefinition* def) {
	return def->theParticleSubType;
}

int    particleDefinition_GetLeptonNumber(particleDefinition* def) {
	return def->theLeptonNumber;
}

int    particleDefinition_GetBaryonNumber(particleDefinition* def) {
	return def->theBaryonNumber;
}

int    particleDefinition_GetPDGEncoding(particleDefinition* def) {
	return def->thePDGEncoding;
}

int    particleDefinition_GetAntiPDGEncoding(particleDefinition* def) {
	return def->theAntiPDGEncoding;
}

//void     particleDefinition_SetAntiPDGEncoding(particleDefinition def, int aEncoding){
//
//}


//int    particleDefinition_GetQuarkContent(particleDefinition def, int flavor){
//
//}
//
//int    particleDefinition_GetAntiQuarkContent(particleDefinition def, int flavor){
//
//}

char   particleDefinition_IsShortLived(particleDefinition* def) {
	return def->fShortLivedFlag;
}

char   particleDefinition_GetPDGStable(particleDefinition* def){
	return def->thePDGStable;
}

void     particleDefinition_SetPDGStable(particleDefinition* def, const char aFlag) {
	def->thePDGStable=aFlag;
}

dREAL particleDefinition_GetPDGLifeTime(particleDefinition* def){
	return def->thePDGLifeTime;
}

void     particleDefinition_SetPDGLifeTime(particleDefinition* def, dREAL aLifeTime) {
	def->thePDGLifeTime=aLifeTime;
}

//dREAL particleDefinition_GetIonLifeTime(particleDefinition def){
//        // Get life time of a generic ion through G4NuclideTable.
//}

//int particleDefinition_GetAtomicNumber(particleDefinition def){
//
//}
//
//int particleDefinition_GetAtomicMass(particleDefinition def){
//        // Get AtomicNumber and AtomicMass
//        // These properties are defined for nucleus
//}

void particleDefinition_DumpTable(particleDefinition* def){
        //  Prints information of data members.
	printf("\n");
	printf("--- G4ParticleDefinition ---\n");
	printf("Particle Name : %s \n", def->theParticleName);
	printf("PDG particle code : %d \n", def->thePDGEncoding);
	printf("[PDG anti-particle code: %d] \n", particleDefinition_GetAntiPDGEncoding(def));
	printf("Mass [GeV/c2] : %f \n", def->thePDGMass/GeV);
	printf("Width : %f \n", def->thePDGWidth/GeV);
	printf("Lifetime [nsec] : %f \n", def->thePDGLifeTime/ns);
	printf("Charge [e]: %f \n", def->thePDGCharge/eplus);
	printf("Spin : %d /2 \n", def->thePDGiSpin);
	printf("Parity : %d \n", def->thePDGiParity );
	printf("Charge conjugation : %d \n", def->thePDGiConjugation);
	printf("Isospin : (I,Iz): ( %d /2, %d /2 ) \n", def->thePDGiIsospin, def->thePDGiIsospin3);
	printf("GParity : %d \n", def->thePDGiGParity);


	if (def->thePDGMagneticMoment != 0.0) {
		printf("MagneticMoment [MeV/T] : %f \n", def->thePDGMagneticMoment/MeV*tesla);
	}
	printf("Quark contents     (d,u,s,c,b,t) : %d, %d, %d, %d, %d, %d \n",
			def->theQuarkContent[0], def->theQuarkContent[1], def->theQuarkContent[2],
			def->theQuarkContent[3], def->theQuarkContent[4], def->theQuarkContent[5]);
	printf("AntiQuark contents      : %d, %d, %d, %d, %d, %d \n",
				def->theAntiQuarkContent[0], def->theAntiQuarkContent[1], def->theAntiQuarkContent[2],
				def->theAntiQuarkContent[3], def->theAntiQuarkContent[4], def->theAntiQuarkContent[5]);
	printf("Lepton number : %d \n", def->theLeptonNumber);
	printf("Baryon number : %d \n", def->theBaryonNumber);
	printf("Particle type : %s [ %s ]\n", def->theParticleType, def->theParticleSubType);

//	  if (   (theParticleTable->GetIonTable()->IsIon(this))
//	      || (theParticleTable->GetIonTable()->IsAntiIon(this)) ) {
//	    G4cout << " Atomic Number : " << GetAtomicNumber();
//	    G4cout << "  Atomic Mass : " << GetAtomicMass()  << G4endl;
//	  }

	  if ( def->fShortLivedFlag ){
		  printf("ShortLived : ON\n");
	  }

	  if ( particleDefinition_IsGeneralIon(def) ) {
		BLURT;
		printf("unhandled situation: ( GetIonLifeTime() ) \n");
		exit(0);

//	    dREAL lftm = GetIonLifeTime();
//	    if(lftm<-1000.)
//	    { G4cout << " Stable : No data found -- unknown" << G4endl; }
//	    else if(lftm<0.)
//	    { G4cout << " Stable : stable" << G4endl; }
//	    else
//	    {
//	      G4cout << " Stable : unstable -- lifetime = " << G4BestUnit(lftm,"Time")
//	             << "\n  Decay table should be consulted to G4RadioactiveDecayProcess."
//	             << G4endl;
//	    }
	  }
	  else{
	    if ( def->thePDGStable ){
	    	printf("Stable : stable\n");
	    } else {
	      if( 0 /*def.theDecayTable != 0*/ ){
//	        theDecayTable->DumpInfo();
	      } else {
	    	  printf("Decay Table is not defined !!\n");
	      }
	    }
	  }
}

void  particleDefinition_SetVerboseLevel(particleDefinition* def, int value){
	def->verboseLevel = value;
}

int particleDefinition_GetVerboseLevel(particleDefinition* def){
        // controle flag for output message
        //  0: Silent
        //  1: Warning message
        //  2: More
}

void   particleDefinition_SetApplyCutsFlag(particleDefinition* def, char flg){
	if(def->theParticleName=="gamma"
	  || def->theParticleName=="e-"
	  || def->theParticleName=="e+"
	  || def->theParticleName=="proton")
	  {
		def->fApplyCutsFlag = flg;
	  }
	  else{
		printf("ParticleDefinition::SetApplyCutsFlag() for %s\n", def->theParticleName);
		printf("becomes obsolete. Production threshold is applied only for gamma, e- ,e+ and proton.\n");
	  }
}

char particleDefinition_GetApplyCutsFlag(particleDefinition* def){
	return def->fApplyCutsFlag;
}

char particleDefinition_IsGeneralIon(particleDefinition* def){
	return def->isGeneralIon;
}

//int particleDefinition_GetInstanceID(particleDefinition def){
//      // Returns the instance ID.
//}

int particleDefinition_FillQuarkContents(particleDefinition* def){
        // Calculates quark and anti-quark contents
        // return value is PDG encoding for this particle.
        // It means error if the return value is deffernt from
        // this->thePDGEncoding.
	  int flavor;
	  for (flavor= 0; flavor<NumberOfQuarkFlavor; flavor++){
		  def->theQuarkContent[flavor]     = 0;
		  def->theAntiQuarkContent[flavor] = 0;
	  }

	  int temp = 0;

	  BLURT;
	  printf("unhandled situation: (G4PDGCodeChecker checker) \n");
	  exit(0);


//	  G4PDGCodeChecker checker;
//	  checker.SetVerboseLevel(verboseLevel);
//
//	  int temp = checker.CheckPDGCode(thePDGEncoding, theParticleType);
//
//	  if ( temp != 0) {
//	    for (flavor= 0; flavor<NumberOfQuarkFlavor; flavor++){
//	      theQuarkContent[flavor]     = checker.GetQuarkContent(flavor);
//	      theAntiQuarkContent[flavor] = checker.GetAntiQuarkContent(flavor);
//	    }
//	    if ((theParticleType == "meson")||(theParticleType == "baryon")) {
//	      // check charge
//	      if (!checker.CheckCharge(thePDGCharge) ){
//		temp = 0;
//		G4Exception( "G4ParticleDefintion::G4ParticleDefintion",
//			  "PART103", JustWarning,
//			  "Inconsistent charge against PDG code ");
//	#ifdef G4VERBOSE
//		if (verboseLevel>0) {
//		  G4cout << "G4ParticleDefinition::FillQuarkContents  : "
//		         << " illegal charge (" << thePDGCharge/eplus
//		         << " PDG code=" << thePDGEncoding <<G4endl;
//		}
//	#endif
//	      }
//	      // check spin
//	      if (checker.GetSpin() != thePDGiSpin) {
//		temp=0;
//		G4Exception( "G4ParticleDefintion::G4ParticleDefintion",
//			  "PART104", JustWarning,
//			  "Inconsistent spin against PDG code ");
//	#ifdef G4VERBOSE
//		if (verboseLevel>0) {
//		  G4cout << "G4ParticleDefinition::FillQuarkContents  : "
//		         << " illegal SPIN (" << thePDGiSpin << "/2"
//		         << " PDG code=" << thePDGEncoding <<G4endl;
//		}
//	#endif
//	      }
//	    }
//	  }
	  return temp;


}

//void particleDefinition_SetParticleSubType(particleDefinition def, const char* subtype){
//
//}

//void particleDefinition_SetAtomicNumber(particleDefinition def, int ){
//
//}
//
//void particleDefinition_SetAtomicMass(particleDefinition def, int ){
//
//}


//-- No longer needed to access to G4IonTable.
//-- Method GetIonLifeTime() itself is kept for compatibility
//-- but moved to icc file as an inlined method.
//G4double G4ParticleDefinition::GetIonLifeTime() const
//{
//  if(!isGeneralIon) return thePDGLifeTime;
//
//  G4IonTable* ionTable =  G4IonTable::GetIonTable();
//  return ionTable->GetLifeTime(this);
//}


//void G4ParticleDefinition::SetParticleDefinitionID(G4int id)
//{
//  if(id<0)
//  {
//    g4particleDefinitionInstanceID = subInstanceManager.CreateSubInstance();
//    G4MT_pmanager = 0;
//  }
//  else
//  {
//    if(isGeneralIon)
//    { g4particleDefinitionInstanceID = id; }
//    else
//    {
//      G4ExceptionDescription ed;
//      ed << "ParticleDefinitionID should not be set for the particles <"
//         << theParticleName << ">.";
//      G4Exception( "G4ParticleDefintion::SetParticleDefinitionID","PART10114",
//                   FatalException,ed);
//    }
//  }
//}


