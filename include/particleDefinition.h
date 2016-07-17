//source/particles/management/include/G4ParticleDefinition.hh
#ifndef PARTICLEDEFINITION
#define PARTICLEDEFINITION 1

#include "debug.h"
#include "type.h"
#include "physicalConstant.h"

// Class Description:
//
// This class contains all the static data of a particle.
// It also has uses a process manager in order to collect
// all the processes this kind of particle can undertake.

typedef struct particleDefinition {
//protected:

    enum {NumberOfQuarkFlavor = 6} nQuarkFlavor;
    int  theQuarkContent[NumberOfQuarkFlavor];
    int  theAntiQuarkContent[NumberOfQuarkFlavor];
    //  the number of quark (minus Sign means anti-quark) contents
    //  The value of flavor is assigned as follows
    //    0:d, 1:u, 2:s, 3:c, 4:b, 5:t

//private:

    //  --- Following values can not be changed
    //  --- i.e. No Setxxxx Methods for them

    char* theParticleName;
    //  The name of the particle.
    //  Each object must have its specific name!!

    //  --- Following member values must be defined with Units

    dREAL thePDGMass;
    //  The mass of the particle, in units of equivalent energy.

    dREAL thePDGWidth;
    //  The decay width of the particle, usually the width of a
    //  Breit-Wigner function, assuming that you are near the
    //  mass center anyway. (in units of equivalent energy)

    dREAL thePDGCharge;
    //  The charge of the particle.(in units of Coulomb)

    //   --- Following members are quantum number
    //       i.e. discrete numbers can be allowded
    //       So, you can defined only by using integer in constructor

    int thePDGiSpin;
    //  The total spin of the particle, also often denoted as
    //  capital J, in units of 1/2.
    dREAL thePDGSpin;
    //  The total spin of the particle, in units of 1.

    int thePDGiParity;
    //  The parity quantum number, in units of 1. If the parity
    //  is not defined for this particle, we will set this to 0.

    int thePDGiConjugation;
    //  This charge conjugation quantum number in units of 1.

    int thePDGiGParity;
    //  The value of the G-parity quantum number.

    int thePDGiIsospin;
    int thePDGiIsospin3;
    //  The isospin and its 3rd-component in units of 1/2.
    dREAL thePDGIsospin;
    dREAL thePDGIsospin3;
    //  The isospin quantum number in units of 1.

    dREAL thePDGMagneticMoment;
    //  The magnetic moment.

    int theLeptonNumber;
    //  The lepton quantum number.

    int theBaryonNumber;
    //  The baryon quantum number.

    char* theParticleType;
    //  More general textual type description of the particle.

    char* theParticleSubType;
    // Textual type description of the particle
    // eg. pion, lamda etc.

    int thePDGEncoding;
    //  The Particle Data Group integer identifier of this particle

    int theAntiPDGEncoding;
    //  The Particle Data Group integer identifier of the anti-particle

    // --- Following members can be changed after construction

    char fShortLivedFlag; //G4bool fShortLivedFlag;
    //  Particles which have true value of this flag
    //  will not be tracked by TrackingManager

    dREAL thePDGStable;
    //  Is an indicator that this particle is stable. It must
    //  not decay. If the user tries to assign a kind of decay
    //  object to it, it will refuse to take it.

    dREAL thePDGLifeTime;
    //  Is related to the decay width of the particle. The mean
    //  life time is given in seconds.

//    G4DecayTable *theDecayTable;
    //  Points DecayTable

// private:

//    G4ParticleTable* theParticleTable;

    int theAtomicNumber;
    int theAtomicMass;

    int verboseLevel;
    char fApplyCutsFlag; //G4bool fApplyCutsFlag;

// protected:
    char isGeneralIon; //G4bool isGeneralIon;

    int g4particleDefinitionInstanceID;
      // This field is used as instance ID.
} particleDefinition;


// Only one type of constructor can be used for G4ParticleDefinition.
// If you want to create new particle, you must set name of the particle
// at construction. Most of members seen as arguments of the constructor
// (except last 3 arguments concerning with decay ) are  "constant"
// and can not be changed later. (No "SET" methods are available)
// Each type of particle must be constructed as a unique object
// of special class derived from G4ParticleDefinition.
// see G4ParticleTypes for detail

extern void constructParticleDefinition(particleDefinition* def,
									const char*  	aName,
									dREAL         	mass,
									dREAL         	width,
									dREAL        	charge,
									int            iSpin,
									int            iParity,
									int            iConjugation,
									int            iIsospin,
									int            iIsospinZ,
									int            gParity,
									const char*  	pType,
									int            lepton,
									int            baryon,
									int            encoding,
									char           stable, //G4bool           stable,
									dREAL         lifetime,
//									G4DecayTable     *decaytable,
									char           shortlived, //G4bool           shortlived = false,
									const char*  	subType,
									int            anti_encoding,// =0,
									dREAL         magneticMoment //= 0.0
									);

     // With the following Getxxxx methods, one can get values
      // for members which can not be changed

extern const char* particleDefinition_GetParticleName(particleDefinition* def);

extern dREAL particleDefinition_GetPDGMass(particleDefinition* def);
extern dREAL particleDefinition_GetPDGWidth(particleDefinition* def);
extern dREAL particleDefinition_GetPDGCharge(particleDefinition* def);

extern dREAL particleDefinition_GetPDGSpin(particleDefinition* def);
extern int    particleDefinition_GetPDGiSpin(particleDefinition* def);
extern int    particleDefinition_GetPDGiParity(particleDefinition* def);
extern int    particleDefinition_GetPDGiConjugation(particleDefinition* def);
extern dREAL particleDefinition_GetPDGIsospin(particleDefinition* def);
extern dREAL particleDefinition_GetPDGIsospin3(particleDefinition* def);
extern int    particleDefinition_GetPDGiIsospin(particleDefinition* def);
extern int    particleDefinition_GetPDGiIsospin3(particleDefinition* def);
extern int    particleDefinition_GetPDGiGParity(particleDefinition* def);

extern dREAL particleDefinition_GetPDGMagneticMoment(particleDefinition* def);
//extern void  particleDefinition_SetPDGMagneticMoment(particleDefinition def, dREAL mageticMoment);
extern dREAL particleDefinition_CalculateAnomaly(particleDefinition*def);
        // Gives the anomaly of magnetic moment for spin 1/2 particles

extern const char* particleDefinition_GetParticleType(particleDefinition* def);
extern const char* particleDefinition_GetParticleSubType(particleDefinition* def);
extern int    particleDefinition_GetLeptonNumber(particleDefinition* def);
extern int    particleDefinition_GetBaryonNumber(particleDefinition* def);

extern int    particleDefinition_GetPDGEncoding(particleDefinition* def);
extern int    particleDefinition_GetAntiPDGEncoding(particleDefinition* def);
//extern void     particleDefinition_SetAntiPDGEncoding(particleDefinition def, int aEncoding);


//extern int    particleDefinition_GetQuarkContent(particleDefinition def, int flavor);
//extern int    particleDefinition_GetAntiQuarkContent(particleDefinition def, int flavor);
        // Returns the number of quark with flavor contained in this particle.
        // The value of flavor is assigned as follows
        // 1:d, 2:u, 3:s, 4:c, 5:b, 6:t

extern char   particleDefinition_IsShortLived(particleDefinition* def);

extern char   particleDefinition_GetPDGStable(particleDefinition* def);
extern void     particleDefinition_SetPDGStable(particleDefinition* def, const char aFlag);

extern dREAL particleDefinition_GetPDGLifeTime(particleDefinition* def);
extern void     particleDefinition_SetPDGLifeTime(particleDefinition* def, dREAL aLifeTime);

//extern dREAL particleDefinition_GetIonLifeTime(particleDefinition def);
        // Get life time of a generic ion through G4NuclideTable.

//extern G4DecayTable* GetDecayTable();
//extern void          SetDecayTable(G4DecayTable* aDecayTable);
        // Set/Get Decay Table
        //   !! Decay Table can be modified !!

//extern G4ProcessManager* GetProcessManager();
//extern void SetProcessManager(G4ProcessManager* aProcessManager);
        // Set/Get Process Manager
        //   !! Process Manager can be modified !!

//extern G4ParticleTable* GetParticleTable();
        // Get pointer to the particle table

//extern int particleDefinition_GetAtomicNumber(particleDefinition def);
//extern int particleDefinition_GetAtomicMass(particleDefinition def);
        // Get AtomicNumber and AtomicMass
        // These properties are defined for nucleus

extern void particleDefinition_DumpTable(particleDefinition* def);
        //  Prints information of data members.

extern void  particleDefinition_SetVerboseLevel(particleDefinition* def, int value);
extern int particleDefinition_GetVerboseLevel(particleDefinition* def);
        // controle flag for output message
        //  0: Silent
        //  1: Warning message
        //  2: More

extern void   particleDefinition_SetApplyCutsFlag(particleDefinition* def, char);//G4bool);
extern char particleDefinition_GetApplyCutsFlag(particleDefinition* def);

extern char particleDefinition_IsGeneralIon(particleDefinition* def);
      // true only if the particle is G4Ions
      // (it means that theProcessManager is same as one for G4GenricIon)

//      G4int operator==(const G4ParticleDefinition &right) const;
//      G4int operator!=(const G4ParticleDefinition &right) const;

//  public :  // without description

//      inline G4ProcessManager* GetMasterProcessManager() const;
      // Returns the process manager master pointer.
//      inline void SetMasterProcessManager(G4ProcessManager* aNewPM);
      //Sets the shadow master pointer (not to be used by user)

extern int particleDefinition_GetInstanceID(particleDefinition* def);
      // Returns the instance ID.

//      static const G4PDefManager& GetSubInstanceManager();
      // Returns the private data instance manager.
// private:
      // --- Shadow of master pointers.

//      G4ProcessManager *theProcessManagerShadow;
      //  Each worker thread can access this field from the master thread
      //  through this pointer.



//      G4PART_DLL static G4PDefManager subInstanceManager;
      // This field helps to use the class G4PDefManager introduced above.

//  protected:

extern int particleDefinition_FillQuarkContents(particleDefinition* def);
        // Calculates quark and anti-quark contents
        // return value is PDG encoding for this particle.
        // It means error if the return value is deffernt from
        // this->thePDGEncoding.

//extern void particleDefinition_SetParticleSubType(particleDefinition def, const char* subtype);

//extern void particleDefinition_SetAtomicNumber(particleDefinition def, int );
//extern void particleDefinition_SetAtomicMass(particleDefinition def, int );

      //  !!!  can not use "copy constructor" nor "default constructor" !!!!
      //
//      G4ParticleDefinition(const G4ParticleDefinition &right);
//      G4ParticleDefinition();
#endif /* particleDefinition */
