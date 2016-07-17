//source/global/management/include/G4ThreeVector.hh (typedef CLHEP::Hep3Vector G4ThreeVector;)
//source/externals/CLHEP/Vector/ThreeVector.h
#ifndef THREEVECTOR
#define THREEVECTOR
#include "debug.h"
#include "type.h"

enum {  threeVector_X=0,
		threeVector_Y=1,
		threeVector_Z=2,
		threeVector_NUM_COORDINATES=3,
		threeVector_SIZE=threeVector_NUM_COORDINATES };
  // Safe indexing of the coordinates when using with matrices, arrays, etc.
  // (BaBar)


typedef struct threeVector{
	dREAL dx;
	dREAL dy;
	dREAL dz;
	// The components.

	/*DLL_API static*/ dREAL tolerance;
	// default tolerance criterion for isNear() to return true.

} threeVector;
/*
	void constructthreeVector(threeVector *self);
  //explicit threeVector(double x);
  void constructthreeVector_xy(threeVector *self, double x, double y);
  void constructthreeVector_xyz(threeVector *self, double x, double y, double z);
  // The constructor.

  //inline threeVector(const threeVector &);
  // The copy constructor.

  inline void destructthreeVector(threeVector *self);
  // The destructor.  Not virtual - inheritance from this class is dangerous.

  inline double operator () (int) const;
  // Get components by index -- 0-based (Geant4)

  inline double & operator () (int);
  // Set components by index.  0-based.

  inline double threeVector_phi(const threeVector *self);
  // The azimuth angle.

  inline double threeVector_theta(const threeVector *self);
  // The polar angle.

  inline double threeVector_cosTheta(const threeVector *self);
  // Cosine of the polar angle.

  inline double threeVector_cos2Theta(const threeVector *self);
  // Cosine squared of the polar angle - faster than cosTheta(). (ZOOM)

  inline double threeVector_mag2(const threeVector *self);
  // The magnitude squared (r^2 in spherical coordinate system).

  inline double threeVector_mag(const threeVector *self);
  // The magnitude (r in spherical coordinate system).

  inline void threeVector_setPhi(threeVector *self, double);
  // Set phi keeping mag and theta constant (BaBar).

  inline void threeVector_setTheta(threeVector *self, double);
  // Set theta keeping mag and phi constant (BaBar).

         void threeVector_setMag(threeVector *self, double);
  // Set magnitude keeping theta and phi constant (BaBar).

  inline double threeVector_perp2(const threeVector *self);
  // The transverse component squared (rho^2 in cylindrical coordinate system).

  inline double threeVector_perp(const threeVector *self);
  // The transverse component (rho in cylindrical coordinate system).

  inline void threeVector_setPerp(threeVector *self, double);
  // Set the transverse component keeping phi and z constant.

  void threeVector_setCylTheta(threeVector *self, double);
  // Set theta while keeping transvers component and phi fixed

  inline double threeVector_perp2(const threeVector *, const threeVector *);
  // The transverse component w.r.t. given axis squared.

  inline double threeVector_perp(const threeVector *, const threeVector *);
  // The transverse component w.r.t. given axis.

  // inline threeVector & operator = (const threeVector &);
  // Assignment.

  inline bool threeVector_OPeq (const threeVector *);
  inline bool threeVector_OPne (const threeVector *);
  // Comparisons (Geant4).

  bool isNear (const threeVector *, const threeVector *, double epsilon=tolerance);
  // Check for equality within RELATIVE tolerance (default 2.2E-14). (ZOOM)
  // |v1 - v2|**2 <= epsilon**2 * |v1.dot(v2)|

  double howNear(const threeVector *, const threeVector & v );
  // sqrt ( |v1-v2|**2 / v1.dot(v2) ) with a maximum of 1.
  // If v1.dot(v2) is negative, will return 1.

  double deltaR( const threeVector *, const threeVector * v);
  // sqrt( pseudorapity_difference**2 + deltaPhi **2 )

  inline threeVector & operator += (const threeVector &);
  // Addition.

  inline threeVector & operator -= (const threeVector &);
  // Subtraction.

  inline threeVector operator - () const;
  // Unary minus.

  inline threeVector & operator *= (double);
  // Scaling with real numbers.

         threeVector & operator /= (double);
  // Division by (non-zero) real number.

  inline threeVector unit() const;
  // Vector parallel to this, but of length 1.

  inline threeVector orthogonal() const;
  // Vector orthogonal to this (Geant4).

  inline double dot(const threeVector &) const;
  // double product.

  inline threeVector cross(const threeVector &) const;
  // Cross product.

  double angle(const threeVector &) const;
  // The angle w.r.t. another 3-vector.

  double pseudoRapidity() const;
  // Returns the pseudo-rapidity, i.e. -ln(tan(theta/2))

  void setEta  ( double p );
  // Set pseudo-rapidity, keeping magnitude and phi fixed.  (ZOOM)

  void setCylEta  ( double p );
  // Set pseudo-rapidity, keeping transverse component and phi fixed.  (ZOOM)

  threeVector & rotateX(double);
  // Rotates the threeVector around the x-axis.

  threeVector & rotateY(double);
  // Rotates the threeVector around the y-axis.

  threeVector & rotateZ(double);
  // Rotates the threeVector around the z-axis.

  threeVector & rotateUz(const threeVector&);
  // Rotates reference frame from Uz to newUz (unit vector) (Geant4).

    threeVector & rotate(double, const threeVector &);
  // Rotates around the axis specified by another threeVector.
  // (Uses methods of HepRotation, forcing linking in of Rotation.cc.)

  threeVector & operator *= (const HepRotation &);
  threeVector & transform(const HepRotation &);
  // Transformation with a Rotation matrix.

// = = = = = = = = = = = = = = = = = = = = = = = =
//
// Esoteric properties and operations on 3-vectors:
//
// 1 - Set vectors in various coordinate systems
// 2 - Synonyms for accessing coordinates and properties
// 3 - Comparisions (dictionary, near-ness, and geometric)
// 4 - Intrinsic properties
// 5 - Properties releative to z axis and arbitrary directions
// 6 - Polar and azimuthal angle decomposition and deltaPhi
// 7 - Rotations
//
// = = = = = = = = = = = = = = = = = = = = = = = =

// 1 - Set vectors in various coordinate systems

  inline void setRThetaPhi  (double r, double theta, double phi);
  // Set in spherical coordinates:  Angles are measured in RADIANS

  inline void setREtaPhi  ( double r, double eta,  double phi );
  // Set in spherical coordinates, but specify peudorapidiy to determine theta.

  inline void setRhoPhiZ   (double rho, double phi, double z);
  // Set in cylindrical coordinates:  Phi angle is measured in RADIANS

  void setRhoPhiTheta ( double rho, double phi, double theta);
  // Set in cylindrical coordinates, but specify theta to determine z.

  void setRhoPhiEta ( double rho, double phi, double eta);
  // Set in cylindrical coordinates, but specify pseudorapidity to determine z.

// 2 - Synonyms for accessing coordinates and properties

  inline double getX() const;
  inline double getY() const;
  inline double getZ() const;
  // x(), y(), and z()

  inline double getR    () const;
  inline double getTheta() const;
  inline double getPhi  () const;
  // mag(), theta(), and phi()

  inline double r       () const;
  // mag()

  inline double rho     () const;
  inline double getRho  () const;
  // perp()

  double eta     () const;
  double getEta  () const;
  // pseudoRapidity()

  inline void setR ( double s );
  // setMag()

  inline void setRho ( double s );
  // setPerp()

// 3 - Comparisions (dictionary, near-ness, and geometric)

  int compare (const threeVector & v) const;
  bool operator > (const threeVector & v) const;
  bool operator < (const threeVector & v) const;
  bool operator>= (const threeVector & v) const;
  bool operator<= (const threeVector & v) const;
  // dictionary ordering according to z, then y, then x component

  inline double diff2 (const threeVector & v) const;
  // |v1-v2|**2

  static double setTolerance (double tol);
  static inline double getTolerance ();
  // Set the tolerance used in isNear() for threeVectors

  bool isParallel (const threeVector & v, double epsilon=tolerance) const;
  // Are the vectors parallel, within the given tolerance?

  bool isOrthogonal (const threeVector & v, double epsilon=tolerance) const;
  // Are the vectors orthogonal, within the given tolerance?

  double howParallel   (const threeVector & v) const;
  // | v1.cross(v2) / v1.dot(v2) |, to a maximum of 1.

  double howOrthogonal (const threeVector & v) const;
  // | v1.dot(v2) / v1.cross(v2) |, to a maximum of 1.

  enum { ToleranceTicks = 100 };

// 4 - Intrinsic properties

  double beta    () const;
  // relativistic beta (considering v as a velocity vector with c=1)
  // Same as mag() but will object if >= 1

  double gamma() const;
  // relativistic gamma (considering v as a velocity vector with c=1)

  double coLinearRapidity() const;
  // inverse tanh (beta)

// 5 - Properties relative to Z axis and to an arbitrary direction

	  // Note that the non-esoteric CLHEP provides
	  // theta(), cosTheta(), cos2Theta, and angle(const threeVector&)

  inline double angle() const;
  // angle against the Z axis -- synonym for theta()

  inline double theta(const threeVector & v2) const;
  // synonym for angle(v2)

  double cosTheta (const threeVector & v2) const;
  double cos2Theta(const threeVector & v2) const;
  // cos and cos^2 of the angle between two vectors

  inline threeVector project () const;
         threeVector project (const threeVector & v2) const;
  // projection of a vector along a direction.

  inline threeVector perpPart() const;
  inline threeVector perpPart (const threeVector & v2) const;
  // vector minus its projection along a direction.

  double rapidity () const;
  // inverse tanh(v.z())

  double rapidity (const threeVector & v2) const;
  // rapidity with respect to specified direction:
  // inverse tanh (v.dot(u)) where u is a unit in the direction of v2

  double eta(const threeVector & v2) const;
  // - ln tan of the angle beween the vector and the ref direction.

// 6 - Polar and azimuthal angle decomposition and deltaPhi

  // Decomposition of an angle within reference defined by a direction:

  double polarAngle (const threeVector & v2) const;
  // The reference direction is Z: the polarAngle is abs(v.theta()-v2.theta()).

  double deltaPhi (const threeVector & v2) const;
  // v.phi()-v2.phi(), brought into the range (-PI,PI]

  double azimAngle  (const threeVector & v2) const;
  // The reference direction is Z: the azimAngle is the same as deltaPhi

  double polarAngle (const threeVector & v2,
					const threeVector & ref) const;
  // For arbitrary reference direction,
  // 	polarAngle is abs(v.angle(ref) - v2.angle(ref)).

  double azimAngle  (const threeVector & v2,
					const threeVector & ref) const;
  // To compute azimangle, project v and v2 into the plane normal to
  // the reference direction.  Then in that plane take the angle going
  // clockwise around the direction from projection of v to that of v2.

// 7 - Rotations

// These mehtods **DO NOT** use anything in the HepRotation class.
// Thus, use of v.rotate(axis,delta) does not force linking in Rotation.cc.

  threeVector & rotate  (const threeVector & axis, double delta);
  // Synonym for rotate (delta, axis)

  threeVector & rotate  (const HepAxisAngle & ax);
  // HepAxisAngle is a struct holding an axis direction and an angle.

  threeVector & rotate (const HepEulerAngles & e);
  threeVector & rotate (double phi,
                        double theta,
                        double psi);
  // Rotate via Euler Angles. Our Euler Angles conventions are
  // those of Goldstein Classical Mechanics page 107.

//protected:
  void setSpherical (double r, double theta, double phi);
  void setCylindrical (double r, double phi, double z);
  double negativeInfinity() const;
*/
#endif /* threeVector */
