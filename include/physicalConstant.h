//source/externals/clhep/include/CLHEP/Units/PhysicalConstants.h

#include "systemUnits.h"

//
//
#define Avogadro  (6.02214179e+23/mole)

//
// c   = 299.792458 mm/ns
// c^2 = 898.7404 (mm/ns)^2
//
#define c_light    (2.99792458e+8 * m/s)
#define c_squared  (c_light * c_light)

//
// h     = 4.13566e-12 MeV*ns
// hbar  = 6.58212e-13 MeV*ns
// hbarc = 197.32705e-12 MeV*mm
//
#define h_Planck       (6.62606896e-34 * joule*s)
#define hbar_Planck    (h_Planck/twopi)
#define hbarc          (hbar_Planck * c_light)
#define hbarc_squared  (hbarc * hbarc)

//
//
//
#define electron_charge  (- eplus) // see SystemOfUnits.
#define e_squared  (eplus * eplus)

//
// amu_c2 - atomic equivalent mass unit
//        - AKA, unified atomic mass unit (u)
// amu    - atomic mass unit
//
#define electron_mass_c2  (0.510998910 * MeV)
#define   proton_mass_c2  (938.272013 * MeV)
#define  neutron_mass_c2  (939.56536 * MeV)
#define           amu_c2  (931.494028 * MeV)
#define              amu  (amu_c2/c_squared)

//
// permeability of free space mu0    = 2.01334e-16 Mev*(ns*eplus)^2/mm
// permittivity of free space epsil0 = 5.52636e+10 eplus^2/(MeV*mm)
//
#define mu0       (4*pi*1.e-7 * henry/m)
#define epsilon0  (1./(c_squared*mu0))

//
// electromagnetic coupling = 1.43996e-12 MeV*mm/(eplus^2)
//
#define elm_coupling            (e_squared/(4*pi*epsilon0))
#define fine_structure_const    (elm_coupling/hbarc)
#define classic_electr_radius   (elm_coupling/electron_mass_c2)
#define electron_Compton_length  (hbarc/electron_mass_c2)
#define Bohr_radius  (electron_Compton_length/fine_structure_const)

#define alpha_rcl2  (fine_structure_cons \
                                   *classic_electr_radius \
                                   *classic_electr_radius)

#define twopi_mc2_rcl2  (twopi*electron_mass_c \
                                             *classic_electr_radius \
                                             *classic_electr_radius)
//
//
//
#define k_Boltzmann  (8.617343e-11 * MeV/kelvin)

//
//
//
#define STP_Temperature  (273.15*kelvin)
#define STP_Pressure     (1.*atmosphere)
#define kGasThreshold    (10.*mg/cm3)

//
//
//
#define universe_mean_density  (1.e-25*g/cm3)
