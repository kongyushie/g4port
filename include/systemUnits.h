//source/externals/clhep/include/CLHEP/Units/SystemOfUnits.h
#pragma once
  static const double     pi  = 3.14159265358979323846;
  static const double  twopi  = 2*3.14159265358979323846;
  static const double halfpi  = 3.14159265358979323846/2;
  static const double     pi2 = 3.14159265358979323846*3.14159265358979323846;

  // 
  // Length [L]
  //
  static const double millimeter  = 1.;                        
  static const double millimeter2 = 1.*1.;
  static const double millimeter3 = 1.*1.*1.;

  static const double centimeter  = 10.*1.;   
  static const double centimeter2 = 10.*10.;
  static const double centimeter3 = 10.*10.*10.;
    
  static const double meter  = 1000.*1.;                  
  static const double meter2 = 1000.*1000.;
  static const double meter3 = 1000.*1000.*1000.;

  static const double kilometer = 1000.*1000.;                   
  static const double kilometer2 = 1000000.*1000000.;
  static const double kilometer3 = 1000000.*1000000.*1000000.;

  static const double parsec = 3.0856775807e+16*(1000.*1.);

  static const double micrometer = 1.e-6 *1000.*1.;             
  static const double  nanometer = 1.e-9 *1000.*1.;
  static const double  angstrom  = 1.e-10*1000.*1.;
  static const double  fermi     = 1.e-15*1000.*1.;

  static const double      barn = 1.e-28*1000.*1000.;
  static const double millibarn = 1.e-3 *(1.e-28*1000.*1000.);
  static const double microbarn = 1.e-6 *(1.e-28*1000.*1000.);
  static const double  nanobarn = 1.e-9 *(1.e-28*1000.*1000.);
  static const double  picobarn = 1.e-12*(1.e-28*1000.*1000.);

  // symbols
  static const double nm  = 1.e-9 *1000.*1.;                        
  static const double um  = 1.e-6 *1000.*1.;                        

  static const double mm  = 1.;                        
  static const double mm2 = 1.*1.;
  static const double mm3 = 1.*1.*1.;

  static const double cm  = 10.*1.;   
  static const double cm2 = 10.*10.;
  static const double cm3 = 10.*10.*10.;

  static const double liter = 1.e+3*(10.*10.*10.);
  static const double  L = (1.e+3*(10.*10.*10.));
  static const double dL = 1.e-1*(1.e+3*(10.*10.*10.));
  static const double cL = 1.e-2*(1.e+3*(10.*10.*10.));
  static const double mL = 1.e-3*(1.e+3*(10.*10.*10.));       
  
  static const double m  = 1000.*1.;                  
  static const double m2 = 1000.*1000.;
  static const double m3 = 1000.*1000.*1000.;

  static const double km  = 1000.*1000.;                   
  static const double km2 = 1000000.*1000000.;
  static const double km3 = 1000000.*1000000.*1000000.;

  static const double pc = 3.0856775807e+16*(1000.*1.);

  //
  // Angle
  //
  static const double radian      = 1.;                  
  static const double milliradian = 1.e-3*1.;
  static const double degree = (3.14159265358979323846/180.0)*1.;

  static const double   steradian = 1.;
  
  // symbols
  static const double rad  = 1.;
  static const double mrad = (1.e-3*1);
  static const double sr   = 1.;
  static const double deg  = (3.14159265358979323846/180.0)*1.;

  //
  // Time [T]
  //
  static const double nanosecond  = 1.;
  static const double second      = 1.e+9 *1.;
  static const double millisecond = 1.e-3 *(1.e+9 *1.);
  static const double microsecond = 1.e-6 *(1.e+9 *1.);
  static const double  picosecond = 1.e-12*(1.e+9 *1.);

  static const double hertz = 1./(1.e+9 *1.);
  static const double kilohertz = 1.e+3*(1./(1.e+9 *1.));
  static const double megahertz = 1.e+6*(1./(1.e+9 *1.));

  // symbols
  static const double ns = 1.;
  static const double  s = (1.e+9 *1.);
  static const double ms = (1.e-3 *(1.e+9 *1.));

  //
  // Electric charge [Q]
  //
  static const double eplus = 1. ;// positron charge
  static const double e_SI  = 1.602176487e-19;// positron charge in coulomb
  static const double coulomb = (1.)/(1.602176487e-19);// coulomb = 6.24150 e+18 * eplus

  //
  // Energy [E]
  //
  static const double megaelectronvolt = 1. ;
  static const double     electronvolt = 1.e-6*1.;
  static const double kiloelectronvolt = 1.e-3*1.;
  static const double gigaelectronvolt = 1.e+3*1.;
  static const double teraelectronvolt = 1.e+6*1.;
  static const double petaelectronvolt = 1.e+9*1.;

  static const double joule = (1.e-6*1)/(1.602176487e-19);// joule = 6.24150 e+12 * MeV

  // symbols
  static const double MeV = 1.;
  static const double  eV = 1.e-6*1.;
  static const double keV = 1.e-3*1.;
  static const double GeV = 1.e+3*1.;
  static const double TeV = 1.e+6*1.;
  static const double PeV = 1.e+9*1.;

  //
  // Mass [E][T^2][L^-2]
  //
  static const double  kilogram = ((1.e-6*1)/(1.602176487e-19))*(1.e+9 *1.)*(1.e+9 *1.)/((1000.*1.)*(1000.*1.));   
  static const double      gram = 1.e-3*(((1.e-6*1)/(1.602176487e-19))*(1.e+9 *1.)*(1.e+9 *1.)/((1000.*1.)*(1000.*1.)));
  static const double milligram = 1.e-3*(1.e-3*(((1.e-6*1)/(1.602176487e-19))*(1.e+9 *1.)*(1.e+9 *1.)/((1000.*1.)*(1000.*1.))));

  // symbols
  static const double  kg = ((1.e-6*1)/(1.602176487e-19))*(1.e+9 *1.)*(1.e+9 *1.)/((1000.*1.)*(1000.*1.));
  static const double   g = 1.e-3*(((1.e-6*1)/(1.602176487e-19))*(1.e+9 *1.)*(1.e+9 *1.)/((1000.*1.)*(1000.*1.)));
  static const double  mg = 1.e-3*(1.e-3*(((1.e-6*1)/(1.602176487e-19))*(1.e+9 *1.)*(1.e+9 *1.)/((1000.*1.)*(1000.*1.))));

  //
  // Power [E][T^-1]
  //
  static const double watt = ((1.e-6*1)/(1.602176487e-19))/((1.e+9 *1.));// watt = 6.24150 e+3 * MeV/ns

  //
  // Force [E][L^-1]
  //
  static const double newton = ((1.e-6*1)/(1.602176487e-19))/(1000.*1.);// newton = 6.24150 e+9 * MeV/mm

  //
  // Pressure [E][L^-3]
  //
#define pascal hep_pascal                          // a trick to avoid warnings 
  static const double hep_pascal = (((1.e-6*1)/(1.602176487e-19))/(1000.*1.))/(1000.*1000.);   // pascal = 6.24150 e+3 * MeV/mm3
  static const double bar        = 100000*((((1.e-6*1)/(1.602176487e-19))/(1000.*1.))/(1000.*1000.)); // bar    = 6.24150 e+8 * MeV/mm3
  static const double atmosphere = 101325*((((1.e-6*1)/(1.602176487e-19))/(1000.*1.))/(1000.*1000.)); // atm    = 6.32420 e+8 * MeV/mm3

  //
  // Electric current [Q][T^-1]
  //
  static const double      ampere = (1.)/(1.602176487e-19)/(1.e+9 *1.); // ampere = 6.24150 e+9 * eplus/ns
  static const double milliampere = 1.e-3*((1.)/(1.602176487e-19)/(1.e+9 *1.));
  static const double microampere = 1.e-6*((1.)/(1.602176487e-19)/(1.e+9 *1.));
  static const double  nanoampere = 1.e-9*((1.)/(1.602176487e-19)/(1.e+9 *1.));

  //
  // Electric potential [E][Q^-1]
  //
  static const double megavolt = 1./1.;
  static const double kilovolt = 1.e-3*(1./1.);
  static const double     volt = 1.e-6*(1./1.);

  //
  // Electric resistance [E][T][Q^-2]
  //
  static const double ohm = (1.e-6*(1./1.))/((1.)/(1.602176487e-19)/(1.e+9 *1.));// ohm = 1.60217e-16*(MeV/eplus)/(eplus/ns)

  //
  // Electric capacitance [Q^2][E^-1]
  //
  static const double farad = ((1.)/(1.602176487e-19))/(1.e-6*(1./1.));// farad = 6.24150e+24 * eplus/Megavolt
  static const double millifarad = 1.e-3*(((1.)/(1.602176487e-19))/(1.e-6*(1./1.)));
  static const double microfarad = 1.e-6*(((1.)/(1.602176487e-19))/(1.e-6*(1./1.)));
  static const double  nanofarad = 1.e-9*(((1.)/(1.602176487e-19))/(1.e-6*(1./1.)));
  static const double  picofarad = 1.e-12*(((1.)/(1.602176487e-19))/(1.e-6*(1./1.)));

  //
  // Magnetic Flux [T][E][Q^-1]
  //
  static const double weber = (1.e-6*(1./1.))*(1.e+9 *1.);// weber = 1000*megavolt*ns

  //
  // Magnetic Field [T][E][Q^-1][L^-2]
  //
  static const double tesla     = (1.e-6*(1./1.))*(1.e+9 *1.)/(1000.*1000.);// tesla =0.001*megavolt*ns/mm2
/*
  static const double gauss     = 1.e-4*tesla;
  static const double kilogauss = 1.e-1*tesla;

  //
  // Inductance [T^2][E][Q^-2]
  //
  static const double henry = weber/((1.)/(1.602176487e-19)/(1.e+9 *1.));// henry = 1.60217e-7*MeV*(ns/eplus)**2

  //
  // Temperature
  //
  static const double kelvin = 1.;

  //
  // Amount of substance
  //
  static const double mole = 1.;

  //
  // Activity [T^-1]
  //
  static const double becquerel = 1./second ;
  static const double curie = 3.7e+10 * becquerel;
  static const double kilobecquerel = 1.e+3*becquerel;
  static const double megabecquerel = 1.e+6*becquerel;
  static const double gigabecquerel = 1.e+9*becquerel;
  static const double millicurie = 1.e-3*curie;
  static const double microcurie = 1.e-6*curie;
  static const double Bq = becquerel;
  static const double kBq = kilobecquerel;
  static const double MBq = megabecquerel;
  static const double GBq = gigabecquerel;
  static const double Ci = curie;
  static const double mCi = millicurie;
  static const double uCi = microcurie;

  //
  // Absorbed dose [L^2][T^-2]
  //
  static const double      gray = joule/kilogram ;
  static const double  kilogray = 1.e+3*gray;
  static const double milligray = 1.e-3*gray;
  static const double microgray = 1.e-6*gray;

  //
  // Luminous intensity [I]
  //
  static const double candela = 1.;

  //
  // Luminous flux [I]
  //
  static const double lumen = candela*steradian;

  //
  // Illuminance [I][L^-2]
  //
  static const double lux = lumen/(1000.*1000.);
*/
  //
  // Miscellaneous
  //
  static const double perCent     = 0.01 ;
  static const double perThousand = 0.001;
  static const double perMillion  = 0.000001;


