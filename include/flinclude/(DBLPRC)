*                                                                      *
*=== dblprc ===========================================================*
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     Copyright (C) 2003-2019:  CERN & INFN                            *
*     All Rights Reserved.                                             *
*                                                                      *
*      Dblprc: included in any routine, machine, mathematical and      *
*              physical constants plus global declarations             *
*                                                                      *
*  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   *
*  !!!! O N   M A C H I N E S   W H E R E   T H E   D O U B L E !!!!   *
*  !!!! P R E C I S I O N   I S   N O T   R E Q U I R E D  R E -!!!!   *
*  !!!! M O V E   T H E   D O U B L E   P R E C I S I O N       !!!!   *
*  !!!! S T A T E M E N T,  S E T   K A L G N M = 1   A N D     !!!!   *
*  !!!! C H A N G E   A L L   N U M E R I C A L   C O N S -     !!!!   *
*  !!!! T A N T S   T O   S I N G L E   P R E C I S I O N       !!!!   *
*  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   *
*                                                                      *
*         Kalgnm = real address alignment, 2 for double precision,     *
*                  1 for single precision                              *
*         Kalch8 = character*8 address alignment wrt the precision     *
*                  defined by Kalgnm (mostly 1 in all situations)      *
*         I2algn = integer*2 address alignment wrt the normal integer  *
*                  precision (mostly 2, 4 for 64 bit integers)         *
*         Anglgb = this parameter should be set equal to the machine   *
*                  "zero" with respect to unit                         *
*         Anglsq = this parameter should be set equal to the square    *
*                  of Anglgb                                           *
*         Axcssv = this parameter should be set equal to the number    *
*                  for which unity is negligible for the machine       *
*                  accuracy                                            *
*         Andrfl = "underflow" of the machine for floating point       *
*                  operation                                           *
*         Avrflw = "overflow"  of the machine for floating point       *
*                  operation                                           *
*         Ainfnt = code "infinite"                                     *
*         Azrzrz = code "zero"                                         *
*         Einfnt = natural logarithm of the code "infinite"            *
*         Ezrzrz = natural logarithm of the code "zero"                *
*         Excssv = natural logarithm of the code number for which      *
*                  unit is negligible                                  *
*         Englgb = natural logarithm of the code "zero" with respect   *
*                  to unit                                             *
*         Onemns = 1- of the machine, it is 1 - 2 x Anglgb             *
*         Onepls = 1+ of the machine, it is 1 + 2 x Anglgb             *
*         Csnnrm = maximum tolerable error on cosine normalization,    *
*                  u**2+v**2+w**2: assuming a typical anglgb relative  *
*                  error on each component we would get 2xanglgb: use  *
*                  4xanglgb to avoid too many normalizations           *
*         Dmxtrn = "infinite" distance for transport (cm)              *
*         Rhflmn = minimal density for Fluka (g/cm^3)                  *
*                                                                      *
*   "Global" declarations:                                             *
*         Lfluka = set to true for a real (full) Fluka run             *
*         Lgbias = set to true for a fully biased run                  *
*         Lgbana = set to true for a fully analogue run                *
*         Lflgeo = set to true when using the standard Fluka geometry  *
*         Loflts = set to true for special off-line testing of speci-  *
*                  fic routines                                        *
*         Lusrin = set to true if the user dependent initialization    *
*                  routine Usrini has been called at least onec        *
*         Lnmgeo = set to true for a name-based geometry input         *
*         Lnminp = set to true for a name-based Fluka input            *
*         Lfrfmt = set to true for a free-format based Fluka input     *
*         Ldpmcr = set to true for allowing dump cores inside Flabrt   *
*----------------------------------------------------------------------*
*                                                                      *
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER ( KALGNM = 2 )
      PARAMETER ( KALCH8 = 1 )
      PARAMETER ( I2ALGN = 2 )
      PARAMETER ( ANGLGB = 5.0D-16 )
      PARAMETER ( ANGLSQ = 2.5D-31 )
      PARAMETER ( AXCSSV = 0.2D+16 )
      PARAMETER ( ANDRFL = 1.0D-38 )
      PARAMETER ( AVRFLW = 1.0D+38 )
      PARAMETER ( AINFNT = 1.0D+30 )
      PARAMETER ( AZRZRZ = 1.0D-30 )
      PARAMETER ( EINFNT = +69.07755278982137 D+00 )
      PARAMETER ( EZRZRZ = -69.07755278982137 D+00 )
      PARAMETER ( EXCSSV = +35.23192357547063 D+00 )
      PARAMETER ( ENGLGB = -35.23192357547063 D+00 )
      PARAMETER ( ONEMNS = 0.999999999999999  D+00 )
      PARAMETER ( ONEPLS = 1.000000000000001  D+00 )
      PARAMETER ( CSNNRM = 2.0D-15 )
      PARAMETER ( DMXTRN = 1.0D+15 * 1.000001378D+00 )
      PARAMETER ( RHFLMN = 1.0D-10 )
*      PARAMETER ( DMXTRN = 1.0D+15 )
*      PARAMETER ( RHFLMN = 1.0D-13 )
*
*======================================================================*
*======================================================================*
*=========                                                   ==========*
*=========    M A T H E M A T I C A L   C O N S T A N T S    ==========*
*=========                                                   ==========*
*======================================================================*
*======================================================================*
*                                                                      *
*   Numerical constants (long integer):                                *
*                                                                      *
*         I8zero = "0" for Integer*8                                   *
*                                                                      *
*   Numerical constants (single precision):                            *
*                                                                      *
*         Zersng = 0                                                   *
*                                                                      *
*   Numerical constants (double precision):                            *
*                                                                      *
*         Zerzer = 0                                                   *
*         Oneone = 1                                                   *
*         Twotwo = 2                                                   *
*         Thrthr = 3                                                   *
*         Foufou = 4                                                   *
*         Fivfiv = 5                                                   *
*         Sixsix = 6                                                   *
*         Sevsev = 7                                                   *
*         Eigeig = 8                                                   *
*         Aninen = 9                                                   *
*         Tenten = 10                                                  *
*         Eleven = 11                                                  *
*         Twelve = 12                                                  *
*         Fiften = 15                                                  *
*         Sixten = 16                                                  *
*         Hlfhlf = 1/2                                                 *
*         Onethi = 1/3                                                 *
*         Onefou = 1/4                                                 *
*         Onefiv = 1/5                                                 *
*         Onesix = 1/6                                                 *
*         Onesev = 1/7                                                 *
*         Oneeig = 1/8                                                 *
*         Twothi = 2/3                                                 *
*         Thrfou = 3/4                                                 *
*         Thrtwo = 3/2                                                 *
*         Fouthr = 4/3                                                 *
*         Pipipi = Circumference / diameter                            *
*         Twopip = 2 x Pipipi                                          *
*         Pip5o2 = 5/2 x Pipipi                                        *
*         Pipisq = Pipipi x Pipipi                                     *
*         Pihalf = 1/2 x Pipipi                                        *
*         R3tovl = 4/3 x Pipipi                                        *
*         Erfa00 = Erf (oo) = 1/2 x square root of pi                  *
*         Sqrtpi = square root of pi                                   *
*         Sqtwpi = square root of 2xpi                                 *
*         Eulero = Eulero's constant                                   *
*         Eulexp = exp ( Eulero )                                      *
*         E1m2eu = exp ( 1 - 2 eulero )                                *
*         Eneper = "e", base of natural logarithm                      *
*         Sqrent = square root of "e"                                  *
*         Sqrtwo = square root of  2                                   *
*         Sqrthr = square root of  3                                   *
*         Sqrfiv = square root of  5                                   *
*         Sqrsix = square root of  6                                   *
*         Sqrsev = square root of  7                                   *
*         Sqrt12 = square root of 12                                   *
*         S2fwhm = 2 x square root of 2 x logarithm of 2               *
*         Twolog = logarithm of 2                                      *
*         Two2o3 = 2^2/3                                               *
*         Tenlog = logarithm of 10                                     *
*         Atnfou = Atan[4]                                             *
*         Atnsix = Atan[6]                                             *
*         Alpemd = 2 / exp(Eulero) / Sqrt[e]                           *
*                                                                      *
*----------------------------------------------------------------------*
*
      INTEGER*8 I8ZERO
      REAL ZERSNG
      PARAMETER ( I8ZERO =      0 )
      PARAMETER ( ZERSNG = 0.E+00 )
      PARAMETER ( ZERZER = 0.D+00 )
      PARAMETER ( ONEONE = 1.D+00 )
      PARAMETER ( TWOTWO = 2.D+00 )
      PARAMETER ( THRTHR = 3.D+00 )
      PARAMETER ( FOUFOU = 4.D+00 )
      PARAMETER ( FIVFIV = 5.D+00 )
      PARAMETER ( SIXSIX = 6.D+00 )
      PARAMETER ( SEVSEV = 7.D+00 )
      PARAMETER ( EIGEIG = 8.D+00 )
      PARAMETER ( ANINEN = 9.D+00 )
      PARAMETER ( TENTEN = 10.D+00 )
      PARAMETER ( ELEVEN = 11.D+00 )
      PARAMETER ( TWELVE = 12.D+00 )
      PARAMETER ( FIFTEN = 15.D+00 )
      PARAMETER ( SIXTEN = 16.D+00 )
      PARAMETER ( HLFHLF = 0.5D+00 )
      PARAMETER ( ONETHI = ONEONE / THRTHR )
      PARAMETER ( ONEFOU = ONEONE / FOUFOU )
      PARAMETER ( ONEFIV = ONEONE / FIVFIV )
      PARAMETER ( ONESIX = ONEONE / SIXSIX )
      PARAMETER ( ONESEV = ONEONE / SEVSEV )
      PARAMETER ( ONEEIG = ONEONE / EIGEIG )
      PARAMETER ( TWOTHI = TWOTWO / THRTHR )
      PARAMETER ( THRFOU = THRTHR / FOUFOU )
      PARAMETER ( THRTWO = THRTHR / TWOTWO )
      PARAMETER ( FOUTHR = FOUFOU / THRTHR )
      PARAMETER ( PIPIPI = 3.141592653589793238462643383279D+00 )
      PARAMETER ( TWOPIP = 6.283185307179586476925286766559D+00 )
      PARAMETER ( PIP5O2 = 7.853981633974483096156608458199D+00 )
      PARAMETER ( PIPISQ = 9.869604401089358618834490999876D+00 )
      PARAMETER ( PIHALF = 1.570796326794896619231321691640D+00 )
      PARAMETER ( R3TOVL = FOUFOU * PIPIPI / THRTHR )
      PARAMETER ( ERFA00 = 0.886226925452758013649083741671D+00 )
      PARAMETER ( SQRTPI = 1.772453850905516027298167483341D+00 )
      PARAMETER ( SQTWPI = 2.506628274631000502415765284811D+00 )
      PARAMETER ( EULERO = 0.577215664901532860606512      D+00 )
      PARAMETER ( EULEXP = 1.781072417990197985236504      D+00 )
      PARAMETER ( EULLOG =-0.5495393129816448223376619     D+00 )
      PARAMETER ( E1M2EU = 0.8569023337737540831433017     D+00 )
      PARAMETER ( ENEPER = 2.718281828459045235360287471353D+00 )
      PARAMETER ( SQRENT = 1.648721270700128146848650787814D+00 )
      PARAMETER ( SQRTWO = 1.414213562373095048801688724210D+00 )
      PARAMETER ( SQRTHR = 1.732050807568877293527446341506D+00 )
      PARAMETER ( SQRFIV = 2.236067977499789696409173668731D+00 )
      PARAMETER ( SQRSIX = 2.449489742783178098197284074706D+00 )
      PARAMETER ( SQRSEV = 2.645751311064590590501615753639D+00 )
      PARAMETER ( SQRT12 = 3.464101615137754587054892683012D+00 )
      PARAMETER ( S2FWHM = 2.354820045030949382023138652919D+00 )
      PARAMETER ( TWOLOG = 0.693147180559945309417232121458D+00 )
      PARAMETER ( TWO2O3 = 1.587401051968199474751705639272D+00 )
      PARAMETER ( TENLOG = 2.302585092994045684017991454684D+00 )
      PARAMETER ( ATNFOU = 1.3258176636680326D+00 )
      PARAMETER ( ATNSIX = 1.4056476493802699D+00 )
      PARAMETER ( ALPEMD = TWOTWO / EULEXP /SQRENT )
*
*======================================================================*
*======================================================================*
*=========                                                   ==========*
*=========       P H Y S I C A L   C O N S T A N T S         ==========*
*=========                                                   ==========*
*======================================================================*
*======================================================================*
*                                                                      *
*   Primary constants:                                                 *
*                                                                      *
*         Clight = speed of light in cm s-1                            *
*         Avogad = Avogadro number                                     *
*         Boltzm = k Boltzmann constant (J K-1)                        *
*         Amelgr = electron mass (g)                                   *
*         Plckbr = reduced Planck constant (erg s)                     *
*         Elccgs = elementary charge (CGS unit)                        *
*         Elcmks = elementary charge (MKS unit)                        *
*         Pvcmu0 = vacuum permeability constant (MKS unit)             *
*         Epsil0 = vacuum dielectric constant (MKS unit)               *
*         Amugrm = Atomic mass unit (g)                                *
*         Ammumu = Muon    mass (amu)                                  *
*         Amprmu = Proton  mass (amu)                                  *
*         Amnemu = Neutron mass (amu)                                  *
*                                                                      *
*   Derived constants:                                                 *
*                                                                      *
*         Alamb0 = Compton wavelength = 2 pi r0 / fsc                  *
*         Alpfsc = Fine structure constant  = e^2/(hbar c) (CGS units) *
*         Note:                                                        *
*          e[MKS]^2/ (4 pi epsilon0 hbar[MKS] c[MKS]) = alpha = 1/137  *
*          e[CGS]^2/ (hbar[CGS] c[CGS]) = alpha = 1/137                *
*          with c[MKs] = 10^-2 c[CGS], hbar[MKS] = 10^-7 hbar[CGS]     *
*         Amelct = Electron mass (GeV) = 10^-16Amelgr Clight^2 / Elcmks*
*         Amugev = Atomic mass unit (GeV) = 10^-16Amugrm Clight^2      *
*                                           / Elcmks                   *
*         Ammuon = Muon    mass (GeV) = Ammumu * Amugev                *
*         Amprtn = Proton  mass (GeV) = Amprmu * Amugev                *
*         Amntrn = Neutron mass (GeV) = Amnemu * Amugev                *
*         Amdeut = Deuteron mass (GeV)                                 *
*         Amtrtn = Triton   mass (GeV) (derived from the CODATA value  *
*                  in kg, and the amu value in kg->GeV)                *
*         Am3hel = 3-He     mass (GeV) (derived from the CODATA value  *
*                  in kg, and the amu value in kg->GeV)                *
*         Amalph = Alpha    mass (GeV) (derived from the excess mass   *
*                  and an (approximate) atomic binding not a really    *
*                  measured constant)                                  *
*         Ampiop = Pion+  mass (GeV)                                   *
*         Ampiom = Pion-  mass (GeV)                                   *
*         Ampio0 = Pion0  mass (GeV)                                   *
*         Cougfm = e^2 (GeV fm) = Elccgs^2 / Elcmks * 10^-7 * 10^-9    *
*                * 10^13 (10^..=erg cm->joule cm->GeV cm->GeV fm       *
*                  it is equal to 0.00144 GeV fm                       *
*         Fscto2 = (Fine structure constant)^2                         *
*         Fscto3 = (Fine structure constant)^3                         *
*         Fscto4 = (Fine structure constant)^4                         *
*         Plabrc = Reduced Planck constant times the light velocity    *
*                  expressed in GeV fm                                 *
*         Rclsel = Classical electron radius (cm) = e^2 / (m_e c^2)    *
*         Bltzmn = k Boltzmann constant in GeV K-1                     *
*         A0bohr = Bohr radius, hbar^2 / ( m_e e^2) (fm) = Plabrc**2   *
*                / Amelct / Cougfm, or equivalently,                   *
*                Plabrc / Alpfsc / Amelct                              *
*         Gfohb3 = Fermi constant, G_f/(hbar c)^3, in GeV^-2           *
*         Gfermi = Fermi constant in GeV fm^3                          *
*         Sin2tw = sin^2 theta_Weinberg                                *
*         Prmgnm = proton  magnetic moment (magneton)                  *
*         Anmgnm = neutron magnetic moment (magneton)                  *
*         S0Thms = sigma_0 Thomson, 8/3 pi r_e^2 (mb)                  *
*         Plfrcs = plasma frequency constant = 28.816 eV               *
*                = Sqrt[ 4 pi alpha hbar^3 / m_e] with hbar in [eV cm] *
*                  and m_e in [eV]                                     *
*                                                                      *
*   Astronomical constants:                                            *
*                                                                      *
*         Rearth = Earth equatorial radius (cm)                        *
*         Auastu = Astronomical Unit       (cm)                        *
*                                                                      *
*   Conversion constants:                                              *
*                                                                      *
*         GeVMeV = from GeV to MeV                                     *
*         eMVGeV = from MeV to GeV                                     *
*         GeV2eV = from GeV to  eV                                     *
*         eV2GeV = from eV  to GeV                                     *
*         alGVMV = from GeV to MeV, log                                *
*         Raddeg = from radians to degrees                             *
*         Degrad = from degrees to radians                             *
*         GeVOmg = from (photon) energy [GeV] in 2pi x frequency [s^-1]*
*         Cmq2mb = from square centimetres to millibarns               *
*         Fmb2ba = from millibarns to barns                            *
*         Bar2mb = from barns to millibarns                            *
*         Fmb2fs = from millibarns to fermi^2                          *
*         Fms2mb = from fermi^2 to millibarns                          *
*         Ba2cmq = from barns to square centimetres                    *
*         Cmq2ba = from square centimetres to barns                    *
*         Fmi2cm = from fermi to cm                                    *
*         Cm2fmi = from cm to fermi                                    *
*                                                                      *
*   Useful constants:                                                  *
*                                                                      *
*         Fertho = constant to be used in the Fermi-Thomas approxima-  *
*                  ted expression for atomic binding energies          *
*         Expebn = exponent to be used in the Fermi-Thomas approxima-  *
*                  ted expression for atomic binding energies          *
*                    B_atomic (Z) = Fertho x Z^Expebn (GeV)            *
*         Bex1hy = 1-H atomic binding energies (GeV)                   *
*         Bex2hy = 2-H atomic binding energy   (GeV) as obtained from  *
*                  the 2-H (nuclear) mass, the 2-H atomic excess mass  *
*                  and the amu and electron masses. Consistent within  *
*                  the p and 2-H mass uncertainties with Bex1hy        *
*         Bex3hy = 3-H atomic binding energy   (GeV) as obtained from  *
*                  the 3-H (nuclear) mass, the 3-H atomic excess mass  *
*                  and the amu and electron masses. Consistent within  *
*         Bex3he = 3-He atomic binding energy  (GeV) as obtained from  *
*                  the 3-He (nuclear) mass, the 3-He atomic excess mass*
*                  and the amu and electron masses. Consistent within  *
*                  the alpha and 3-He mass uncertainties with Bex4he   *
*         Bex4he = 4-He atomic binding energy  (GeV) as obtained from  *
*                  the 4-He (nuclear) mass, the 4-He atomic excess mass*
*                  and the amu and electron masses. Consistent within  *
*                  the alpha and 3-He mass uncertainties with Bex4he   *
*         Bexc12 = Fermi-Thomas approximated expression for 12-C ato-  *
*                  mic binding energies (GeV)                          *
*         Bhfzxx = Hartree-Fock atomic binding energies for Z=xx, from *
*                  Adndt18 243 (1976) and Density-Functional Theory of *
*                  Atoms and Molecules, R.G.Parr and W.Yang            *
*                  Please note that the original values are in units   *
*                  of 2 x Bex1hy                                       *
*         Amunmu = difference between the atomic and nuclear mass units*
*         Amuc12 = "Nuclear" mass unit = 1/12 M_nucl (12-C),           *
*                  M_nucl (12-C) = M_atom (12-C) - 6 m_e + B_atom(12-C)*
*         Ex1hmv = 1-H (atomic) excess mass in MeV,                    *
*                  [(Amprtn - Amugev) + Amelct - Bex1hy] x GeVMeV      *
*         Exnemv = Neutron excess mass in MeV,                         *
*                  [Amntrn - Amugev] x GeVMeV                          *
*                                                                      *
*----------------------------------------------------------------------*
*
*  Clight is exact:
      PARAMETER ( CLIGHT = 2.99792458         D+10 )
      PARAMETER ( AVOGAD = 6.0221367          D+23 )
      PARAMETER ( BOLTZM = 1.380658           D-23 )
      PARAMETER ( AMELGR = 9.1093897          D-28 )
      PARAMETER ( PLCKBR = 1.05457266         D-27 )
      PARAMETER ( ELCCGS = 4.8032068          D-10 )
      PARAMETER ( ELCMKS = 1.60217733         D-19 )
      PARAMETER ( AMUGRM = 1.6605402          D-24 )
      PARAMETER ( AMMUMU = 0.113428913        D+00 )
      PARAMETER ( AMPRMU = 1.007276470        D+00 )
      PARAMETER ( AMNEMU = 1.008664904        D+00 )
      PARAMETER ( EPSIL0 = 8.854187817        D-12 )
*     PARAMETER ( ALPFSC = 1.D+00 / 137.035989561D+00 )
*  Pvcmu0 is exact:
      PARAMETER ( PVCMU0 = TWOTWO * TWOPIP * 1.D-07 )
*  Epsilon0 is exact since c has been fixed:
*     PARAMETER ( FSCTO2 = ALPFSC * ALPFSC )
*     PARAMETER ( FSCTO3 = FSCTO2 * ALPFSC )
*     PARAMETER ( FSCTO4 = FSCTO3 * ALPFSC )
*    It is important to set the electron mass exactly with the same
*    rounding as in the mass tables, so use the explicit expression
*     PARAMETER ( AMELCT = 1.D-16 * AMELGR * CLIGHT * CLIGHT / ELCMKS )
*    It is important to set the amu mass exactly with the same
*    rounding as in the mass tables, so use the explicit expression
*     PARAMETER ( AMUGEV = 1.D-16 * AMUGRM * CLIGHT * CLIGHT / ELCMKS )
*    It is important to set the muon,proton,neutron masses exactly with
*    the same rounding as in the mass tables, so use the explicit
*    expression
*     PARAMETER ( AMMUON = AMMUMU * AMUGEV )
*     PARAMETER ( AMPRTN = AMPRMU * AMUGEV )
*     PARAMETER ( AMNTRN = AMNEMU * AMUGEV )
*     PARAMETER ( RCLSEL = ELCCGS * ELCCGS / CLIGHT / CLIGHT / AMELGR )
*     PARAMETER ( BLTZMN = BOLTZM / ELCMKS * 1.D-09 )
      PARAMETER ( ALPFSC = 7.2973530791728595 D-03 )
      PARAMETER ( FSCTO2 = 5.3251361962113614 D-05 )
      PARAMETER ( FSCTO3 = 3.8859399018437826 D-07 )
      PARAMETER ( FSCTO4 = 2.8357075508200407 D-09 )
      PARAMETER ( PLABRC = 0.197327053        D+00 )
      PARAMETER ( AMELCT = 0.51099906         D-03 )
      PARAMETER ( AMUGEV = 0.93149432         D+00 )
      PARAMETER ( AMMUON = 0.105658389        D+00 )
      PARAMETER ( AMPRTN = 0.93827231         D+00 )
      PARAMETER ( AMNTRN = 0.93956563         D+00 )
      PARAMETER ( AMDEUT = 1.87561339         D+00 )
      PARAMETER ( AMTRTN = 2.808921005        D+00 )
      PARAMETER ( AM3HEL = 2.808391482        D+00 )
*      PARAMETER ( AMALPH = 3.727379239        D+00 )
      PARAMETER ( AMALPH = 3.72738025692891   D+00 )
      PARAMETER ( AMPIOP = 1.3957018          D-01 )
      PARAMETER ( AMPIOM = AMPIOP                  )
      PARAMETER ( AMPIO0 = 1.349766           D-01 )
      PARAMETER ( COUGFM = ELCCGS * ELCCGS / ELCMKS * 1.D-07 * 1.D+13
     &                   * 1.D-09 )
*      PARAMETER ( RCLSEL = 2.8179403267       D-13 )
      PARAMETER ( RCLSEL = 2.8179409183694872 D-13 )
      PARAMETER ( ALAMB0 = TWOTWO * PIPIPI * RCLSEL / ALPFSC )
*      PARAMETER ( BLTZMN = 8.6173324          D-14 )
      PARAMETER ( BLTZMN = 8.617385           D-14 )
      PARAMETER ( A0BOHR = PLABRC / ALPFSC / AMELCT )
*      PARAMETER ( GFOHB3 = 1.1663787          D-05 )
      PARAMETER ( GFOHB3 = 1.16639            D-05 )
      PARAMETER ( GFERMI = GFOHB3 * PLABRC * PLABRC * PLABRC )
*      PARAMETER ( SIN2TW = 0.23126            D+00 )
*      PARAMETER ( SIN2TW = 0.23155            D+00 )
*      PARAMETER ( PRMGNM = 2.792847356        D+00 )
*      PARAMETER ( ANMGNM =-1.9130427          D+00 )
      PARAMETER ( SIN2TW = 0.2319             D+00 )
      PARAMETER ( PRMGNM = 2.792847386        D+00 )
      PARAMETER ( ANMGNM =-1.91304275         D+00 )
      PARAMETER ( REARTH = 6.378140           D+08 )
      PARAMETER ( AUASTU = 1.4959787066       D+13 )
      PARAMETER ( GEVMEV = 1.0                D+03 )
      PARAMETER ( EV2GEV = 1.0                D-09 )
      PARAMETER ( GEV2EV = 1.0                D+09 )
      PARAMETER ( EMVGEV = 1.0                D-03 )
      PARAMETER ( CMQ2MB = 1.0                D+27 )
      PARAMETER ( FMB2BA = 1.0                D-03 )
      PARAMETER ( BAR2MB = 1.0                D+03 )
      PARAMETER ( FMB2FS = 1.0                D-01 )
      PARAMETER ( FMS2MB = 1.0                D+01 )
      PARAMETER ( BA2CMQ = 1.0                D-24 )
      PARAMETER ( CMQ2BA = 1.0                D+24 )
      PARAMETER ( FMI2CM = 1.0                D-13 )
      PARAMETER ( CM2FMI = 1.0                D+13 )
      PARAMETER ( AKG2GR = 1.0                D+03 )
      PARAMETER ( ALGVMV = 6.90775527898214   D+00 )
      PARAMETER ( RADDEG = 180.D+00 / PIPIPI )
      PARAMETER ( DEGRAD = PIPIPI / 180.D+00 )
      PARAMETER ( GEVGRA = ELCMKS * GEV2EV * AKG2GR )
      PARAMETER ( GEVOMG = CLIGHT * 1.D+13 / PLABRC )
      PARAMETER ( S0THMS = EIGEIG / THRTHR * PIPIPI * RCLSEL * RCLSEL
     &                   * CMQ2MB )
*  Plasma frequency constant, to be recomputed every time alpha, m_e,
*  hbar or N_av change, it requires a Sqrt ergo it cannot be
*  parametrized explicitly as a function of the above constants
*     PARAMETER ( PLFRCS = SQRT ( ALPFSC * FOUFOU * PIPIPI
*    &                   * ( PLABRC * FMI2CM )**3 / AMELCT * AVOGAD )
*    &                   * GEV2EV )
      PARAMETER ( PLFRCS = 2.8815937226648276D+01 )
*  Old Fermi-Thomas parametrization of atomic binding energies:
*     PARAMETER ( FERTHO = 15.73       D-09 )
*     PARAMETER ( EXPEBN = 7.D+00 / 3.D+00  )
*     PARAMETER ( BEXC12 = FERTHO * 65.41634134195703D+00 )
*  New Fermi-Thomas parametrization of atomic binding energies:
      PARAMETER ( FERTHO = 14.33       D-09 )
      PARAMETER ( EXPEBN = 2.39        D+00 )
      PARAMETER ( BEX1HY = 13.61       D-09 )
      PARAMETER ( BEX2HY = 14.188      D-09 )
      PARAMETER ( BEX3HY = 14.818      D-09 )
      PARAMETER ( BEX3HE = 81.346      D-09 )
      PARAMETER ( BEX4HE = 77.246      D-09 )
      PARAMETER ( BHFZ02 = TWOTWO * BEX1HY * 2.8615  D+00 )
      PARAMETER ( BHFZ03 = TWOTWO * BEX1HY * 7.4290  D+00 )
      PARAMETER ( BHFZ04 = TWOTWO * BEX1HY * 14.569  D+00 )
      PARAMETER ( BHFZ05 = TWOTWO * BEX1HY * 24.530  D+00 )
      PARAMETER ( BEXC12 = FERTHO * 72.40715579499394D+00 )
      PARAMETER ( BHFZ06 = TWOTWO * BEX1HY * 37.645  D+00 )
      PARAMETER ( BHFZ07 = TWOTWO * BEX1HY * 54.303  D+00 )
      PARAMETER ( BHFZ08 = TWOTWO * BEX1HY * 74.801  D+00 )
      PARAMETER ( BHFZ09 = TWOTWO * BEX1HY * 99.471  D+00 )
      PARAMETER ( BHFZ10 = TWOTWO * BEX1HY * 128.55  D+00 )
      PARAMETER ( BHFZ18 = TWOTWO * BEX1HY * 526.82  D+00 )
      PARAMETER ( BHFZ36 = TWOTWO * BEX1HY * 2752.0  D+00 )
      PARAMETER ( BHFZ54 = TWOTWO * BEX1HY * 7232.1  D+00 )
*      PARAMETER ( AMUNMU = HLFHLF * AMELCT - BHFZ06 / 12.D+00 )
      PARAMETER ( AMUNMU = HLFHLF * AMELCT - BEXC12 / 12.D+00 )
      PARAMETER ( AMUC12 = AMUGEV - AMUNMU )
*  p and n excess mass in MeV:
      PARAMETER ( EX1HMV = ( ( AMPRTN - AMUGEV ) + AMELCT - BEX1HY )
     &                       * GEVMEV )
      PARAMETER ( EXNEMV = ( AMNTRN - AMUGEV ) * GEVMEV )
*  Old MeV units:
      PARAMETER ( AMEMEV = GEVMEV * AMELCT )
*  Infinite and "zero" T1/2:
      PARAMETER ( T12INF = 1.D+30 )
      PARAMETER ( T12ZER = 1.D-15 )
*
      LOGICAL LFLUKA, LGBIAS, LGBANA, LFLGEO, LOFLTS, LUSRIN,
     &        LUSRGL, LNMGEO, LNMINP, LFRFMT, LDMPCR
      LOGICAL LLDFLR, LDPBDY
      COMMON / GLOBAL / LFLUKA, LGBIAS, LGBANA, LFLGEO, LOFLTS, LUSRIN,
     &                  LUSRGL, LNMGEO, LNMINP, LFRFMT, LDMPCR,
     &                  LLDFLR, LDPBDY,
     &                  KFLGEO, KFLDNR
      COMMON / GLOBCH / CRVRFL
      CHARACTER*8 CRVRFL
      SAVE / GLOBAL /, / GLOBCH /
