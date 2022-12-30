
      SUBROUTINE DPMINI( Lesdpm, Loudpm, Ifdpm , Pinp  , Pbeam ,
     &                   Apfrmx, Iproa , Iproz , Idpmvr, Idpmhk,
     &                   Idpmfs, Ihflev )
 
C***********************************************************************
C                                                                      *
C     Version September 2001      by   Stefan Roesler                  *
C                                                                      *
C     Last change  14-Feb-22      by   Alfredo Ferrari                 *
C                                                                      *
C     This subroutine is part of the FLUKA interface to DPMJET 3,      *
C     new version (labelled 3.10 for convenience) after A.Fedynitch    *
C     work                                                             *
C     Initialization of DPMJET 3 event generation.                     *
C                                                                      *
C                                                                      *
C     argument list (FLUKA control card: codewd = dpmjet)              *
C                                                                      *
C          LESDPM  dpmjet input follows at the bottom of Fluka input   *
C                  file and a second start card is required in this    *
C                  case                                                *
C          LOUDPM  logical unit number for dpmjet output               *
C          IFDPM   this flag indicates if (and at what level?) we      *
C                  expect dpmjet output                                *
C                  for dpmjet-3 just on/off                            *
C          PINP    initial FLUKA beam momentum to indicate maximum     *
C                  available energy in the considered system           *
C          PBEAM   from FLUKA, used only to initialize histograms      *
C         APFRMX   from FLUKA, average (p/n) maximum Fermi momentum    *
C         IPROA/Z  from FLUKA, projectile A and Z                      *
C          IDPMVR  output flag indicating service is rendered by       *
C                  dpmjet-3                                            *
C          IDPMHK  carry over the dpmet(2/3) event common block size   *
C                  to enable important cross check since we have in    *
C                  principle 3 actual copies of it                     *
C          IDPMFS  string fusion flag                                  *
C          IHFLEV  using Fluka hadron-hadron event generator for       *
C                  (re)interactions flag                               *
C                                                                      *
C***********************************************************************
C
 
C  Do we trust removing??
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     SAVE
 
      DOUBLE PRECISION epn , Pinp , Pbeam, Apfrmx, xdumb , xlim1 ,
     &                 xlim2 , xlim3
      INTEGER ibin , idp , IDPmev , Idpmfs , Idpmhk , Idpmvr , Ifdpm , 
     &        iglau , iglaub , IHEhad , IHEnuc , IHIjpr , IHMapr , 
     &        IHMata , IHflev, IProa  , IProz  , Lesdpm , Loudpm ,
     &        ninp   , npchar , npmass
      INTEGER ntchar , ntmass
 
C emulsion treatment
      INCLUDE 'inc/dtcomp'
 
      INCLUDE 'inc/dtevt1'
C event flag
      INCLUDE 'inc/dtevno'
C Glauber formalism: flags and parameters for statistics
      INCLUDE 'inc/dtglgp'
C flags for input different options
      INCLUDE 'inc/dtflg1'
C flags for diffractive interactions (DTUNUC 1.x)
      INCLUDE 'inc/dtflg3'
C
      INCLUDE 'inc/pydat1'
C flags from from/for Fluka
#if defined(FLDOTINCL) && defined(FOR_FLUKA)
      INCLUDE 'inc/dtflka12ca'
#else
      INCLUDE 'inc/dtflka'
#endif
C histogram indices for Fluka-interface related statistics
C     CHARACTER*72 HEADER
      DIMENSION xdumb(40)
      COMMON /DTFLHX/ IHMapr , IHMata , IHIjpr , IHEnuc , IHEhad , 
     &                IDPmev
 
C
C Redirect Dpmjet and Pythia output and turn on/off output
C
      Idpmhk = NMXHKK
      Idpmvr = 310
      LPRi = Ifdpm
 
      LOUt = Loudpm
      MSTu(11) = Loudpm
C
C Flag for special settings needed to run the code as event
C  generator in Fluka (do not change !)
C
      ITRspt = 1
      IEMul = 0
      IFUsion = Idpmfs
      epn    = Pinp             !---> should we  * 1.5
      AMXpfr = APFrmx           !
      npmass = IPRoa            !     
      npchar = IPRoz            !---> last three are passed by FLUKA
C
C Special settings if no Dpmjet input follows the Fluka input
C
      IF ( Lesdpm.EQ.0 ) THEN
         ninp = -1
      ELSE
         ninp = 0
      END IF
 
C
C Initialization of Dpmjet
C
Cc
Cc ---------------------------------------------------------------------
 
      iglaub = 0
 
      CALL DT_INIT(ninp,epn,npmass,npchar,ntmass,ntchar,idp,iglau)
      CALL DT_STATIS(1)

      IFLEVG = IHFLEV
      IF ( ABS (IFLEVG) .LT. 1000 ) THEN
         ISIngd = 0
      ELSE
         ISIngd = IFLEVG / 1000
         IFLEVG = ABS ( MOD ( IFLEVG, 1000 ) )
      END IF
      IF ( IFLEVG .GT. 0 ) THEN
         EHADLO = EHFLLO
         EHADHI = EHFLHI
         EHADTH = -99.0D+00
      END IF
 
Cc ---------------------------------------------------------------------
Cc
C     CALL DT_DTUINI(NINP,EPN,NPMASS,NPCHAR,NTMASS,NTCHAR,IDP,IEMU)
      NEVent = 0
C
 
 
 
      IF ( LPRi.LT.3 ) RETURN
C
C Initialization of histograms
C
      IDPmev = 0
C  mass number of projectile and target nuclei
      xlim1 = 0.5D0
      xlim2 = 200.5D0
      xlim3 = 0D0
      ibin = INT(xlim2-xlim1)
      CALL DT_NEWHGR(xlim1,xlim2,xlim3,xdumb,ibin,IHMapr)
      CALL DT_NEWHGR(xlim1,xlim2,xlim3,xdumb,ibin,IHMata)
C  index of projectile hadrons
      CALL DT_NEWHGR(xlim1,xlim2,xlim3,xdumb,ibin,IHIjpr)
C  energy of projectile nuclei and hadrons
      xlim1 = 5.0D0
      xlim2 = PBEam
C     write(0,*) ' -pbeam-',pbeam
      xlim3 = 0.D0
      ibin = -200
      CALL DT_NEWHGR(xlim1,xlim2,xlim3,xdumb,ibin,IHEnuc)
      CALL DT_NEWHGR(xlim1,xlim2,xlim3,xdumb,ibin,IHEhad)
C
 
C=== End of subroutine Dpmini =========================================*
      END SUBROUTINE
