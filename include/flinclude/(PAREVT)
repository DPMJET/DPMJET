*                                                                      *
*=== parevt ===========================================================*
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     Copyright (C) 2003-2019:  CERN & INFN                            *
*     All Rights Reserved.                                             *
*                                                                      *
*     Created on 28-March-90 by A. Ferrari, INFN - Milan               *
*                                                                      *
*           included in:                                               *
*                                                                      *
*               flukam                                                 *
*               ancnuc                                                 *
*               anncan                                                 *
*               bdnopt                                                 *
*               bimsel                                                 *
*               bselsc                                                 *
*               clbmjt                                                 *
*               corevt                                                 *
*               corrin                                                 *
*               deflts                                                 *
*               difevv                                                 *
*               difpro                                                 *
*               difevv                                                 *
*               difpro                                                 *
*               distnu                                                 *
*               eldsev                                                 *
*               epilog                                                 *
*               eventv                                                 *
*               evevap                                                 *
*               evxtes                                                 *
*               evvini                                                 *
*               ferevv                                                 *
*               ferhav                                                 *
*               frmbrk                                                 *
*               frbkjp                                                 *
*               ferevv                                                 *
*               ferhav                                                 *
*               fiprou                                                 *
*               frmbrk                                                 *
*               frbkjp                                                 *
*               geta                                                   *
*               hdncin                                                 *
*               hdncgl                                                 *
*               hdncev                                                 *
*               kaonuc                                                 *
*               kbanuc                                                 *
*               kbsnuc                                                 *
*               kaskad                                                 *
*               kasneu                                                 *
*               incini                                                 *
*               lpncln                                                 *
*               muoabs                                                 *
*               muprab                                                 *
*               ncdtrd                                                 *
*               nucdcy                                                 *
*               nucnuc                                                 *
*               nucevv                                                 *
*               nuncln                                                 *
*               nwisel                                                 *
*               peaini                                                 *
*               peanut                                                 *
*               phnsch                                                 *
*               phncev                                                 *
*               phnvev                                                 *
*               piabth                                                 *
*               pioabs                                                 *
*               pionuc                                                 *
*               pmprab                                                 *
*               pnansg                                                 *
*               pnpisg                                                 *
*               ponfin                                                 *
*               precll                                                 *
*               prepre                                                 *
*               rcfset                                                 *
*               rsncli                                                 *
*               rstsel                                                 *
*               seafsl                                                 *
*               stglau                                                 *
*               stprnc                                                 *
*               umofin                                                 *
*               usrsrn                                                 *
*                                                                      *
*           frdiff = fraction of the total number of events which are  *
*                    diffractive ones (obsolete, no longer used)       *
*           fsprd0 = asymptotic reduction factor for cascade particle  *
*                    energy slope parameter                            *
*           fshpfn = modifying factor for bamjet momentum sharing jet  *
*                    function                                          *
*           rn1gsc = correlation factor between shower and grey part-  *
*                    ticle multiplicities                              *
*           rn2gsc = correlation factor between shower and grey part-  *
*                    ticle multiplicities                              *
*           rnswtc = random particle id switch factor                  *
*           dpower = power exponent used in the binomial distributions *
*                    for grey particles                                *
*           lpower = flag to signal if dpower is different from the    *
*                    default option (dpower=2)                         *
*           ldiffr = flag for particles undergoing diffractive events  *
*           linctv = flag to switch on/off cascade nucleons and exci-  *
*                    tation energy (it is important for testing the    *
*                    generator                                         *
*           levprt = flag for nuclear evaporation                      *
*           ldeexg = flag for nuclear deexcitation gammas              *
*           lgdhpr = flag for geometry dependent calculations in the   *
*                    preequilibrium part                               *
*           lpreex = flag for explicit treatment of the first inte-    *
*                    ctions in the preequilibrium part                 *
*           lhlfix = flag for "freezing" the hole depth of the holes   *
*                    produced in the explicitly treated interactions   *
*                    in the preequilibrium part                        *
*           lprfix = flag for "freezing" the already sampled particle  *
*                    exciton energy in the preequilibrium part         *
*           lparwv = flag to apply the Heisenberg principle according  *
*                    to the reduced DeBroglie wavelength in the calcu- *
*                    lation of the impact parameter and in computing   *
*                    the particle position after interactions          *
*           lsngch = logical flag for allowing single chain events     *
*           lschdf = logical flag for allowing single chain events     *
*                    in diffractive interactions                       *
*           lbmprq = logical flag for allowing preequilibrium in BME   *
*           lhadri = logical flag for indicating whether the resonance *
*                    model was used or not for primary interactions    *
*           lnucri = logical flag for indicating whether the Nucriv    *
*                    model was used or not for the current interaction *
*           lpeanu = logical flag for indicating whether the Peanut    *
*                    model was used or not for the current interaction *
*           levbme = logical flag for indicating whether the BME       *
*                    model was used or not for the current interaction *
*           levrqm = logical flag for indicating whether the rQMD      *
*                    model was used or not for the current interaction *
*           lphdrc = logical flag for indicating whether a photonuc-   *
*                    lear interaction is a "direct" one or not         *
*           latmss = logical flag for indicating whether to use ato-   *
*                    mic masses during evaporation                     *
*           lismrs = logical flag for indicating whether isomer calcu- *
*                    lations are possible                              *
*           lchdcy = logical flag for charm particle decays            *
*           lchdcr = logical flag for charm particle decays occured in *
*                    the current collision                             *
*           lmlccr = logical flag true when the current collisions is  *
*                    a multiple (Glauber like) collision               *
*           lrvkin = logical flag for interactions occuring in reverse *
*                    kinematics                                        *
*           lvp2xx = logical flag for interactions where vph2xx has    *
*                    been called                                       *
*           lv2xnw = logical flag for using vph2xx, set overall        *
*           lnwv2x = logical flag for using vph2xx, set event-by-event *
*           levfin = logical flag for event flag initialization        *
*                                                                      *
*----------------------------------------------------------------------*
*
      PARAMETER ( ETHSEA = 1.0D+00 )
*
      LOGICAL LDIFFR, LINCTV, LEVPRT, LHEAVY, LDEEXG, LGDHPR, LPREEX,
     &        LHLFIX, LPRFIX, LPARWV, LPOWER, LSNGCH, LSCHDF, LBMPRQ,
     &        LHADRI, LNUCRI, LPEANU, LEVBME, LEVRQM, LPHDRC, LATMSS,
     &        LISMRS, LCHDCY, LCHDCR, LMLCCR, LRVKIN, LVP2XX, LV2XNW,
     &        LNWV2X, LEVFIN
      COMMON / PAREVT / DPOWER, FSPRD0, FSHPFN, RN1GSC, RN2GSC, RNSWTC,
     &                  LDIFFR (NALLWP),LPOWER, LINCTV, LEVPRT, LHEAVY,
     &                  LDEEXG, LGDHPR, LPREEX, LHLFIX, LPRFIX, LPARWV,
     &                  LSNGCH, LSCHDF, LBMPRQ, LHADRI, LNUCRI, LPEANU,
     &                  LEVBME,	LEVRQM, LPHDRC, LATMSS, LISMRS, LCHDCY,
     &                  LCHDCR, LMLCCR, LRVKIN, LVP2XX, LV2XNW, LNWV2X,
     &                  LEVFIN
      SAVE / PAREVT /
