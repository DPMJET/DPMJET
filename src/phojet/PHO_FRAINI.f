
      SUBROUTINE PHO_FRAINI(Idefau)
C***********************************************************************
C
C     initialization of fragmentation packages
C      (currently LUND JETSET)
C
C     initialization for JETSET call in DTUNUC 1.04 (J.R. 6/93)
C                      changed to work in PHOJET   (R.E. 1/94)
C
C     input:  IDEFAU    0  no hadronization at all
C                       1  do not touch any parameter of JETSET
C                       2  default parameters kept, decay length 10mm to
C                          define stable particles
C                       3  load tuned parameters for JETSET 7.3
C             neg. value:  prevent strange/charm hadrons from decaying
C
C***********************************************************************
      IMPLICIT NONE
      DOUBLE PRECISION def19 , def2 , def21 , def41 , def42 , EPS,
     &                 def1 , def3 , def5 , def7  
      INTEGER idef12 , idefab , Idefau , kc, idxsta, iunstab, NUNSTAB, i
      SAVE 
 
      PARAMETER (EPS=1.D-10)
 
C  input/output channels
      INCLUDE 'inc/poinou'
 
      INCLUDE 'inc/pyjets'
 
      INCLUDE 'inc/pydat1'
 
      INCLUDE 'inc/pydat2'
 
      INCLUDE 'inc/pydat3'
 
      INTEGER PYCOMP
      EXTERNAL PYCOMP

      DIMENSION idxsta(40)
C          K0s   pi0  lam   alam  sig+  asig+ sig-  asig- tet0  atet0
C          tet- atet-  om-  aom-   D+    D-    D0    aD0   Ds+   aDs+
C          etac lamc+alamc+sigc++ sigc+ sigc0asigc++asigc+asigc0 Ksic+
C         Ksic0 aKsic+aKsic0 sig0 asig0
      DATA idxsta/310 , 111 , 3122 , -3122 , 3222 , -3222 , 3112 , 
     &     -3112 , 3322 , -3322 , 3312 , -3312 , 3334 , -3334 , 411 , 
     &     -411 , 421 , -421 , 431 , -431 , 441 , 4122 , -4122 , 4222 , 
     &     4212 , 4112 , -4222 , -4212 , -4112 , 4232 , 4132 , -4232 , 
     &     -4132 , 3212 , -3212 , 5*0/
 
      PARAMETER (NUNSTAB=11)
      DIMENSION iunstab(NUNSTAB)
      DATA iunstab/4132 , 4232 , 4332 , 4312 , 4322 , 4324 , 4214 , 
     &     4224 , 4314 , 4114 , 10421/
 
      idefab = ABS(Idefau)
 
      IF ( idefab.EQ.0 ) THEN
         IF ( LPRi.GT.4 ) WRITE (LO,'(/1X,A)')
     &         'PHO_FRAINI: hadronization switched off'
         RETURN
      END IF
C  defaults, 
      def1 = PARj(1)
      def2 = PARj(2)
      def3 = PARj(3)
      def5 = PARj(5)
      def7 = PARj(7)
      def19 = PARj(19)
      def41 = PARj(41)
      def42 = PARj(42)
      def21 = PARj(21)
      idef12 = MSTj(12)
 
C  declare stable particles
      IF ( idefab.GE.2 ) MSTj(22) = 2
 
C  load optimized parameters
      IF ( idefab.GE.3 ) THEN
         PARj(1) = 0.09D0
         PARj(2) = 0.22D0
         PARj(3) = 0.9D0
         PARj(5) = 0.1D0
         PARj(7) = 0.95D0
         PARj(21) = 0.42D0
         PARj(41) = 0.3D0
         PARj(42) = 1.0D0
         MSTj(12)= 2
      END IF
C  Force pythia to follow the MDCy decay settings
      IF ( idefab.GE.4 ) THEN
         MSTj(21) = 1
         MSTj(22) = 1
      ENDIF
C prevent particles from decaying
      IF ( Idefau.LT.0) THEN
         DO i = 1 , 35
            kc = PYCOMP(idxsta(i))
            MDCy(kc,1) = 0
         END DO
      END IF
C prevent some charmed baryons which don't have a BAMJET code to be
C set as stable
      DO i = 1 , NUNSTAB
         kc = PYCOMP(iunstab(i))
         MDCy(kc,1) = 1
         kc = PYCOMP(-iunstab(i))
         MDCy(kc,1) = 1
      END DO
 
      IF ( LPRi.GT.4 ) WRITE (LO,99010) Idefau , def2 , PARj(2) , 
     &                        idef12 , MSTj(12) , def19 , PARj(19) , 
     &                        def41 , PARj(41) , def42 , PARj(42) , 
     &                        def21 , PARj(21)
99010 FORMAT (/' PHO_FRAINI: fragmentation initialization ISWMDL(6)',
     &        I3/,' --------------------------------------------------',
     &        /,5X,
     &        'parameter description               default / current',/,
     &        5X,'PARJ( 2) strangeness suppression : ',2F7.3,/,5X,
     &        'MSTJ(12) popcorn                 : ',2I7,/,5X,
     &        'PARJ(19) popcorn                 : ',2F7.3,/,5X,
     &        'PARJ(41) Lund a                  : ',2F7.3,/,5X,
     &        'PARJ(42) Lund b                  : ',2F7.3,/,5X,
     &        'PARJ(21) sigma in pt distribution: ',2F7.3,/)
 
      END SUBROUTINE
