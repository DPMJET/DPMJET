C* This program demonstrates how to generate photo-hadronic events
C* with phojet without input cards.
C* Usage:
C*    ./bin/pho_hadronic 10000 (for 10k event)
C*
C* Anatoli Fedynitch, ICRR, (2021)

      PROGRAM PHOUPC
C**********************************************************************
C
C   example program calling PHOJET photon flux routines without cards
C
C**********************************************************************
      IMPLICIT NONE

      DOUBLE PRECISION ee , p1 , p2 , pcm , PHO_PMASS, pm1 , pm2 , s , 
     &                 sigcur , sigmax , sqs , ZERO
      INTEGER IARGC , id , irej , itry , k , neve, IDUM, i
      CHARACTER*15 PHO_PNAME

      EXTERNAL PHO_PMASS
      EXTERNAL PHO_PNAME
      EXTERNAL PYDATA

      SAVE 
 
      PARAMETER (ZERO=0.D0)
 
C  event debugging information
      INCLUDE 'inc/podebg'
 
C  standard particle data interface
      INCLUDE 'inc/poevt1'

C  photon flux kinematics and cuts
      INCLUDE 'inc/pofcut'

C  extension to standard particle data interface (PHOJET specific)
      INCLUDE 'inc/poevt2'
 
      DIMENSION p1(4) , p2(4)
      CHARACTER*72 title
      CHARACTER*32 arginp
 
C  *********** hp compiler settings ******************************
C     ON DOUBLE PRECISION UNDERFLOW IGNORE
C     ON REAL UNDERFLOW IGNORE
C  **************************************************************
 
C  number of events
      IF ( IARGC().GT.0 ) THEN
         CALL GETARG(1,arginp)
         READ (arginp,*) neve
      ELSE
         neve = 1
      END IF
C  general initialization of PHOJET data structures (mandatory)
C  (-2 means that no steering file is expected)
      irej = 20
      CALL PHO_INIT(-2,6,irej)

C Set (real) photon for Side 1 (left). Last parameter is the
C virtuality in case you need virtual photons
      CALL PHO_SETPAR(1, 22, 0, 0.0D0)
C Set proton target on Side 2
      CALL PHO_SETPAR(2, 2212, 0, 0.0D0)

C Photon from the left (max. 100 GeV initialization)
      p1(:) = (/0.0D0, 0.0D0, 100.0D0, 100.0D0/)
C Proton at rest
      p2(:) = (/0.0D0, 0.0D0, 0.0D0, PHO_PMASS(2212,1)/)

      CALL PHO_EVENT(-1,p1,p2,sigcur,irej)

C Photon from the left (10 GeV simulation energy)
      p1(:) = (/0.0D0, 0.0D0, 10.0D0, 10.0D0/)

      DO i = 0, neve
         CALL PHO_EVENT(1,p1,p2,sigcur,irej)
C        Print the event to study the details
         CALL PHO_PREVNT(1)
C        Loop through the particles and do the typical things with HEPEVT common blocks.
C        (See poevt1 include file)
         DO k = 0, NHEp
            IF (ISThep(k).EQ.1) THEN
               write(6, *) IDhep(k), PHO_PNAME(IDhep(k), 1), PHEp(:, k)
            ENDIF
         END DO
      END DO
C  some (optional) output of PHOJET-internal statistics
      CALL PHO_EVENT(-2,p1,p2,sigcur,irej)
 
      END PROGRAM
