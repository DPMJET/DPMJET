      PROGRAM DPMTEST

      IMPLICIT DOUBLE PRECISION(A-H,O-Z)

C event history
      INCLUDE 'inc/dtevt1'
C extended event history
      INCLUDE 'inc/dtevt2'
C particle properties (BAMJET index convention)
      INCLUDE 'inc/dtpart'

      INCLUDE 'inc/dtflka'
      DIMENSION IDPLIST(6)
      DATA IDPLIST /2212,2112,-321,211,-211, 22/

      SAVE

      LOUT = 6
      LPRI = 20

C  general initialization 
      ECM = 50.D0
      NCASES = -1
      NPMASS = 1
      NPCHAR = 1
      NTMASS = 1
      NTCHAR = 1
      IDPIDX = 1
      IDP = IDPLIST(IDPIDX)
      IGLAU = 1
      ELAB = (ECM**2-AAM(1)**2-AAM(1)**2)/(2.0D0*AAM(1))
      ELABI = 1.D10
      PLAB = SQRT( (ELAB+AAM(1))*(ELAB-AAM(1)) )
C       write (LOUT,*) "initialization at",epn
C       CALL PHO_SETPDF(2212,idum, 10770, 0, 10, 0, -1)
      WRITE(LOUT,*) NCASES,ELAB,NPMASS,NPCHAR,NTMASS,NTCHAR,IDP,IGLAU
      Call dt_init(NCASES,ELABI,NPMASS,NPCHAR,NTMASS,NTCHAR,IDP,IGLAU)

C C     Debug ACTPDF
C       IDEB(2) = 30
C C     Debug SETPDF
C       IDEB(80) = 30
C C     Debug multiparticle
C       IDEB(90) = 30
C C     Debug PHO_PDF
C       ideb(1) = 25
C C     Debug PHO_HARINI   
C       ideb(LOUT6) = 25      
C C     Debug PHO_HARREM
C       ideb(28) = 20
C C     Debug PHO_CHAN2A
C       ideb(86) = 20
C     Debug PHO_CHAN2A
C       ideb(87) = 10
C     Debug PHO_REGPAR
C       ideb(55) = 40      
        
C       Do KD = 1, 100
C         ideb(kd) = 25
C       End Do
C  generate events
 54   CONTINUE  
      IDP = IDPLIST(IDPIDX)
      WRITE(LOUT,*) 'Particle ID:', IDP, IDT_ICIHAD(IDP)

C  number of events
      NEVE = 1000
      ITRY = 0
      IREJ = 0
      KKMAT = -1
      DO 100 K=1,NEVE
        NEVHKK = K
 55     CONTINUE
        ITRY = ITRY+1

C       Random choice of projectile        
        IRPRO = 1 + INT(DT_RNDM(DUM)*5.5D0)
        CALL dt_kkinc(NPMASS,NPCHAR,NTMASS,NTCHAR,
     &     IDT_ICIHAD(IDPLIST(IRPRO)),ELAB, KKMAT, IREJ)

        If (mod(k,1000)==0) Write(LOUT,*)  NEVENT, 'events simulated' 
        write(LOUT,*) "nevent,bimpac,IREJ", NEVENT, bimpac,IREJ
        IF(IREJ.NE.0) GOTO 55
C       Print event stack        
      !   CALL DT_EVTOUT(1)

C       event loop (Print particles stack)
        ICH = 0
        IBAR = 0
        DO I=1,NHKK
          IF(ISTHKK(I).EQ.1) THEN
            NEVENT = K
!             Write(LOUT,*) idhkk(I)
            ich = ich + ipho_chr3(idhkk(I),1)/3
            ibar = ibar + ipho_bar3(idhkk(I),1)/3
            PX = phkk(1,I)
            PY = phkk(2,I)
            PZ = phkk(3,I)
            EE = phkk(4,I)
            AM = phkk(5,I)
C            WRITE(LOUT,'(3I6,5E12.3)') nevhkk, idhkk(I),idbam(I),
C     &        PX, PY, PZ, EE, AM

          ENDIF
        ENDDO
 100  CONTINUE
      Call DT_STATIS(2)

      IDPIDX = IDPIDX + 1
       IF (IDPIDX.LE.LENIDPL) GOTO 54

      END
