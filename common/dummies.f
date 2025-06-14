#ifndef FOR_FLUKA

      SUBROUTINE EXPLOD ( NPEXPL, AMEXPL, ETOTEX, ETEXPL, PXEXPL,
     &              PYEXPL, PZEXPL )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      END

      SUBROUTINE GLAUBR(PPROJ,UMO,IBPROJ,IT,IP,info,Barr)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      END

      INTEGER FUNCTION MCIHAD(IDPDG)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      MCIHAD = IDT_ICIHAD(IDPDG)
      RETURN
      END

      INTEGER FUNCTION MPDGHA(I)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      MPDGHA = IDT_IPDGHA(I)
      RETURN
      END
#endif

      DOUBLE PRECISION FUNCTION JLL_SAMDSDT()
      IMPLICIT NONE
      RETURN
      END

      SUBROUTINE JLL_SET
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      RETURN
      END

#if defined(FOR_FLUKA)      
      DOUBLE PRECISION FUNCTION CLFLEV (WEE)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CALL EVEVAP (WEE)
      CLFLEV = WEE
      RETURN
      END
#endif
