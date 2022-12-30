
      SUBROUTINE DT_XHOINI
C     SUBROUTINE DT_PHOINI
 
      IMPLICIT NONE
      SAVE 
 
#if defined(FLDOTINCL) && defined(FOR_FLUKA)
      INCLUDE 'inc/dtflka12ca'
#else
      INCLUDE 'inc/dtflka'
#endif
 
      END SUBROUTINE
