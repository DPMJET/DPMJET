
      SUBROUTINE DT_XTIME
 
      IMPLICIT NONE
      SAVE 
 
#if defined(FLDOTINCL) && defined(FOR_FLUKA)
      INCLUDE 'inc/dtflka12ca'
#else
      INCLUDE 'inc/dtflka'
#endif
 
      CHARACTER dat*9 , tim*11
 
      dat = '         '
      tim = '           '
C     CALL GETDAT(IYEAR,IMONTH,IDAY)
C     CALL GETTIM(IHOUR,IMINUT,ISECND,IHSCND)
 
C     CALL DATE(DAT)
C     CALL TIME(TIM)
C     WRITE(LOUT,1000) DAT,TIM
C1000 FORMAT(/,2X,'Date: ',A9,3X,'Time: ',A11,/)
 
      END SUBROUTINE
