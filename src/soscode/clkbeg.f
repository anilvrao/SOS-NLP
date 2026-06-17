
      SUBROUTINE CLKBEG(NCLOCK)
C
C ======================================================================
C     CLKBEG===>clkbeg   J.T. BETTS
C ======================================================================
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C       PURPOSE:  BEGIN THE TIMING FOR CLOCK NUMBER NCLOCK
C
C       INPUT:
C
C         NCLOCK  NUMBER OF THE CLOCK TO BE STARTED
C
C                 NOTE IT IS ASSUMED THAT
C                 0 < NCLOCK LEQ MAXCLK
C                 THIS CONDITION IS NOT CHECKED!
C
      PARAMETER (MAXCLK=40)
C
      DOUBLE PRECISION INTCLK,INTOLD
      COMMON /CLKCOM/ CLKOLD(MAXCLK),CLKVAL(MAXCLK),INTCLK(MAXCLK),
     $    INTOLD(MAXCLK),IERCLK(MAXCLK)
C
      CALL XDSLT1(CLKOLD(NCLOCK),WLBEGN)
      IERCLK(NCLOCK) = +1
C
      RETURN
      END
