
      SUBROUTINE RITLYN(IOUNIT,JC,A,NC)
C
C ======================================================================
C     RITLYN===>ritlyn   J.T. BETTS
C ======================================================================
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      DIMENSION A(5),JC(5)
      CHARACTER(LEN=106) LINE,BLANK
      DATA BLANK(1:106) / ' '/
C
      LINE = BLANK
      LINE(3:3) = '*'
      LINE(106:106) = '*'
      LBEG = 11
      DO I = 1,NC
        WRITE(LINE(LBEG:LBEG+3),'(I4)') JC(I)
        LINE(LBEG+4:LBEG+4) = ','
        WRITE(LINE(LBEG+5:LBEG+15),'(1PG11.4)') A(I)
        LBEG = LBEG + 18
      enddo
      WRITE(IOUNIT,1001) LINE
C
 1001 FORMAT(A106)
C
      RETURN
      END
