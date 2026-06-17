
      SUBROUTINE PRFLTR(FTITLE,FHFLTR,LNFLTR,MXFLTR)
C
C ======================================================================
C     PRFLTR===>prfltr   J.T. BETTS
C ======================================================================
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C         PURPOSE:  NLP FILTER PRINT ROUTINE 
C
C         INPUT:
C
C            FHFLTR ARRAY WITH THREE COLUMNS, AND LNFLTR ROWS
C                   CONTAINING CURRENT FILTER VALUES
C            LNFLTR LENGTH OF THE CURRENT FILTER
C            MXFLTR MAXIMUM LENGTH OF THE FILTER
C        
      DIMENSION FHFLTR(MXFLTR,3)
C
      CHARACTER(LEN=3) FTITLE
      CHARACTER(LEN=39) TITLE
C
C-------------------------------------------------------------
      INCLUDE '../commons/NLPSPR.CMN'
C-------------------------------------------------------------
C
      IF(LNFLTR.EQ.0) RETURN
C
      WRITE(IPUNLP,1002)
      TITLE = '    Filter: f(x)           |c(x)|     '
      TITLE(1:3) = FTITLE
      WRITE(IPUNLP,1003) TITLE
      WRITE(IPUNLP,1002)
      DO I = 1,ABS(LNFLTR)
        WRITE(IPUNLP,1001) I,FHFLTR(I,1),FHFLTR(I,2)
      ENDDO
      WRITE(IPUNLP,1002)
C
 1001 FORMAT(T3,'*',T33,'|',I3,' |',SP,1PE14.6,
     $    ' |',SP,1PE14.6,' |',T106,'*')
 1002 FORMAT(T3,'*',T33,'--------------------------------------'
     $    ,T106,'*')
 1003 FORMAT(T3,'*',T33,A39
     $    ,T106,'*')
C
      RETURN
      END
