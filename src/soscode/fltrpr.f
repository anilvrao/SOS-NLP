


      SUBROUTINE FLTRPR(FTITLE,PENMU,PENRHO,FHFLTR,LNFLTR,MXFLTR)
C
C ======================================================================
C     FLTRPR===>FLTRPR   J.T. BETTS
C ======================================================================
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C         PURPOSE:  NLP FILTER PRINT ROUTINE 
C
C         INPUT:
C
C            PENMU  BARRIER PARAMETER
C            PENRHO RELAXATION PENALTY PARAMETER
C            FHFLTR ARRAY WITH FOUR COLUMNS, AND LNFLTR ROWS
C                   CONTAINING CURRENT FILTER VALUES
C            LNFLTR LENGTH OF THE CURRENT FILTER
C            MXFLTR MAXIMUM LENGTH OF THE FILTER
C        
      DIMENSION FHFLTR(MXFLTR,4)
C
      CHARACTER(LEN=3)  FTITLE
      CHARACTER(LEN=39) TITLE
C
C-------------------------------------------------------------
      INCLUDE '../commons/NLPSPR.CMN'
      INCLUDE '../commons/BARNLP.CMN'
C-------------------------------------------------------------
C
      IF(LNFLTR.EQ.0) RETURN
C
      WRITE(IPUNLP,1002)
      TITLE = '    Filter: B(y)           |c(y)|     '
      TITLE(1:3) = FTITLE
      WRITE(IPUNLP,1003) TITLE
      WRITE(IPUNLP,1002)
      DO I = 1,ABS(LNFLTR)
        FLTROB = FHFLTR(I,1) - PENMU*FHFLTR(I,3) + PENRHO*FHFLTR(I,4)
        WRITE(IPUNLP,1001) I,FLTROB,FHFLTR(I,2)
      enddo
      WRITE(IPUNLP,1002)
C
 1001 FORMAT(T3,'*',T33,'|',I3,' |',SP,1PE14.6,
     $    ' |',SP,1PE14.6,' |',T106,'*')
 1002 FORMAT(T3,'*',T33,'--------------------------------------'
     $    ,T106,'*')
 1003 FORMAT(T3,'*',T33,A39,T106,'*')
C
      RETURN
      END
