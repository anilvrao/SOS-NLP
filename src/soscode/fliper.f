

      SUBROUTINE FLIPER(RELAX,
     $     NVARN,NONZCN,MSUBBN,NONZBN,NONZWN,
     $     NVARR,NONZCR,MSUBBR,NONZBR,NONZWR,
     $     NVAR,NONZC,MSUBB,NONZB,NONZW)
C
C ======================================================================
C     FLIPER===>fliper   J.T. BETTS
C ======================================================================
C
C     PURPOSE:
C        CHANGE FROM NORMAL TO RELAXATION OR RELAXATION TO NORMAL
C
C     INPUT:
C        RELAX     = TRUE    SET TO RELAXATION MODE
C                  = FALSE   SET TO NORMAL MODE
C                  ---NORMAL MODE INPUTS---
C        NVARN     NUMBER OF VARIABLES 
C        NONZCN    NUMBER OF NONZEROES IN CMAT 
C        MSUBBN    NUMBER OF BOUNDS 
C        NONZBN    NUMBER OF NONZEROES IN BMAT 
C        NONZWN    NUMBER OF NONZEROES IN WMAT 
C                  ---RELAXATION MODE INPUTS---
C        NVARR     NUMBER OF VARIABLES 
C        NONZCR    NUMBER OF NONZEROES IN CMAT 
C        MSUBBR    NUMBER OF BOUNDS 
C        NONZBR    NUMBER OF NONZEROES IN BMAT 
C        NONZWR    NUMBER OF NONZEROES IN WMAT 
C
C     OUTPUT:
C        NVAR      NUMBER OF VARIABLES 
C        NONZC     NUMBER OF NONZEROES IN CMAT 
C        MSUBB     NUMBER OF BOUNDS 
C        NONZB     NUMBER OF NONZEROES IN BMAT 
C        NONZW     NUMBER OF NONZEROES IN WMAT 
C
C     ******************************************************************
C  ARGUMENT LIST:
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL RELAX
C
C     ******************************************************************
C
      IF(RELAX) THEN
        NVAR  = NVARR 
        NONZC = NONZCR
        MSUBB = MSUBBR
        NONZB = NONZBR   
        NONZW = NONZWR
      ELSE
        NVAR  = NVARN 
        NONZC = NONZCN
        MSUBB = MSUBBN
        NONZB = NONZBN   
        NONZW = NONZWN
      ENDIF
C
      RETURN
      END
