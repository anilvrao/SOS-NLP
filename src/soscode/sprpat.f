

      SUBROUTINE SPRPAT(IMODE,IROW,JCOL,NCOL,NROW,NONZA,IPU,TITLE)
C
C ======================================================================
C     SPRPAT===>sprpat   J.T. BETTS
C ======================================================================
C
      IMPLICIT DOUBLE PRECISION  (A-H,O-Z)
C
C         PURPOSE:
C
C                 DISPLAY THE SPARSITY PATTERN OF A MATRIX.  The
c                 pattern is displayed in rectangular panels with
c                 a maximum width of 100 characters, i.e. the 
c                 rectangles are of size (nrow x mxwid) where
c                 mxwid = max(ncol,100).  When ncol > 100, multiple
c                 panels are displayed.
C
C         ARGUMENTS:
C
C         IMODE   INTEGER FLAG DEFINING STORAGE SCHEME (SEE JCOL)
C                 AND DISPLAY FORMAT. IF  
C                 IMODE > 0   ONLY DISPLAY ROWS WITH NONZEROS
C                 IMODE < 0   DISPLAY ALL ROWS 
C
C         IROW    ROW INDEX OF NONZERO ELEMENT (NONZA)
C
C         JCOL    COLUMN STORAGE INDICATOR
C                 (|IMODE|=1) JCOL CONTAINS THE COLUMN START ARRAY (NCOL+1)
C                 (|IMODE|=2) JCOL CONTAINS THE COLUMN NUMBER (NONZA)
C
C         NCOL    NUMBER OF COLUMNS
C
C         NROW    NUMBER OF ROWS
C
C         NONZA   NUMBER OF NONZEROS
C
C         IPU     OUTPUT UNIT NO.
C                 IPU < 0   DISPLAY PATTERN ON UNIT NO. |IPU|
C                 IPU > 0   SAVE PATTERN FILE "SPRPAT.FIL" ON UNIT NO. IPU
C
C         TITLE   CHARACTER VARIABLE TITLE FOR OUTPUT
C
C-----------------------------------------------------------------------
C
      DIMENSION IROW(NONZA),JCOL(*)
C
      CHARACTER(LEN=100) TENDIG,ONEDIG,ROWPRT,BLANK,LYNE,HEADER
      CHARACTER(LEN=*) TITLE
      LOGICAL ROWPRN
      DATA TENDIG/'         1         2         3         4         5   
     $      6         7         8         9         0'/
      DATA BLANK(1:100) / ' '/
C
      ROWPRN = IMODE.LT.0
      IMODPR = ABS(IMODE)
      ONEDIG = REPEAT('1234567890',10)
      LYNE = REPEAT('=',100)
      LNTTLE = LEN(TITLE)
C
C-----------------------------------------------------------------------
C
C         DISPLAY SPARSITY PATTERN
C
      IF(IPU.LT.0) THEN
        IPUINT = -IPU
C
C         COUNT THE NUMBER OF FULL WIDTH PARTITIONS:
C         100 COLUMNS PER PARTITION
C
        NPART = NCOL/100
C
C         COUNT THE LENGTH OF THE LAST PARTITION
C
        LASTPR = NCOL - 100*NPART
C
        IF(LASTPR.GT.0) NPART = NPART + 1
C
        WRITE(IPUINT,2001) BLANK(1:100)
        WRITE(IPUINT,2001) LYNE(1:100)
        WRITE(IPUINT,2001) BLANK(1:100)
        WRITE(IPUINT,2001) TITLE(1:LNTTLE)
C
        DO K = 1,NPART
C
C        FIRST AND LAST COLUMN IN PARTITION K
C
          LFIRST = 1 + (K-1)*100 
          LLAST = LFIRST + 99
          LLAST = MIN(LLAST,NCOL)
C
C        WIDTH OF PARTITION
C
          LWIDTH = LLAST-LFIRST+1
C
          HEADER = BLANK
          HEADER(1:11) = '<--- COLUMN'
          WRITE(HEADER(12:21),'(I10)') LFIRST
          CALL HHCPRS(HEADER,' ',' ')
          WRITE(IPUINT,2001) BLANK(1:LWIDTH)
          WRITE(IPUINT,2001) LYNE(1:LWIDTH)
          WRITE(IPUINT,2001) BLANK(1:LWIDTH)
          WRITE(IPUINT,2001) HEADER(1:100)
          WRITE(IPUINT,2001) TENDIG(1:LWIDTH)
          WRITE(IPUINT,2001) ONEDIG(1:LWIDTH)
C
          DO I = 1,NROW
            ROWPRT = BLANK
C
            IF(IMODPR.EQ.1) THEN
                JK = 0
                DO JR = LFIRST,LLAST
                  JK = JK + 1
                  DO J = JCOL(JR),JCOL(JR+1)-1
                    IF(IROW(J).EQ.I) ROWPRT(JK:JK) = 'X'
                  enddo
                enddo
            ELSEIF(IMODPR.EQ.2) THEN
              DO II = 1,NONZA
                IF(IROW(II).EQ.I) THEN
                  JCL = JCOL(II)
                  IF(LFIRST.LE.JCL.AND.JCL.LE.LLAST) THEN
                    JK = JCL - LFIRST + 1
                    ROWPRT(JK:JK) = 'X'
                  ENDIF
                ENDIF
              enddo
            ENDIF
C
            IF(ROWPRT.NE.BLANK) THEN
              WRITE(IPUINT,2002) I,ROWPRT
            ELSEIF(ROWPRN) THEN
              WRITE(IPUINT,2002) I,ROWPRT
            ENDIF
C
          enddo
C
        enddo
C
      ELSE
        IPUINT = IPU
        OPEN(IPUINT,FILE='SPRPAT.FIL',STATUS='UNKNOWN')
        WRITE(IPUINT,2001) '-----'
        WRITE(IPUINT,2001) TITLE(1:LNTTLE)
        WRITE(IPUINT,2004) NONZA,NROW,NCOL
        IF(IMODPR.EQ.1) THEN
          WRITE(IPUINT,2003) ( (IROW(II),JC, 
     $       II = JCOL(JC),JCOL(JC+1)-1 ), JC = 1,NCOL )
        ELSEIF(IMODPR.EQ.2) THEN
          WRITE(IPUINT,2003) (IROW(II),JCOL(II),II = 1,NONZA)
        ENDIF
      ENDIF
C
      RETURN
 2001 FORMAT(T7,A)
 2002 FORMAT(T2,I4,1X,A100)
 2003 FORMAT(2X,2I10)
 2004 FORMAT(2X,3I10)
      END
