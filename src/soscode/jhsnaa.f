
      INTEGER FUNCTION JHSNAA( IPU, IERMF, MFERV, KTERV, MXERRS )
C ======================================================================
C     JHSNAA===>KTERMP     P.D. FRANK
C ======================================================================
C
      INTEGER I, IPU
      INTEGER IERMF, MXERRS
      INTEGER MFERV(MXERRS), KTERV(MXERRS)
C
C
      loop: DO I=1, MXERRS
        IF ( MFERV(I) .EQ. IERMF ) THEN
C
C         JHSNAA <-- BARRIER KT ROUTINE ERROR CODE THAT CORRESPONDS
C                    TO THE MULTIFRONTAL ERROR CODE IERMF
C
          JHSNAA = KTERV(I)
          exit loop
        ENDIF
C
        IF ( I .EQ. MXERRS ) THEN
C         NON-RECOVERABLE PROBLEM, INTERNAL ERROR MAP FOR BARKT
C         ROUTINE INCORRECT.
C
          WRITE(IPU,*) ' JHSNAA: NO BARKT ERROR MAP FOR MULTI-FRONTAL',
     $                 ' ERROR NUMBER =', IERMF, ' RUN STOPPED!'
C
          RETURN
        ENDIF
C
      enddo loop
C
C
      RETURN
      END
