
      SUBROUTINE STORCH(IOPT,NONZ,NCOL,JCOL,IWORK,LNJCOL,IER)
C
C ======================================================================
C     STORCH===>storch   J.T. BETTS
C ======================================================================
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C         PURPOSE:
C
C             This utility routine is designed to convert between
c             a sparse triple format and a compressed (row or column) 
c             format.  THE ROUTINE IS CODED FOR COLUMN CONVERSION -- 
C             ROW CONVERSION IS ACHIEVED BY INPUTING NROW INSTEAD OF 
C             NCOL, AND IROW INSTEAD OF JCOL.
C
C         INPUT/OUTPUT:
C
C             IOPT   CONVERSION OPTION CODE
C                    = 1  CONVERT JCOL TO JCOLST (triple to column compressed)
c                    on input jcol is of length nonz, and contains the column 
c                    indices.  on output, jcol is of length ncol+1, and 
c                    contains the column start locations.
C                    = 2  CONVERT JCOLST TO JCOL (column compressed to triple)
c                    on input jcol is of length ncol+1, and contains the column
c                    start locations.  on output jcol is of length nonz, and 
c                    contains the column indices.
C             NONZ   total number of nonzeros, i.e. the NUMBER OF ELEMENTS 
c                    IN THE JCOL ARRAY.
C             NCOL   NUMBER OF COLUMNS
C             JCOL   INTEGER COLUMN INDEX VECTOR (LNJCOL)
C             IWORK  INTEGER WORK ARRAY (NCOL)
C             LNJCOL LENGTH OF JCOL
C             IER    ERROR RETURN CODE
c                    = 0     normal return
c                    = -100  LNJCOL < max(nonz,NCOL+1)
c                    = -200  ncol < 1
c                    = -300  nonz < 1
C
      DIMENSION JCOL(LNJCOL),IWORK(NCOL)
C
      IER = 0
C
      IF(LNJCOL.LT.MAX(NONZ,NCOL+1)) THEN
        IER = -100
        JCOL(1) = JHMCON(1)
        RETURN
      ELSEIF(NCOL.LT.1) THEN
        IER = -200
        JCOL(1) = JHMCON(1)
        RETURN
      ELSEIF(NONZ.LT.1) THEN
        IER = -300
        JCOL(1) = JHMCON(1)
        RETURN
      ENDIF
C
      IF(IOPT.EQ.1) THEN
C
C                    CONVERT JCOL TO JCOLST
C
C          COUNT THE NUMBER OF NONZEROS IN EACH COLUMN
C
        IWORK(1:NCOL) = 0
C
        DO I = 1,NONZ
          II = JCOL(I)
          IWORK(II) = IWORK(II) + 1
        ENDDO
C
        JCOL(1:NONZ) = 0
C
        JCOL(1) = 1
        DO I = 2,NCOL+1
          JCOL(I) = JCOL(I-1) + IWORK(I-1)
        ENDDO
C
      ELSE 
C
C                    CONVERT JCOLST TO JCOL
C
C          COUNT THE NUMBER OF NONZEROS IN EACH COLUMN
C
        IWORK(1:NCOL) = 0
C
        DO I = 1,NCOL
          IWORK(I) = JCOL(I+1) - JCOL(I)
        ENDDO
C
        JCOL(1:NONZ) = 0
C
        J = 1
        DO I = 1,NCOL
          JU = J + IWORK(I) - 1
          IF(IWORK(I).GT.0) JCOL(J:JU) = I
          J = J + IWORK(I)
        ENDDO
C
      ENDIF
C
      RETURN
      END SUBROUTINE STORCH


      SUBROUTINE STORHL(NDIM,MCON,QPOPTN,ALGNAM,NHOLD,NIHOLD)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C         THIS ROUTINE COMPUTES AN ESTIMATE FOR THE LENGTH
C         OF THE REAL AND INTEGER HOLD ARRAYS (NHOLD, NIHOLD)
C         BASED ON AN EMPERICAL FIT TO A SET OF TEST PROBLEMS
C
      CHARACTER(LEN=6) QPOPTN,ALGNAM
C
      DIMENSION HLDDN1(10), HINDN1(10) 
      DIMENSION SLDDN1(10), SINDN1(10) 
      DIMENSION BLDDN1(10), BINDN1(10) 
C
      DATA HLDDN1 / 
     $    0.471045181313D+01, 0.439474626000D+02, 0.230989966512D+01,
     $   -0.790455969992D-02, 0.269142344284D+02, 0.571288112403D-01,
     $   -0.110380479038D-03,-0.228589464354D+00, 0.120243059526D-01,
     $   -0.245960553187D-02 /
      DATA HINDN1 / 
     $    0.225116580315D+03, 0.481850880430D+01,-0.731354222711D+00,
     $    0.757599037483D-01, 0.765065220660D+00, 0.114920158733D+00,
     $   -0.575088541652D-02, 0.715387891493D+00,-0.101568957525D+00,
     $    0.292754115704D-01 /
C
      DATA SLDDN1 / 
     $    0.103688987047D+04,-0.717008607390D+02, 0.570189231399D+02,
     $   -0.217853580462D+01,-0.333688034843D+02, 0.315957114859D+02,
     $   -0.413685959333D+00,-0.761378378284D+01, 0.146615571943D+01,
     $   -0.931771296490D-01 /
      DATA SINDN1 / 
     $    0.251959204346D+03, 0.288930969632D+02, 0.925812321545D+01,
     $   -0.304966612814D+00, 0.478688017497D+02, 0.146457646729D+01,
     $   -0.248707632301D-01, 0.908556151242D+01, 0.232698972443D+00,
     $   -0.332042349688D-01 /
C
      DATA BLDDN1 / 
     $    0.479892553185D+04, -0.721647397662D+03,  0.507241879755D+02,
     $    0.151530020030D+01,  0.318247019097D+03, -0.125308284317D+02,
     $    0.204187622154D+00,  0.541975136971D+02, -0.281776373324D+01,
     $   -0.347120388906D+00  /
      DATA BINDN1 / 
     $    0.264852167772D+03,  0.201410467779D+01,  0.280658075516D+01,
     $   -0.574092936854D-01,  0.157306439491D+02,  0.350693130609D+00,
     $   -0.630833908370D-02,  0.249142948485D+01,  0.632933037273D-01,
     $    0.207566379256D-01  /
C
      IF(ALGNAM.EQ.'SPRNLP') THEN
C
        IF(QPOPTN.EQ.'DENSE ') THEN
C
C         DENSE QP
C
          RNHOLD = SURFUN(HLDDN1,NDIM,MCON)
          NHOLD = NINT(RNHOLD) + 1
C
          RNIHLD = SURFUN(HINDN1,NDIM,MCON)
          NIHOLD = NINT(RNIHLD) + 1
C
        ENDIF
C
        IF(QPOPTN.EQ.'SPARSE') THEN
C
C         SPARSE QP
C
          RNSOLD = SURFUN(SLDDN1,NDIM,MCON)
          NHOLD = NINT(RNSOLD) + 1
C
          RNISLD = SURFUN(SINDN1,NDIM,MCON)
          NIHOLD = NINT(RNISLD) + 1
C
        ENDIF
C
      ELSEIF(ALGNAM.EQ.'BARNLP') THEN
C
C         SPARSE BARRIER
C
        RNBOLD = SURFUN(BLDDN1,NDIM,MCON)
        NHOLD = NINT(RNBOLD) + 1
C
        RNIBLD = SURFUN(BINDN1,NDIM,MCON)
        NIHOLD = NINT(RNIBLD) + 1
C
      ENDIF
C
      RETURN
      END
