      SUBROUTINE HJRSUN(ISEED,KMAX,N,K)
 
* VERSION FOR REAL MANTISSAS WITH LESS THAN 46 BITS
 
*------------------------
* INTERFACE SPECIFICATION
*------------------------
      INTEGER  ISEED,KMAX,N,K(*)
 
*================================================================
* CODE GENERATED FROM MASTER SOURCE FILE HJRSUN
*
* PURPOSE: STANDARD RANDOM NUMBER GENERATOR FOR UNIFORMLY
*          DISTRIBUTED INTEGERS IN 1..KMAX.  THE OUTPUT
*          WILL BE THE SAME ON ALL COMPUTERS IF KMAX .LE. 32768.
*
*
* REVISION HISTORY:
*
*----------------------------------------------------------------
* PARAMETERS:
*
* KMAX    INPUT    UPPER BOUND OF OUTPUT VALUES.  IF KMAX .LT. 1,
*                  1 IS USED FOR THE BOUND.
*
* N       INPUT    NUMBER OF RANDOM VARIABLES TO BE GENERATED.
*
* ISEED   INOUT    RANDOM SEED.
*                  (USUALLY SET ONCE AND LEFT ALONE.)
*
* K       OUTPUT   ARRAY OF N UNIFORMLY DISTRIBUTED INTEGERS IN
*                  1..MAX(KMAX,1).
*                  (OUTPUT IN K(1:N), NO OUTPUT IF N .LE. 0.)
*
*----------------------------------------------------------------
* METHOD:
*
* THE UNDERLYING INTEGER SEQUENCE IS
*
*      ISEED = MULT*ISEED MOD BASE
*
*      BASE = (2**31)-1 = 2147483647 = PRIME
*      MULT =  7**5     = 16807      = PRIMITIVE ELEMENT OF BASE
*      1 .LE. ISEED .LE. BASE-1
*
* THE OUTPUT IS  K = INT( KMAX*(ISEED/2**31) ) + 1.
*
* REFERENCES:
*
* STEPHEN K. PARK AND KEITH W. MILLER, 'RANDOM NUMBER GENERATORS:
* GOOD ONES ARE HARD TO FIND', COMM. ACM., 31(10), OCTOBER 1988,
* PP. 1192-1201.
*
* DONALD E. KNUTH, SEMINUMERICAL ALGORITHMS - THE ART OF COMPUTER
* PROGRAMMING - VOL. 2 (SECOND EDITION), ADDISON-WESLEY, 1981,
* PP. 114-115.
*
*================================================================
 
      INTEGER           I
      DOUBLE PRECISION  BASE,BASEM,H,MULT,SEED
      PARAMETER ( BASE = 2147483647.0D0, BASEM = 2147483648.0D0,
     A            MULT =      16807.0D0 )
 
      INTRINSIC  INT,MAX,MOD
C
 
*----------------------------------------------------------------
 
      H    = MAX(KMAX,1)
      SEED = ISEED
      DO I = 1, N
         SEED = MOD(MULT*SEED,BASE)
         K(I) = INT( H*(SEED/BASEM) ) + 1
      ENDDO
      ISEED = INT(SEED)
 
      RETURN
      END
