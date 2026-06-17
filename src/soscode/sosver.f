      SUBROUTINE SOSVER( CODEV, LIBRV )


C SOSVER returns information about the current version.
C
C This information will be printed on the banner page.
C
C USAGE
C
C    CHARACTER(LEN=60) CODEV, LIBRV
C    CALL SOSVER( CODEV, LIBRV )
C
C CODEV will contain the current version number.
C LIBRV will contain information about the library build.
C
C The contents of the local variable CODEV0 will be updated by the
C person responsible for the source-code maintenance.	 It
C must be changed when the version number changes.
C
C The contents of the local variable LIBRV0 will be updated by the
C person responsible for a library build.  This update is performed
C on the local copy of SOSVER used for the build.
C

      CHARACTER(LEN=60) CODEV, LIBRV, CODEV0, LIBRV0

C     >>>>> Center the character strings about column 30; J. Betts

      DATA CODEV0
     +  /'                     SOS Version 2025.02              '/
C         123456789012345678901234567890123456789012345678901234567890
C                                      ^
      DATA LIBRV0
     +  /'     Compiled with SUSE Linux 10.2 and Intel ifort 10.2'/
C         123456789012345678901234567890123456789012345678901234567890
C                                      ^

C e.g.,  'Compiled with SunOS 3.2 and Fortran 77 4.0' 
C         123456789012345678901234567890123456789012345678901234567890

      CODEV = CODEV0
      LIBRV = LIBRV0

      END
