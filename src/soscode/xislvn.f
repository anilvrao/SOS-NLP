      subroutine xislvn (  libver )
 
c     ==================================================================
c     ==================================================================
c     ====  xislvn -- version number for bcslib-ext                 ====
c     ==================================================================
c     ==================================================================
c
c     purpose
c     -------
c
c     returns the version number of bcslib-ext for printing purposes
c     to all of the initialization and print statistics subroutines
c     (x???in and x???ps).
c
c     created         19-oct-99   -- rgg --
c     modified        30-jun-00   -- dkw -- new version 4.0.1
c                     14-dec-00   -- dkw -- new version 4.0.2
c                     15-nov-01   -- dkw -- changed for Release 4.1
c                     18-mar-04   -- dkw -- changed for Release 4.2
c                     21-feb-07   -- dkw -- changed for Release 4.3
c
c     output arguments
c     ----------------
c
c     libver      ch  character(len=55) variable holding current version
c                     number and date of the current version.
c
c     ==================================================================
 
c     --------------
c     ... parameters
c     --------------

      character(len=55)    libver 
 
c     ==================================================================

      libver = 'BCSLIB-EXT Release 4.3  ( 10 Mar 2007 )                '

c     ==================================================================

      return
      end
