      subroutine   xdslt2   ( cpbegn, wlbegn, cputim, waltim )
c
c     ==================================================================
c     ====  xdslt2 -- compute elapsed times                         ====
c     ==================================================================
c     ==================================================================
c
c     purpose -- compute elapsed time since beginning of some process
c
c     input parameters --
c
c             cpbegn -- beginning cpu time in some system format
c             wlbegn -- beginning wall clock time in some system format
c
c     updated parameters --
c
c             cputim -- elapsed cpu time in seconds
c             waltim -- elapsed wall clock time in seconds
c
c     created       ... july 12, 1988
c     last modified ... oct. 12, 1988
c
c     ==================================================================
 
c     --------------
c     ... parameters
c     --------------
 
      double precision    cpbegn, wlbegn, cputim, waltim
 
c     -------------------------
c     ... local variables (sun)
c     -------------------------
 
      real                timvec (2), prcstm
 
c     ---------------------------
c     ... subprogram called (sun)
c     ---------------------------
 
      real                etime
 
cgnu      external            etime
 
c     ==================================================================
 
c     -----------------------------------------------------------------
c     ... sun version uses unix utilities  etime and time.
c
c         function  etime  returns three results.  the explicit result
c             of the function is the sum of the following two numbers.
c             timvec (1) is elapsed time exclusive of calls to system
c                    functions (such as i/o or servicing page
c                    faults) at time of call -- that is, cpu time
c             timvec (2) is elapsed time spent servicing system calls
c         in a multiuser environment  timvec(2)  reflects the
c         volatile portion of the net cpu time.  there is no way
c         to separate the paging cost from the user instigated
c         costs.  the resolution is 1/60 of a second
c
c         function  time  returns the wall clock time, in seconds
c         after 00:00 1 jan 1970.  note that the resolution of the
c         wall clock time is in seconds, not 1/60 of a second
c     -----------------------------------------------------------------
 
      prcstm = etime (timvec)
      cputim = timvec(1) - cpbegn
c     waltim = time  ()  - wlbegn
      waltim = 0.
 
      return
      end
