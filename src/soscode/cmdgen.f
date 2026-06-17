      subroutine   CMDGEN   ( job0, msglvl,
     1                         n, nclin, nmoved, iter, numinf,
     2                        istate, bl, bu, featol, featlu, x,
     3                        iPrint, iSumm,
     4                        tolinc, idegen, kdegen, ndegen,
     5                        itnfix, nfix )
c     ==================================================================
c     ==================================================================
c     ====  cmdgen / CMDGEN -- degeneracy resolution strategy       ====
c     ====                     initialization                       ====
c     ==================================================================
c     ==================================================================
     
      character(len=1)   job, job0

      integer            msglvl, n     , nclin , nmoved, iter  , numinf,
     1                   iPrint, iSumm , idegen, kdegen, ndegen, itnfix

      integer            nfix (2) 

      double precision   tolinc
      
      integer            istate (n+nclin)
      
      double precision   x (n)
      
      double precision   bl (n+nclin), bu (n+nclin)
      
      double precision   featol (n+nclin), featlu (n+nclin)

c     ==================================================================
c
c     derived from qpopt version 1.0
c     last modification -- 26-March-1996
c     
c         Cmdgen is based on 
c
c         19-Apr-1988. Original version based on MINOS routine m5dgen.
c         09-Apr-1994. Expand frequency allowed to expand. This allows
c                      small initial values of kdegen.
c         28-Jul-1994. Current version.
c     
c     CMDGEN / cmdgen performs most of the manoeuvres associated with
c     degeneracy.  the degeneracy-resolving strategy operates in the
c     following way.
c
c     Over a cycle of iterations, the feasibility tolerance featol
c     increases slightly (from tolx0 to tolx1 in steps of tolinc).
c     this ensures that all steps taken will be positive.
c
c     After kdegen consecutive iterations, variables within
c     featol of their bounds are set exactly on their bounds and x is
c     recomputed to satisfy the general constraints in the working set.
c     Featol is then reduced to tolx0 for the next cycle of iterations.
c
c     featlu  is the array of user-supplied feasibility tolerances.
c     featol  is the array of current feasibility tolerances.
c
c     If job = 'i', cmdgen initializes the feasibility control
c     parameters:
c
c     (tolx0   is the minimum (scaled) feasibility tolerance.)
c     (tolx1   is the maximum (scaled) feasibility tolerance.)
c     tolinc  is the scaled increment to the current featol.
c     idegen  is the expand frequency. It is the frequency of resetting
c             featol to (scaled) tolx0.
c     (kdegen  (specified by the user) is the initial value of idegen.)
c     ndegen  counts the number of degenerate steps (incremented
c             by cmchzr/hdqpcz).
c     (itnfix  is the last iteration at which a job 'e' or 'o' entry
c             caused an x to be put on a constraint.)
c     nfix(j) counts the number of times a job 'o' entry has
c             caused the variables to be placed on the working set,
c             where j=1 if infeasible, j=2 if feasible.
c
c     tolx0*featlu and tolx1*featlu are both close to the feasibility
c     tolerance featlu specified by the user.  (They must both be less
c     than featlu.)
c
c
c     If job = 'e',  cmdgen has been called after a cycle of kdegen
c     iterations.  Constraints in the working set are examined to see if
c     any are off their bounds by an amount approaching featol.  Nmoved
c     returns how many.  If nmoved is positive,  x  is moved onto the
c     constraints in the working set.  It is assumed that the calling
c     routine will then continue iterations.
c
c
c     If job = 'o',  cmdgen is being called after a subproblem has been
c     judged optimal, infeasible or unbounded.  Constraint violations
c     are examined as above.
c
c     ==================================================================
     
      integer            is, j, maxfix
      
      double precision   d

      double precision   hdmcon

      external           hdmcon
      
      double precision   zero, point6, tolx0, tolx1
      
      parameter        ( zero   =  0.0d0,
     1                   point6 =  6.0d-1, 
     2                   tolx0  =  5.0d-1,
     3                   tolx1  = 99.0d-2 )

c     ==================================================================
      job = job0(1:1)
     
      nmoved = 0
     
      if  ( job .eq. 'i'  .or. job .eq. 'I' )  then
         
c        ---------------------------------------------------------------
c        Job = 'Initialize'.
c        Initialize at the start of each linear problem.
c        kdegen  is the expand frequency      and
c        featlu  are the user-supplied feasibility tolerances.
c        They are not changed.
c        ---------------------------------------------------------------

         ndegen   = 0
         itnfix   = 0
         nfix (1) = 0
         nfix (2) = 0

         idegen = kdegen
         if  ( kdegen .lt. 9 999 999 )  then
            tolinc = (tolx1 - tolx0) / idegen
         else
            tolinc = zero
         end if

         do j = 1, n+nclin
            featol(j) = tolx0 * featlu(j)
         enddo

      else
         
c        ---------------------------------------------------------------
c        Job = 'end of cycle' or 'optimal'.
c        Initialize local variables maxfix.
c        ---------------------------------------------------------------
         maxfix = 2

         if  ( job .eq. 'o'  .or. job .eq. 'O' )  then
            
c           ------------------------------------------------------------
c           Job = 'optimal'.
c           Return with nmoved = 0 if the last call was at the same
c           iteration,  or if there have already been maxfix calls with
c           the same state of feasibility.
c           ------------------------------------------------------------

           if  ( itnfix .eq. iter  )  return
           
           if  ( numinf .gt.   0    )  then
               j = 1
            else
               j = 2
            end if

            if  ( nfix (j).ge. maxfix )  then
               return
            endif
            
            nfix (j) = nfix (j) + 1
         end if

c        Increase the expand frequency.
c        Reset featol to its minimum value.

         idegen = idegen + 10
         if  ( kdegen .lt. 9 999 999 )  then
            tolinc = (tolx1 - tolx0) / idegen
            idegen = idegen + iter
         else
            tolinc = zero
         end if

         do j = 1, n+nclin
            featol(j) = tolx0 * featlu(j)
         enddo

c        Count the number of times a variable is moved a nontrivial
c        distance onto its bound.

         itnfix = iter

         do j = 1, n
            
            is     = istate(j)
            if  ( is .gt. 0  .and.  is .lt. 4 )  then            
              if  ( is .eq. 1 )  then
                  d   = abs(x(j) -  bl(j))
               else 
                  d   = abs(x(j) -  bu(j))
               end if

c              ... compare d to eps**.6
               
               if  ( d .gt. hdmcon (5)**point6 )  then
                  nmoved = nmoved + 1
               endif
               
            end if
            
         enddo

        if  ( nmoved .gt. 0  .and.  msglvl .ge. 5 )  then

c           Some variables were moved onto their bounds.

           if  ( iPrint .gt. 0)  write( iPrint, 1000 ) iter, nmoved
           if  ( iSumm  .gt. 0)  write( iSumm , 1000 ) iter, nmoved
         end if
      end if

 1000 format(' Itn', i6, ' --', i7,
     $       '  variables moved to their bounds.')

c     end of CMDGEN (cmdgen)
      end
