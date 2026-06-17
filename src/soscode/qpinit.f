      subroutine   QPINIT   ( nerror, msglvl, iPrint, start, 
     1                        liwork, lwork, litotl, lwtotl,
     2                        n, nclin, 
     3                        istate, named, names,
     4                        bigbnd, bl, bu )

c     ==================================================================
c     ==================================================================
c     ====  QPINIT /                                                ====
c     ====  qpinit -- data check and initialization                 ====
c     ==================================================================
c     ==================================================================

      integer            nerror, msglvl, iPrint,
     1                   liwork, lwork, litotl, lwtotl,
     2                   n, nclin

      double precision   bigbnd
      
      character(len=4)   start
      
      character(len=16)  names(*)
      
      logical            named
      
      integer            istate (n+nclin)
      
      double precision   bl (n+nclin), bu (n+nclin)

c     ==================================================================
c     derived from qpopt version 1.0 cminit
c     last modification -- 24-July-1996
c     
c         First version written by PEG,  15-Nov-1990.
c         This version of cminit dated   13-Jul-94.
c     
c     QPINIT / qpinit   checks the data.
c
c     ==================================================================

      integer            j, k ,l

      logical            ok

      double precision   b1, b2

      double precision   zero
      
      parameter        ( zero = 0.0d0 )
                   
      character(len=5)   id (3)
      
      data                id(1)   ,  id(2)   ,  id(3)
     $                 / 'varbl'  , 'lncon'  , 'nlcon'   /
                 
      nerror = 0

c     ==================================================================

c     -----------
c     Check nclin
c     -----------

      if  ( nclin .lt. 0 )  then
         nerror = nerror + 1
         if  ( iPrint .gt. 0) write (iPrint, 1200) nclin
         return
      end if

c     ------------------------------------------------------------------
c     Check that there is enough workspace to solve the problem.
c     ------------------------------------------------------------------

      ok     = litotl .le. liwork  .and.  lwtotl .le. lwork
      
      if  ( .not. ok)  then
         
         nerror = nerror + 1          
         if  ( iPrint .gt. 0 )  then
            write (iPrint, 1100) liwork, lwork, litotl, lwtotl
            write (iPrint, 1110)
         end if
         
      else
         
     1if  ( msglvl .gt. 0)  then
         if  ( iPrint .gt. 0 )  then
            write (iPrint, 1100) liwork, lwork, litotl, lwtotl
         end if
         
      end if
      
c     ------------------------------------------------------------------
c     Check the bounds on all variables and constraints.
c     ------------------------------------------------------------------

      do j = 1, n+nclin
         
         b1     = bl(j)
         b2     = bu(j)      
         ok     = b1 .lt. b2  .or. 
     1            b1 .eq. b2  .and.  abs(b1) .lt. bigbnd

         if  ( .not. ok)  then
            
            nerror = nerror + 1
            if  ( j .gt. n+nclin)  then
               k = j - n - nclin
               l = 3
            else if  ( j .gt. n)  then
               k = j - n
               l = 2
            else
               k = j
               l = 1
            end if

            if  ( iPrint .gt. 0 )  then
               if  ( named )  then
                  if  ( b1 .eq. b2 )  then
                     write (iPrint, 1310) names(j), b1, bigbnd
                  else
                     write (iPrint, 1315) names(j), b1, b2
                  end if
               else 
                  if  ( b1 .eq. b2 )  then
                     write (iPrint, 1300) id(l), k, b1, bigbnd
                  else
                     write (iPrint, 1305) id(l), k, b1, b2
                  end if
               end if
            end if
            
         end if                    
      enddo
         
c     ------------------------------------------------------------
c     Check  istate settings for consistency in warm or hot starts
c     ------------------------------------------------------------

      if  ( start .eq. 'warm'  .or.  start .eq. 'hot ' )  then
         
         do j = 1, n+nclin
            
            ok     = istate (j) .ge. 0   .and.   istate (j) .le. 3
            if  ( .not. ok)  then
               
               nerror = nerror + 1
               if  ( iPrint .gt. 0)   write (iPrint, 1500) j, istate (j)

            else

               if  ( istate (j) .eq. 0 )  then

c                 ... inactive state cannot match equality constraint

                  ok = bl (j) .lt. bu (j)
               else
     1         if  ( istate (j) .eq. 1 )  then

c                 ... fixed on lower bound must have lower bound
                  
                  ok = bl (j) .gt. -bigbnd

               else
     1         if  ( istate (j) .eq. 2 )  then

c                 ... fixed on upper bound must have upper bound
                  
                  ok = bu (j) .lt. bigbnd

               else
     1         if  ( istate (j) .eq. 3  )  then

                  ok = bl (j) .eq. bu (j)

               endif

               if  ( .not. ok)  then
               
                  nerror = nerror + 1
                  if  ( iPrint .gt. 0)   then
                     write (iPrint, 1600) j, istate (j), bl (j), bu (j)
                  endif

               endif
               
            end if
            
         enddo

      end if

      return

 1100 format(/ ' Workspace provided is     iw(', i8,
     1         '),  w(', i8, ').' /
     2         ' To solve problem we need  iw(', i8,
     3         '),  w(', i8, ').')
      
 1110 format(/ ' XXX  Not enough workspace to solve problem.')
      
 1200 format (/' XXX  nclin  is out of range...', i10)
      
 1300 format(/ ' XXX  The equal bounds on  ', a5, i3,
     1         '  are infinite.   Bounds =', g16.7,
     2         '  bigbnd =', g16.7)
      
 1305 format(/ ' XXX  The bounds on  ', a5, i3,
     1     '  are inconsistent.   bl =', g16.7, '   bu =', g16.7)
      
 1310 format(/ ' XXX  The equal bounds on  ', a16,
     1         '  are infinite.   Bounds =', g16.7,
     2         '  bigbnd =', g16.7)
      
 1315 format(/ ' XXX  The bounds on  ', a16,
     1     '  are inconsistent.   bl =', g16.7, '   bu =', g16.7)
      
 1500 format(/ ' XXX  Component', i5, '  of  istate  is out of',
     1         ' range...', i10)

 1600 format(/ ' XXX  Component', i5, 
     1   '  of  istate  is inconsistent, with the bounds.  state = ',
     2   i2,'   bl =', g16.7,
     3   '   bu =', g16.7 )
      
c     end of QPINIT (qpinit)
      
      end
