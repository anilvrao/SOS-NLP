      subroutine   xdslq2   (n, key, datar )
c
c     ==================================================================
c     ==================================================================
c     ====  xdslq3 - sort an integer array into ascending order     ====
c     ====           a real dependent array is also reordered       ====
c     ==================================================================
c     ==================================================================
c
c     created         19-dec-96   -- rgg --
c
c     purpose
c     -------
c
c         sort an integer array into ascending order.  the real   
c         dependent array are reordered in the same manner 
c         as the key array.
c
c     parameters
c     ----------
c
c     input ...
c
c         n       i   number of elements to sort
c
c     updated ...
c
c         key     i   array to be sorted
c         datar   d   dependent data array
c
c     method
c     ------
c
c         quicksort, as modified by bob sedgewick.  see.CMNications
c         of the acm, october, 1978.
c
c         program is a direct translation into fortran of sedgewick's
c         program 2, which is non-recursive, ignores files of length
c         less than 'tiny' during partitioning, and uses median of
c         three partitioning.
c
c     ==================================================================
 
c     --------------
c     ... parameters
c     --------------
 
      integer             n
 
      integer             key (*)

      double precision    datar (*)
 
c     --------------------------------------------------
c     ... hybrid  quicksort / insertion sort  parameters
c     --------------------------------------------------
 
      integer             stklen, tiny
 
      parameter         ( stklen = 200, tiny = 9 )
 
c     -------------------
c     ... local variables
c     -------------------
 
      integer             i     , ip1   , j     , jm1   , left  ,
     1                    middle, llen  , right , rlen  , top
 
      logical             done
 
      integer             k,      v
 
      integer             stack (stklen)

      double precision    dr
 
c     ==================================================================
 
      if  ( n .le. 1 )  go to 3000
 
      top   = 1
      left  = 1
      right = n
      done  = (n .le. tiny)
 
c     -----------------------------------------------------------
c     quicksort -- partition the file until no subfile remains of
c     length greater than 'tiny'
c     -----------------------------------------------------------
 
c     -------------------------
c     ... while not done do ...
c     -------------------------
 
  100 continue
      if  ( .not. done )  then
 
c         ---------------------------------------------------------
c         ... find median of left, right and middle elements of
c             current subfile, which is  key(left), ..., key(right)
c             (correction to cacm article incorporated)
c             median  value is moved to position  left
c             least   value is moved to position  left+1
c             biggest value is moved to position  right
c         ---------------------------------------------------------
 
          middle = ( left + right ) / 2
 
          k              = key (middle)
          dr             = datar (middle)
          key (middle)   = key (left)
          datar (middle) = datar (left)
          key (left)     = k
          datar (left)   = dr
 
          if  ( key(left+1) .gt. key(right) )  then
              k              = key (left+1)
              dr             = datar (left+1)
              key (left+1)   = key (right)
              datar (left+1) = datar (right)
              key (right)    = k
              datar (right)  = dr
          endif
 
          if  ( key(left) .gt. key(right) )  then
              k             = key (left)
              dr            = datar (left)
              key (left)    = key (right)
              datar (left)  = datar (right)
              key (right)   = k
              datar (right) = dr
          endif
 
          if  ( key (left+1) .gt. key (left) )  then
              k              = key (left+1)
              dr             = datar (left+1)
              key (left+1)   = key (left)
              datar (left+1) = datar (left)
              key (left)     = k
              datar (left)   = dr
          endif
 
          v = key (left)
 
c         ------------------------------------------------------------
c         ... v is now the median value of the three keys.  now move
c             from the left and right ends simultaneously, exchanging
c             keys and data until all keys less than  v  are packed to
c             the left, all keys larger than  v  are packed to the
c             right.
c         ------------------------------------------------------------
 
          i = left+1
          j = right
 
c         -------------------------------------
c         loop
c             repeat i = i+1 until key(i) >= v;
c             repeat j = j-1 until key(j) <= v;
c         exit if j < i;
c             << exchange keys i and j >>
c         end
c         -------------------------------------
 
  200     continue
  300     continue
              i  = i + 1
              if  ( key(i) .lt. v )  go to 300
 
  400     continue
              j = j - 1
              if  ( key(j) .gt. v )  go to 400
 
          if  ( j .gt. i )  then
              k         = key (i)
              dr        = datar (i)
              key (i)   = key (j)
              datar (i) = datar (j)
              key (j)   = k
              datar (j) = dr
              go to 200
          endif
 
c         ------------------------------------------
c         ... quicksort  while  loop finished.
c             move median back into proper position.
c         ------------------------------------------
 
          k            = key (left)
          dr           = datar (left)
          key (left)   = key (j)
          datar (left) = datar (j)
          key (j)      = k
          datar (j)    = dr
 
 
c         --------------------------------------------------------
c         ... we have now partitioned the file into two subfiles,
c             one is (left ... j-1)  and the other is (i...right).
c             process the smaller next.  stack the larger one.
c         --------------------------------------------------------
 
          llen = j-left
          rlen = right - i + 1
          if  ( max (llen, rlen) .le. tiny )  then
 
c             -------------------------------------------------------
c             ... both subfiles are tiny, so unstack next larger file
c             -------------------------------------------------------
 
              if  ( top .eq. 1 )  then
                  done = .true.
              else
                  top   = top - 2
                  left  = stack (top)
                  right = stack (top+1)
              endif
 
c             ---------------------------------------
c             ... else one or both subfiles are large
c             ---------------------------------------
 
          else
     1    if  ( min (llen, rlen) .le. tiny )  then
 
c             ------------------------------------
c             ... one subfile is small, one large.
c                 ignore the small one
c             ------------------------------------
 
              if  ( llen .lt. rlen )  then
                  left = i
              else
                  right = j - 1
              endif
 
c         ------------------------------------------------------
c         ... else both are larger than tiny.  one must stacked.
c         ------------------------------------------------------
 
          else
              if  ( llen .lt. rlen )  then
                  stack (top)   = i
                  stack (top+1) = right
                  right         = j-1
              else
                  stack (top)   = left
                  stack (top+1) = j-1
                  left          = i
              endif
 
              top = top + 2
 
          endif
 
          go to 100
 
      endif
 
c     ------------------------------------------------------------
c     insertion sort the entire file, which consists of a list
c     of 'tiny' subfiles, locally out of order, globally in order.
c     ------------------------------------------------------------
 
c     -------------------------------------------------------
c     ... insertion sort ... for i := n-1 step -1 to 1 do ...
c     -------------------------------------------------------
 
      i   = n - 1
      ip1 = n
 
 2100 continue
          if  ( key (i) .gt. key (ip1) )  then
 
c             ---------------------------------------------
c             ... out of order ... move up to correct place
c             ---------------------------------------------
 
              k   = key (i)
              dr  = datar (i)
              j   = ip1
              jm1 = i
 
c             ------------------------------------------------
c             ... repeat ... until 'correct place for k found'
c             ------------------------------------------------
 
 2200         continue
                  key (jm1)   = key (j)
                  datar (jm1) = datar (j)
                  jm1         = j
                  j           = j + 1
 
                  if  ( j .gt. n       )  go to 2300
                  if  ( key (j) .lt. k )  go to 2200
 
 2300         continue
              key (jm1)   = k
              datar (jm1) = dr
 
          endif
 
          ip1 = i
          i   = i - 1
          if  ( i .gt. 0 )  go to 2100
 
 3000 continue
      return
 
      end
