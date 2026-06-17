      subroutine   xislq2   (n, key, error)
 
c
c     q u i c k s o r t
c
c         in the style of the cacm paper by bob sedgewick, october 1978
c
c         sorts integers into descending order.
c
c     input:
c         n    -- number of elements to be sorted
c         key  -- an array of length  n  containing the values
c                 which are to be sorted
c
c     output:
c         key  -- will be arranged so that values are in descending
c                 order
c         error -- 0 unless  n <= 0  or  n  is truly gigantic
c                  (signalling an input error).
c
c-----------------------------------------------------------------------
 
c     -------------------------------
c     ... global variable declaration
c     -------------------------------
 
      integer                 error, n
 
      integer                 key (*)
 
c     ------------------------------
c     ... local variable declaration
c     ------------------------------
 
      integer                 stklen, tiny
 
      parameter               ( stklen = 50,
     1                          tiny   =  9 )
 
      integer                 i     , ip1   , j     , jm1   , k     ,
     1                        left  , llen  , right , rlen  , top   ,
     2                        v
 
      logical                 done
 
      integer                 stack (stklen)
 
c-----------------------------------------------------------------------
 
c     ------------------------------------------------------------------
c     ... program is a direct translation into fortran of sedgewick's
c         program 2, which is non-recursive, ignores files of length
c         less than 'tiny' during partitioning, and uses median of three
c         partitioning.
c     ------------------------------------------------------------------
 
      error = 0
 
      if  (n .eq. 1)  return
      if  (n .lt. 1)  go to 6000
 
      top = 1
      left = 1
      right = n
      done = (n .le. tiny)
 
c     ------------------------------------------------------------
c     ... quicksort -- partition the file until no subfile remains
c         of length greater than 'tiny'
c     ------------------------------------------------------------
 
c     -------------------------
c     ... while not done do ...
c     -------------------------
 
  100 continue
      if  (done)  go to 2000
 
c         -------------------------------------------------------------
c         ... find median of left, right and middle elements of current
c             subfile, which is  key(left), ..., key(right)
c             (correction to cacm article incorporated)
c         -------------------------------------------------------------
 
          k = key ((left+right)/2)
          key ((left+right)/2) = key (left)
          key (left) = k
 
          if  ( key(left+1) .lt. key(right) ) then
              k = key (left+1)
              key (left+1) = key (right)
              key (right) = k
          endif
 
          if  ( key(left) .lt. key(right) )  then
              k = key (left)
              key (left) = key (right)
              key (right) = k
          endif
 
          if  ( key (left+1) .lt. key (left) )  then
              k = key (left+1)
              key (left+1) = key (left)
              key (left) = k
          endif
 
          v = key (left)
 
c         -----------------------------------------------------------
c         ... v is now the median value of the three keys.  now move
c             from the left and right ends simultaneously, exchanging
c             keys until all keys greater than  v  are packed to the
c             left, all keys less than  v  are packed to the right.
c         -----------------------------------------------------------
 
          i = left+1
          j = right
 
c         -------------------------------------
c         loop
c             repeat i = i+1 until key(i) <= v;
c             repeat j = j-1 until key(j) >= v;
c         exit if j < i;
c             << exchange keys i and j >>
c         end
c         -------------------------------------
 
  500     continue
  600     continue
              i  = i + 1
              if  ( key(i) .gt. v )  go to 600
 
  700      continue
              j = j - 1
              if  ( key(j) .lt. v )  go to 700
 
          if  (j .ge. i)  then
              k = key (i)
              key (i) = key (j)
              key (j) = k
              go to 500
          endif
 
          k = key (left)
          key (left) = key (j)
          key (j) = k
 
c         --------------------------------------------------------
c         ... we have now partitioned the file into two subfiles,
c             one is (left ... j-1)  and the other is (i...right).
c             process the smaller next.  stack the larger one.
c         --------------------------------------------------------
 
          llen = j-left
          rlen = right - i + 1
          if  ( max (llen, rlen) .gt. tiny )  go to 1100
 
c             -------------------------------------------------------
c             ... both subfiles are tiny, so unstack next larger file
c             -------------------------------------------------------
 
              if  (top .eq. 1)  go to 900
                  top = top - 2
                  left = stack (top)
                  right = stack (top+1)
                  go to 1000
 
  900             continue
                  done = .true.
 
 1000             continue
                  go to 1700
 
c             ---------------------------------------
c             ... else one or both subfiles are large
c             ---------------------------------------
 
 1100     continue
          if  (min (llen, rlen) .le. tiny)  then
 
c             ----------------------------------------------------------
c             ... one subfile is small, one large.  ignore the small one
c             ----------------------------------------------------------
 
              if  ( llen .le. rlen )  then
                  left = i
              else
                  right = j - 1
              endif
              go to 1700
          endif
 
c         ------------------------------------------------------
c         ... else both are larger than tiny.  one must stacked.
c         ------------------------------------------------------
 
          if  ( top+1 .gt. stklen )  go to 6000
          if  ( llen .le. rlen )  then
              stack (top) = i
              stack (top+1) = right
              right = j-1
          else
              stack (top) = left
              stack (top+1) = j-1
              left = i
          endif
 
          top = top + 2
 
 1700 continue
      go to 100
 
c     ---------------------------------------------------------------
c     ... insertion sort the entire file, which consists of a list of
c         'tiny' subfiles, locally out of order, globally in order.
c     ---------------------------------------------------------------
 
c     -------------------------------------------------------
c     ... insertion sort ... for i := n-1 step -1 to 1 do ...
c     -------------------------------------------------------
 
 2000 continue
      i = n - 1
      ip1 = n
 
 2100     continue
          if  ( key (i) .le. key (ip1) )  then
 
c             ---------------------------------------------
c             ... out of order ... move up to correct place
c             ---------------------------------------------
 
              k = key (i)
              j = ip1
              jm1 = i
 
c             ------------------------------------------------
c             ... repeat ... until 'correct place for k found'
c             ------------------------------------------------
 
 2200         continue
                  key (jm1) = key (j)
                  jm1 = j
                  j = j + 1
                  if  ( j .gt. n )  go to 2300
                  if  (key (j) .gt. k)  go to 2200
 
 2300         continue
              key (jm1) = k
 
          endif
          ip1 = i
          i = i - 1
          if  ( i .gt. 0 )  go to 2100
 
      return
 
c-----------------------------------------------------------------------
 
 6000 continue
      error = 1
      return
 
c-----------------------------------------------------------------------
 
      end
