      subroutine   xislbs   ( index , lenind, lindex, loc )
c
c
c     ==================================================================
c     ====  xislbs -- bisection search of a sorted list for a       ====
c     ====            specific entry.                               ====
c     ==================================================================
c     ==================================================================
c
c     purpose
c     -------
c
c     xislbs performs a bisection search through the array lindex
c     looking for the entry index.  entries in lindex are assumed 
c     to be sorted in ascending order and to be unique.  if loc is
c     nonzero the locations loc+1 and loc-1 are checked first to
c     see if they are the location.
c
c     created         06-jan-98   -- rgg --
c     last modified   
c
c     input arguments
c     ---------------
c
c     index       i   the entry to find 
c     lenind      i   the length of array lindex to be searched
c     lindex      i   the array to be searched.  the entries in lindex 
c                     are assumed to be sorted in ascending order.
c
c     input/output arguments
c     ----------------------
c
c     loc         i   if nonzero on entry, loc represents the last
c                     entry found and is used to start the current
c                     search.  on output, loc is the location of
c                     index in the array lindex.  if index is not
c                     found in lindex, loc is set to 0.
c
c-----------------------------------------------------------------------
 
c     -------------------------------
c     ... global variable declaration
c     -------------------------------
 
      integer                 index , lenind, loc
 
      integer                 lindex (*)
 
c     ------------------------------
c     ... local variable declaration
c     ------------------------------
 
      integer                 fstind, kloc  , lstind
 
c-----------------------------------------------------------------------
c.debug
c     write(6,'("on entry to xislbs - loc, lenind = ", 2i8)')
c    1                                 loc, lenind
c     call xislp3 ( 'lindex', lenind, lindex, 6 )
c.debug

c     ------------------
c     ... initialization
c     ------------------

      if ( lenind .le. 0 ) then
          loc = 0
          return
      end if

      loc = max ( min ( loc, lenind ), 0 ) 

      fstind = 0
      lstind = 0

c     ------------------------
c     ... check location loc+1
c     ------------------------

      kloc = min ( loc+1, lenind )

      if ( index .lt. lindex(kloc) ) then
          lstind = kloc
      else if ( index .eq. lindex(kloc) ) then
          loc = kloc
          return
      else
          fstind = kloc
      end if

c     ------------------------
c     ... check location loc-1
c     ------------------------

      kloc = max ( loc-1, 1 )

      if ( index .lt. lindex(kloc) ) then
          lstind = kloc
      else if ( index .eq. lindex(kloc) ) then
          loc = kloc
          return
      else
          fstind = kloc
      end if

c     -------------------------------------------
c     ... check end points if not already checked
c     -------------------------------------------

      if ( fstind .eq. 0 ) then
          if ( index .lt. lindex(1) ) then
              loc = 0
              return
          else if ( index .eq. lindex(1) ) then
              loc = 1
              return
          else
              fstind = 1
          end if
      end if

      if ( lstind .eq. 0 ) then
          if ( index .lt. lindex(lenind) ) then
              lstind = lenind
          else if ( index .eq. lindex(lenind) ) then
              loc = lenind
              return
          else
              loc = 0
              return
          end if
      end if

c     -------------------------
c     ... bisection search mode
c     -------------------------

  100 continue
      kloc = ( fstind + lstind ) / 2
c.debug
c     write(6,'("index, fstind, lstind, kloc = ", 4i8)')
c    1            index, fstind, lstind, kloc 
c     write(6,'("lindex(fstind)              = ",  i8)')
c    1            lindex(fstind)       
c     write(6,'("lindex(kloc  )              = ",  i8)')
c    1            lindex(kloc  )       
c     write(6,'("lindex(lstind)              = ",  i8)')
c    1            lindex(lstind)       
c.debug

      if ( kloc .eq. fstind ) then
          loc = 0
          return
      end if

      if ( index .lt. lindex(kloc) ) then
          lstind = kloc
      else if ( index .eq. lindex(kloc) ) then
          loc = kloc
          return
      else
          fstind = kloc
      end if

      go to 100
 
c-----------------------------------------------------------------------
 
      end
