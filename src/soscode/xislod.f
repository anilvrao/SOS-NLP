      subroutine xislod( nsuper, parsup, fstsup, brosup, snroot, fsize ,
     1                   usize , isize , snmson, list  , val   , rstore,
     2                   istore, stkstr, stksti                        )
 
c
c  purpose -- to obtain the optimal supernode forest reorganization
c             with respect to in-core stack storage.
c             this code is based on the report cs-85-02 of joseph liu
c             at york university, with the following changes:
c              1) the supernode tree is reordered
c
c  created            -- 04-feb-87, cca
c  last modifications -- 04-feb-87, cca
c                        18-feb-87, cca, sliding enforced
c                        31-may-87, cca,
c                           massive revisions made to allow for more
c                           general supernode partitions
c                        02-nov-90, mlc, changed list & val dimensioning
c                           to * instead of snmson
c
c  input variables --
c
c      nsuper -- number of supernodes
c      parsup -- supernode parent vector
c      fstsup -- supernode first son vector
c      brosup -- supernode brother vector
c      fsize  -- size of the frontal matrices
c      usize  -- size of the update matrices
c      isize  -- number of supernodal indices for each supernode
c      snmson -- maximum number of sons for a supernode
c
c  working variables --
c
c      list   -- list array to hold supernode numbers of the sons,
c                has dimension snmson
c      val    -- value array to hold comparison values of the sons,
c                has dimension snmson
c      rstore -- working storage array to hold the minimum in-core
c                storage for a node for the subtree rooted at
c                itself.
c      istore -- working storage array to hold the minimum in-core
c                integer storage for a node for the subtree rooted at
c                itself.
c
c  output variables --
c
c      parsup -- updated supernode parent vector
c      fstsup -- updated supernode first son vector
c      brosup -- updated supernode brother vector
c      stkstr -- maximum stack storage
c      stksti -- maximum integer stack storage
c
c  =====================================================================
 
      integer   nsuper, parsup(*), fstsup(*),
     1          brosup(*), snroot, usize(*), fsize(*),
     2          isize(*), snmson
 
      integer   list(*),  val(*), rstore(*),
     1          istore(*)
 
      integer   stksti, stkstr
 
      integer   i     , isum  , isuper, j     , nlist , root  , rsum  ,
     1          sonsup, temp
 
c  =====================================================================
 
c  ------------------------------------------------------------------
c  walk through the supernode tree in the given postordered traversal
c  ------------------------------------------------------------------
 
      isuper = snroot
 
   10 continue
 
      sonsup = fstsup(isuper)
      if ( sonsup .ne. 0 ) then
          isuper = sonsup
          go to 10
      endif
 
   20 continue
 
c  --------------------
c  visiting a supernode
c  --------------------
 
c  ------------------------------------------------
c  get the list of supernode sons of this supernode
c  ------------------------------------------------
 
      nlist  = 0
      sonsup = fstsup(isuper)
 
   30 continue
      if ( sonsup .ne. 0 ) then
          nlist       = nlist + 1
          list(nlist) = sonsup
          val(nlist)  = max ( fsize(isuper), rstore(sonsup) )
     1                              - usize(sonsup)
          sonsup      = brosup(sonsup)
          go to 30
      endif
 
c  ------------------------------
c  if there is more than one son,
c  reorder the sons
c  ------------------------------
 
      if ( nlist .gt. 1 ) then
 
c         ------------------------------
c         if there is more than one son,
c         reorder the sons
c         ------------------------------
 
          do i = 1,nlist-1
              do j = nlist,i+1,-1
                  if ( val(j) .gt. val(j-1) ) then
                      temp      = list(j)
                      list(j)   = list(j-1)
                      list(j-1) = temp
                      temp      = val(j)
                      val(j)    = val(j-1)
                      val(j-1)  = temp
                  endif
              enddo
          enddo
 
c         ----------------------------------
c         the sons are now ordered,
c         redo the fstsup and brosup vectors
c         ----------------------------------
 
          fstsup(isuper) = list(1)
          do i = 2,nlist
              brosup(list(i-1)) = list(i)
          enddo
          brosup(list(nlist)) = 0
 
      endif
 
c  ---------------------------------
c  compute the maximum storage value
c  to appear from the sons
c  ---------------------------------
 
      if ( nlist .eq. 0 ) then
          rstore(isuper) = fsize(isuper)
          istore(isuper) = isize(isuper)
      else
          sonsup         = fstsup(isuper)
          rstore(isuper) = max ( fsize(isuper), rstore(sonsup) )
          rsum           = usize(sonsup)
          istore(isuper) = istore(sonsup)
          isum           = isize(sonsup)
   70     continue
          sonsup         = brosup(sonsup)
          if ( sonsup .gt. 0 ) then
              rstore(isuper) = max ( rstore(isuper),
     1                  max ( fsize(isuper), rstore(sonsup) ) + rsum )
              rsum           = rsum + usize(sonsup)
              istore(isuper) = max ( istore(isuper),
     1                               istore(sonsup) + isum )
              isum           = isum + isize(sonsup)
              go to 70
          endif
      endif
 
c  -----------------------------
c  proceed to the next supernode
c  -----------------------------
 
      if ( brosup(isuper) .ne. 0 ) then
          isuper = brosup(isuper)
          go to 10
      endif
      if ( parsup(isuper) .ne. 0 ) then
          isuper = parsup(isuper)
          go to 20
      endif
 
c  -----------------------------------
c  determine the maximum stack storage
c  -----------------------------------
 
      stkstr = rstore(snroot)
      stksti = istore(snroot)
      root   = brosup(snroot)
   80 continue
      if ( root .ne. 0 ) then
          stkstr = max ( stkstr, rstore(root) )
          stksti = max ( stksti, istore(root) )
          root   = brosup(root)
          go to 80
      endif
 
c  =====================================================================
 
      return
      end
