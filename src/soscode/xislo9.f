      subroutine xislo9( nsuper, parsup, fstsup, brosup, snroot, snsize,
     1                   fellow, nsupr2, xsup  , nodlst, xsup2 , nodls2)
 
c
c  purpose -- to obtain the node list for the relaxed supernode
c             partition.
c
c  created            -- 20-jul-87, cca
c  last modifications -- 20-jul-87, cca
c
c  input variables --
c
c      nsuper -- number of strict supernodes
c      parsup -- parent supernode vector
c      fstsup -- first son supernode vector
c      brosup -- brother supernode vector
c      snsize -- number of nodes in each supernode
c      fellow -- link vector for strict supernodes merged into
c                a relaxed supernode.
c      xsup   -- pointers into the node list for the strict supernodes
c      nodlst -- the node list for the strict supernodes
c      nsupr2 -- number of relaxed supernodes
c
c  output variables --
c
c      xsup2  -- pointers into the node list for the relaxed supernodes
c      nodls2 -- the node list for the relaxed supernodes
c
c  =====================================================================
 
      integer   nsuper, parsup(*), fstsup(*), brosup(*),
     1          snroot, snsize(*), fellow(*), xsup(*),
     2          nodlst(*), nsupr2
 
      integer   xsup2(*), nodls2(*)
 
      integer   i     , isuper, j     , jsuper, ksuper, start
 
c  =====================================================================
 
c  --------------
c  initialization
c  --------------
 
      jsuper   = 0
      start    = 1
      xsup2(1) = 1
 
c  -----------------------------------------------------------
c  perform a post-order traversal of the strict supernode tree
c  -----------------------------------------------------------
 
      isuper = snroot
 
   10 continue
 
      if ( fstsup(isuper) .gt. 0 ) then
          isuper = fstsup(isuper)
          go to 10
      endif
 
   20 continue
 
c  ----------------------------
c  visiting a relaxed supernode
c  ----------------------------
 
      jsuper          = jsuper + 1
      start           = start + snsize(isuper)
      xsup2(jsuper+1) = start
      j               = start - 1
 
c     -----------------------------------------
c     get the nodes of the top strict supernode
c     -----------------------------------------
 
      do i = xsup(isuper+1)-1,xsup(isuper),-1
          nodls2(j) = nodlst(i)
          j         = j - 1
      enddo
 
c     ------------------------------------------------------------
c     loop through the other strict supernodes and get their nodes
c     ------------------------------------------------------------
 
      ksuper = fellow(isuper)
   40 continue
      if ( ksuper .gt. 0 ) then
          do i = xsup(ksuper),xsup(ksuper+1)-1
              nodls2(j) = nodlst(i)
              j         = j - 1
          enddo
          ksuper = fellow(ksuper)
          go to 40
      endif
 
c  ------------------------------------------------
c  proceed to the next strict supernode in the tree
c  ------------------------------------------------
 
      if ( brosup(isuper) .ne. 0 ) then
          isuper = brosup(isuper)
          go to 10
      endif
 
      if ( parsup(isuper) .ne. 0 ) then
          isuper = parsup(isuper)
          go to 20
      endif
 
c  =====================================================================
 
      return
      end
