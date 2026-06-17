      subroutine xisloe( neqns , nsuper, xsup  , nodlst, perm  , parsup,
     1                   fstsup, brosup, snroot, nassmb, mstack, mmerge,
     2                   isize , xlindx, perm2 , xsup2 , ncomp , xrowls,
     3                   rowlst   )
 
c
c  purpose -- to perform a post-order traversal of the supernodal
c             elimination forest, obtaining
c               1.) the global nodal permutation vector,
c               2.) the new supernodal partition vector
c               3.) the number of assemblies (sons) of each supernode
c               4.) the pointer vector xlindx for the symbolic
c                   factorization
c
c  created            -- 01-jun-87, cca
c  last modifications -- 01-jun-87, cca
c  last modifications -- june 17, 1988, bwp
c  last modifications -- dec. 18, 1997, rgg, added sorting to perserve order
c                                            within compressed nodes
c
c  input variables
c
c      neqns  -- number of equations
c      nsuper -- number of supernodes
c      xsup   -- pointers into the list of nodes for each supernode
c      nodlst -- the list of nodes for each supernode
c                note : on output, nodlst(node) = node for all node
c      perm   -- the initial permutation vector
c      parsup -- supernodal elimination forest parent vector
c      fstsup -- supernodal elimination forest first son vector
c      brosup -- supernodal elimination forest brother vector
c      snroot -- root of the supernodal elimination forest
c      isize  -- number of supernodal indices for each supernode
c      ncomp  -- number of compressed nodes
c      xrowls -- pointer into rowlst
c      rowlst -- list of rows for each compressed node.
c
c  output variables
c
c      perm   -- new nodal permutation vector
c      xsup   -- new pointers into the list of nodes for each supernode,
c      nassmb -- number of sons (and thus assemblies) for each supernode
c      mstack -- maximum number of stacked update matrices
c      mmerge -- maximum number of merges for the symbolic factorization
c      xlindx -- pointer vector for symbolic factorization
c
c  working variables
c
c      perm2  -- temporary permutation vector
c      xsup2  -- temporary pointer vector
c      temp   -- temporary vector of length neqns
c
c  subprograms called
c
c  =====================================================================
 
      integer   neqns, nsuper, xsup(*), nodlst(*), perm(*),
     1          parsup(*), fstsup(*), brosup(*),
     2          snroot, xlindx(*), isize(*),
     3          ncomp, xrowls(*), rowlst(*)
 
      integer   nassmb(*), mstack, mmerge
 
      integer   perm2(*), xsup2(*)
 
      integer   i     , isuper, jsuper, newnod, nsons  , nstack, oldnod,
     1          son

      integer   error, icomp, j, jbgn, jcol, jend, k, newcol
 
c  =====================================================================
 
c  ----------------------------------------
c  copy perm and xsup into perm2 and xsup2,
c  zero out nassmb
c  ----------------------------------------
 
      perm2(1:neqns) = perm(1:neqns)
      xsup2(1:nsuper+1) = xsup(1:nsuper+1)
      xsup(1:nsuper+1) = 0
      nassmb(1:nsuper) = 0
 
c  --------------------------------------------------------
c  perform the given post-order traversal of the supernodal
c  elimination forest, numbering nodes on the fly, creating
c  the new xsup2 and nassmb vectors for the supernodes
c  --------------------------------------------------------
 
      nstack = 0
      mstack = 1
      mmerge = 0
      jsuper = 1
      newnod = 1
      isuper = snroot
 
      xlindx(1) = 1
 
   10 continue
      son = fstsup(isuper)
      if ( son .gt. 0 ) then
          isuper = son
          go to 10
      endif
 
   20 continue
 
c  -------------------------------------------------------------
c  visiting a supernode, number the nodes in it and update xsup2
c  -------------------------------------------------------------
 
      xlindx(jsuper+1) = xlindx(jsuper) + isize(isuper)
      xsup(jsuper)     = newnod
      do i = xsup2(isuper),xsup2(isuper+1)-1
          oldnod       = nodlst(i)
          perm(newnod) = perm2(oldnod)
          newnod       = newnod + 1
      enddo
 
c  ---------------------------------------
c  determine the number of sons for jsuper
c  ---------------------------------------
 
      nsons = 0
      son   = fstsup(isuper)
   40 continue
      if ( son .gt. 0 ) then
          nsons = nsons + 1
          son   = brosup(son)
          go to 40
      endif
      nassmb(jsuper) = nsons
 
      nstack = nstack - nsons
      if ( parsup(isuper) .ne. 0 ) nstack = nstack + 1
      mstack = max ( mstack, nstack )
      mmerge = max ( mmerge, nsons + xsup2(isuper+1) - xsup2(isuper) )
 
c  ---------------------------
c  increment supernode counter
c  ---------------------------
 
      jsuper = jsuper + 1
 
c  -------------------------------
c  proceed onto the next supernode
c  -------------------------------
 
      if ( brosup(isuper) .ne. 0 ) then
          isuper = brosup(isuper)
          go to 10
      endif
      if ( parsup(isuper) .gt. 0 ) then
          isuper                 = parsup(isuper)
          go to 20
      endif
 
c  -----------------------
c  tie off the end of xsup
c  -----------------------
 
      xsup(nsuper +1) = neqns + 1
 
c  ----------------
c  recompute nodlst
c  ----------------
 
      do i = 1, neqns
          nodlst(i) = i
      enddo

c  ---------------------------------------------------------
c  for each compressed node sort the permutation so that the 
c  order is preserved.
c  note that perm2 and xsup2 are all free temporary vectors
c  of length neqns.
c  first generate inverse perm (old-to-new) into perm2.
c  ---------------------------------------------------------
c.debug
c     call xislp3 ( 'in xisloe after 50 perm', neqns, perm, 6 )
c.debug

      call xislog ( neqns, perm, perm2 )

c  --------------------------------------------------------------
c  for each compressed node extract the old-to-new and new-to-old
c  permutation, sort, and then redistribute
c  --------------------------------------------------------------

      do icomp = 1, ncomp

          jbgn = xrowls(icomp)
          jend = xrowls(icomp+1) - 1
 
          k    = 0

          do j = jbgn, jend
              k             = k + 1
              jcol          = rowlst(j)
              newcol        = perm2(jcol)
              xsup2(k)      = newcol
c.debug
c             if ( perm(newcol) .ne. jcol ) then
c                 write(6,'("oops - jcol, perm(newcol) = ", 2i8)')
c    1                               jcol, perm(newcol)
c             end if
c.debug
          enddo
c.debug
c     write(6,'("before sort in xisloe - icomp, k = ", 2i8)')
c    1                                    icomp, k
c     call xislp3 ( 'subset of old-to-new perm', k, xsup2, 6 )
c.debug

          call xislq1 ( k, xsup2, error )
c.debug
c     write(6,'("after  sort in xisloe - icomp, k = ", 2i8)')
c    1                                    icomp, k
c     call xislp3 ( 'subset of old-to-new perm', k, xsup2, 6 )
c.debug

          k    = 0

          do j = jbgn, jend
              k             = k + 1
              jcol          = rowlst(j)
              newcol        = xsup2(k)  
              perm (newcol) = jcol    
          enddo

      enddo
c.debug
c     call xislp3 ( 'in xisloe after 120 perm', neqns, perm, 6 )
c     call xislog ( neqns, perm, perm2 )
c     call xislp3 ( 'in xisloe after 120 invp', neqns, perm2, 6 )
c     call xislp3 ( 'in xisloe after 120 xrowls', ncomp+1, xrowls, 6 )
c     call xislp3 ( 'in xisloe after 120 rowlst', neqns, rowlst, 6 )
c.debug
 
c  =====================================================================
 
      return
      end
