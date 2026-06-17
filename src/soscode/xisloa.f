      subroutine xisloa( neqns , parent, nsuper, xsup  , nodlst, tsldon,
     1                   parsup, fstsup, brosup, snroot, nsntre, sndpth,
     2                   snmson                                        )
 
c
c  purpose -- to determine the supernode elimination tree
c
c  created               03-feb-87, cca
c  last modifications -- 03-feb-87, cca
c                        31-may-87, cca,
c                           modified to allow more general supernode
c                           partition and to return depth and msons
c
c  input variables --
c
c      neqns  -- number of equations
c      parent -- nodal parent vector
c      nsuper -- number of supernodes
c      xsup   -- pointers into the node list for the supernodes
c      nodlst -- list of nodes in each supernode
c
c  working variable --
c
c      tsldon -- inverse of nodlst, gives the supernode for each node
c
c  output variables --
c
c      parsup -- parent vector for the supernodes
c      fstsup -- first son vector for the supernodes
c      brosup -- brother vector for the supernodes
c      snroot -- root of supernode elimination forest
c      nsntre -- number of supernode trees
c      sndpth -- depth of supernode elimination forest
c      snmson -- maximum number of sons in the supernode elimination
c                forest
c
c  subprograms called --
c
c  =====================================================================
 
      integer   neqns, parent(*), nsuper, xsup(*), nodlst(*)
 
      integer   tsldon(*)
 
      integer   parsup(*), fstsup(*), brosup(*),
     1          snroot, nsntre, sndpth, snmson
 
      integer   depth , isuper, j     , jsuper, node  , nsons , parnod,
     1          parsnd, sonsup
 
c  =====================================================================
 
c  -----------------------------------
c  zero out the supernode tree vectors
c  -----------------------------------
 
      do isuper = 1,nsuper
        parsup(isuper) = 0
        fstsup(isuper) = 0
        brosup(isuper) = 0
      enddo
 
c  ----------------------
c  fill the tsldon vector
c  ----------------------
 
      do isuper = 1,nsuper
          do j = xsup(isuper),xsup(isuper+1)-1
              node         = nodlst(j)
              tsldon(node) = isuper
          enddo
      enddo
 
c  ---------------------------------------------
c  loop through the nodes, linking the supernode
c  of each node to the supernode of its parent
c  ---------------------------------------------
 
      do node = 1,neqns
          isuper = tsldon(node)
          parnod = parent(node)
          if ( parnod .ne. 0 ) then
              jsuper = tsldon(parnod)
              if ( jsuper .ne. isuper ) then
                  brosup(isuper) = fstsup(jsuper)
                  fstsup(jsuper) = isuper
                  parsup(isuper) = jsuper
              endif
          endif
      enddo
 
c  ------------------------------------------------------------
c  link together the roots of the supernodal elimination forest
c  ------------------------------------------------------------
 
      snroot = 0
      nsntre = 1
      do isuper = 1,nsuper
          parsnd = parsup(isuper)
          if ( parsnd .eq. 0 ) then
              if ( snroot .eq. 0 ) then
                  snroot = isuper
              else
                  brosup(isuper) = snroot
                  snroot         = isuper
                  nsntre         = nsntre + 1
              endif
          endif
      enddo
 
c  -------------------------------------------------------------------
c  perform a post-order traversal of the supernode elimination forest
c  to determine the depth of the forest and the maximum number of sons
c  -------------------------------------------------------------------
 
      depth  = 1
      sndpth = 1
      snmson = 0
      isuper = snroot
 
   50 continue
 
      if ( fstsup(isuper) .gt. 0 ) then
          isuper = fstsup(isuper)
          depth  = depth + 1
          go to 50
      endif
 
c  --------------------
c  supernode leaf found
c  --------------------
 
      sndpth = max ( sndpth, depth )
 
c  --------------------
c  visiting a supernode
c  --------------------
 
   60 continue
 
      nsons  = 0
      sonsup = fstsup(isuper)
   70 continue
      if ( sonsup .gt. 0 ) then
          nsons  = nsons + 1
          sonsup = brosup(sonsup)
          go to 70
      endif
 
      snmson = max ( snmson, nsons )
 
c  -----------------------------
c  proceed to the next supernode
c  -----------------------------
 
      if ( brosup(isuper) .gt. 0 ) then
          isuper = brosup(isuper)
          go to 50
      endif
 
      if ( parsup(isuper) .eq. 0 ) then
          return
      else
          isuper = parsup(isuper)
          depth  = depth - 1
          go to 60
      endif
 
c  =====================================================================
 
      end
