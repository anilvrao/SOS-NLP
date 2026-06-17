      subroutine xdslob( neqns , lnzcol, nsuper, xsup  , parsup, fstsup,
     1                   brosup, nsons , nodlst, nzla  , msizsn,
     1                   mxnind, nsnind, fnzlf , mfront, mupdmt,
     2                   mxtri , mxrect, mxtotf, fsize , usize , isize ,
     3                   opscnt                                        )
 
c
c  purpose -- to analyze the supernodes and return information
c             which is needed for the optimal traversal and
c             symbolic factorization.
c
c  created            -- 25-jun-87, cca
c  last modifications -- 25-jun-87, cca
c                        13-jun-89, rgg
c                        01-sep-98, rgg, 32 bit integer mods and
c                                        removed nsnnzl(dup. of fnzlf)
c
c  --------------------------------------------------------
c  - note : this routine produces statistics for a relaxed
c  -        supernodal partition where nodes lie in a chain
c  --------------------------------------------------------
c
c  input variables
c
c      neqns  -- number of equations
c      lnzcol -- number of nonzeroes in each column of the factor
c      nsuper -- number of supernodes
c      xsup   -- pointers into the supernode list of nodes
c      parsup -- parent of super node
c      fstsup -- first son of super node
c      brosup -- brother of super node
c      nodlst -- list of nodes in each supernode
c
c  working storage
c
c      nsons  -- number of sons for a supernode
c
c  output variables
c
c      msizsn -- maximum number of nodes in a supernode
c      mxnind -- maximum number of supernodal indices for a supernode,
c                needed for the symbolic factorization
c      nsnind -- total number of supernodal indices
c      fnzlf  -- number of nonzeroes in the lower factor
c      mfront -- maximum size of a frontal matrix
c      mupdmt -- maximum size of a stacked update matrix
c      mxtri  -- maximum number of entries in the triangle
c      mxrect -- maximum number of entries in the rectangle
c      mxtotf -- maximum number of entries in the trapezoid
c      fsize  -- sizes of the frontal matrices for each supernode
c      usize  -- sizes of the update matrices for each supernode
c      isize  -- number of supernodal indices for each supernode
c      opscnt -- vector containing operation count information
c
c  =====================================================================
 
      integer            neqns, lnzcol(*), nsuper, xsup(*),
     1                   parsup(*), fstsup(*), brosup(*),
     2                   nsons (*), nodlst(*)
 
      integer            msizsn, mxnind, nsnind, mfront,
     1                   mupdmt, mxrect, mxtotf, mxtri , fsize(*),
     2                   usize(*), isize(*),   nzla
 
      double precision   fnzlf , opscnt(8)
 
      integer            fstnod, isuper, lfront, lstnod, lupdmt, nnode ,
     1                   nzrect, nztri,  nson  , i     , j
 
      double precision   fdops , fsops , fvdops, fvsops, sdops , ssops ,
     1                   svdops, svsops, fnode , ffront, fupdmt
 
c  =====================================================================
 
c  -------------------------------------------------
c  loop through the supernodes, gathering statistics
c  -------------------------------------------------
 
      do i = 1, nsuper
 
          nson = 0
          j = fstsup(i)
 
   20     continue
          if ( j .ne. 0 ) then
             nson = nson + 1
             j = brosup(j)
             go to 20
          end if
 
          nsons(i) = nson
 
      enddo
 
      msizsn = 0
      mxnind = 0
      nsnind = 0
      fnzlf  = 0.
      mxrect = 0
      mxtotf = 0
      mxtri  = 0
      mfront = 0
      mupdmt = 0
      fsops  = 0.
      fvsops = 0.
      sdops  = neqns
      ssops  = 0.
      svdops = 2.*dble(neqns-nsuper)
      svsops = 2.*dble(neqns)
      fdops  = dble(nzla)
      fvdops = dble(neqns)
      do isuper = 1,nsuper
 
          fstnod = nodlst(xsup(isuper))
          lstnod = nodlst(xsup(isuper+1)-1)
          nnode  = xsup(isuper+1) - xsup(isuper)
          lupdmt = lnzcol(lstnod)
          lfront = lupdmt + nnode
          nson   = nsons(isuper)
 
          nztri  = (nnode*(nnode-1))/2
          mxtri  = max ( mxtri, nztri+nnode )
          nzrect = lupdmt*nnode
          fnzlf  = fnzlf + nztri + nzrect
          mxrect = max ( mxrect, nzrect )
          mxtotf = max ( mxtotf, nzrect+nztri+nnode )
          msizsn = max ( msizsn, nnode  )
          mxnind = max ( mxnind, lfront )
          nsnind = nsnind + lupdmt
          fsize(isuper) = (lfront*(lfront+1))/2
          mfront        = max ( mfront, fsize(isuper) )
          usize(isuper) = (lupdmt*(lupdmt+1))/2 + 3
          mupdmt        = max ( mupdmt, usize(isuper) )
          isize(isuper) = lupdmt

          fnode  = nnode 
          fupdmt = lupdmt
          ffront = lfront
 
          fdops  = fdops  + ( 2*ffront*fnode - fnode**2 )
     1            + ((ffront**2+2)*fnode)
     2            - (fnode*(fnode+1)*(2*ffront+1)/2)
     3            + (fnode*(fnode+1)*(2*fnode+1))/6.
          fsops  = fsops  + 0.5*(fupdmt*(fupdmt+1))
          fvsops = fvsops + fupdmt
          sdops  = sdops  + 2.*(fnode*(fnode-1))
          ssops  = ssops  + 4.*(fnode*fupdmt)
 
      enddo
 
      fvdops = fvdops + fnzlf
      opscnt(1) = fdops
      opscnt(2) = fvdops
      opscnt(3) = fsops
      opscnt(4) = fvsops
      opscnt(5) = sdops
      opscnt(6) = svdops
      opscnt(7) = ssops
      opscnt(8) = svsops
 
c  =====================================================================
 
      return
      end
