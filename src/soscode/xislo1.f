      subroutine   xislo1   ( ncomp,  xrowls, rowlst, cmpmap,
     1                        neqns , xadj  , adjncy, invp  , perm  ,
     1                        nsuper, xsup  , lnzcol, sparnt, delta ,
     2                        dhead , qsize , llist , marker, maxint,
     3                        nofsub, nadj  , neqp1 , xcliqu, clqsiz )
 
c
c     created by joseph w. h. liu, york university
c     last modification   apr. 11, 1986
c     modification        apr. 1988  (bwp)
c     modification        jun. 1988  (bwp)
c     modification        apr. 1993  (cca)
c     modification        dec. 1996  (rgg)
c
c     purpose: this routine implements the minimum degree
c         algorithm.  it makes use of the implicit represent-
c         ation of elimination graphs by quotient graphs, and
c         the notion of indistinguishable nodes.
c         it also implements the modifications by
c         multiple elimination and minimum external degree.
c
c         the first set of modifications increases the efficiency
c         of the code by preventing the traversal of some nodes during
c         degree update.  during the  degree update of node v, the code
c         computes the contribution of neighbors not yet appearing in
c         one of v's elimination cliques with a list length, rather than
c         traversing them.
c
c         the second set of modifications include recording the
c         supernode partition, the pointer array for the in-core
c         multifrontal symbolic factorization, and the supernodal
c         elimination tree parent vector.
c
c         the third and fourth set of modifications were to 
c         convert to using a compressed adjacency structure.  the
c         third set did an inplace compression.  the fourth set
c         assumes a compression done prior to the ordering.
c
c         caution: the adjacency structure will be destroyed.
c
c     input parameters:
c         ncomp  - number of compressed nodes
c         (xrowls, rowlst) - list of rows for each compressed node
c         cmpmap - map from row numbers to compressed nodes
c         neqns  - number of equations.
c         (xadj, adjncy) - the adjacency structure.
c         delta  - tolerance value for multiple elimination.
c         maxint - maximum machine representable integer (any
c                  smaller estimate will do.) for marking nodes.
c
c     output parameters:
c         perm   - the minimum degree ordering.
c         invp   - the inverse of perm.
c         nofsub - an upper bound on the number of nonzero
c                  subscripts for the compressed storage scheme.
c         nsuper - number of supernodes.
c         xsup   - supernode partition
c         lnzcol - number of nonzeros in each column of the cholesky
c                  factor l, excluding the diagonal entry
c         sparnt - parent vector for the supernodal elimination tree
c
c     working parameters:
c         dhead  - vector for head of degree lists.
c         invp   - used temporarily for degree forward link.
c         perm   - used temporarily for degree backward link.
c         qsize  - vector for size of supernodes.
c         llist  - vector for temporary linked lists.
c         marker - a temporary marker vector.
c         xcliqu - pointer to the first clique stored in each
c                  uneliminated node's adjacency set
c         clqsiz - elimination clique sizes
c
c     program subroutines:
c         xisloh, xislo2, xislo6, xislo4, xislo5
c
c     ==================================================================
 
c     --------------
c     ... parameters
c     --------------
 
      integer             ncomp , neqns , delta , maxint, nofsub, 
     .                    nadj  , neqp1 , nsuper
 
      integer             xrowls (*), rowlst (*), cmpmap (*),
     1                    adjncy (*), invp   (*), perm   (*),
     1                    dhead  (*), qsize  (*), llist  (*),
     2                    marker (*), xcliqu (*), clqsiz (*),
     3                    xsup   (*), lnzcol (*), sparnt (*)
 
      integer             xadj (*)
 
c     -------------------
c     ... local variables
c     -------------------
 
      integer             ehead, i, mdeg, mdlmt, mdnode, nextmd, num,
     1                    tag , istart, istop, nabor, isuper, nisola,
     2                    nsupp1, numnzs, fstcol, lstcol, icol, k1

c.debug
c     integer             k2
c.debug
 
c.d   integer             msglvl
c.d   parameter           (msglvl = 4)
 
      external            xisloh, xislo2, xislo4, xislo5, xislo6
 
c     ==================================================================
 
      mdnode = 0
      num    = 1
      nofsub = 0
 
c     ======================================
      nsuper    = 0
 
      do i = 1, neqns
          lnzcol(i) = 0
          sparnt(i) = 0
      enddo
c     ======================================
 
c     ----------------------------------------------
c     initialization for min deg algorithm.
c     num counts the number of ordered nodes plus 1.
c     here is where a 1993/cca modification was made.
c     the previous call was to xislo3.
c     ----------------------------------------------
c
c     call xislo3 ( neqns, xadj, adjncy, dhead, invp, perm,
c    1              qsize, llist, marker, nadj, neqp1, xcliqu,
c    2              clqsiz )
 
      call xisloh( ncomp , xrowls, rowlst, cmpmap,
     1             neqns , xadj,   adjncy, dhead,  invp,   perm,
     2             maxint, qsize , llist,  
     3             marker, xcliqu, clqsiz )
 
c     ------------------------------------------------------
c     eliminate all isolated nodes.
c     note : here is where a modification was made.
c            the code used to assume that all isolated nodes
c            were singletons. this need not be true if the
c            graph is compressed. therefore, the lines
c                nisola = nisola + qsize(mdnode)
c                num = num + 1
c            were replaced by
c                nisola = nisola + qsize(mdnode)
c                num = num + qsize(mdnode)
c     ------------------------------------------------------
 
      nisola = 0
      nextmd = dhead(1)
 
  100     continue
          if ( nextmd .gt. 0 ) then
              mdnode = nextmd
              nextmd = invp(mdnode)
              k1     = cmpmap(mdnode)
c.d           if ( msglvl .gt. 1 ) then
c.d              write(*,*) 'eliminating isolated node ', mdnode
c.d           endif
c.debug
c             k2 = rowlst(xrowls(k1))
c             if ( k2 .ne. mdnode ) then
c                 write(6,'("oops no. 01 - mdnode, k1, k2 = ", 3i8)')
c    1                                      mdnode, k1, k2
c             end if
c.debug
              nisola = nisola + qsize(k1)
 
c             =============================
              nsuper           = nsuper + 1
              xsup(nsuper)     = num
              lnzcol(num)      = qsize(k1) - 1
c             =============================
 
              marker(mdnode) = maxint
              invp(mdnode)   = - num
              num            = num + qsize(k1)
              goto 100
 
          endif
c     ---------------------------------------------
c     search for node of the minimum degree.
c          mdeg is the current min degree;
c          tag is used to facilitate marking nodes.
c     ---------------------------------------------
 
c.d       if ( msglvl .gt. 0 ) then
c.d          write(*,*) nisola, ' isolated nodes eliminated'
c.d       endif
          if ( num .gt. neqns )  goto 800
              tag      = 1
              dhead(1) = 0
              mdeg     = 2
  300         continue
              if ( dhead(mdeg) .le. 0 )  then
                  mdeg = mdeg + 1
                  goto 300
              endif
 
c             -------------------------------------------------
c             use value of delta to set up mdlmt, which governs
c             when a degree update is to be performed.
c             -------------------------------------------------
 
              mdlmt = mdeg + delta
              ehead = 0
  500         continue
              mdnode = dhead(mdeg)
              if ( mdnode .le. 0 ) then
                  mdeg = mdeg + 1
                  if ( mdeg .gt. mdlmt )  goto 700
                      goto 500
 
              endif
c             ----------------------------------------
c             remove mdnode from the degree structure.
c             ----------------------------------------
 
              nextmd      = invp(mdnode)
              dhead(mdeg) = nextmd
              k1          = cmpmap(mdnode)
              if ( nextmd .gt. 0 )  perm(nextmd) = - mdeg
c.d           if ( msglvl .gt. 1 ) then
c.d               write(*,*) 'eliminating node ', mdnode,
c.d  1                       ', mdeg = ', mdeg, ', size = ',
c.d  2                       qsize(k1), ', num = ', num
c.d           endif
 
c             ===========================================
              nsuper           = nsuper + 1
              xsup(nsuper)     = num
              lnzcol(num)      = mdeg + qsize(k1) - 2
c             ===========================================
 
              clqsiz(mdnode) = mdeg - 1
 
              invp(mdnode) = - num
              nofsub       = nofsub + mdeg + qsize(k1) - 2
              if (num+qsize(k1).gt.neqns) goto 800
 
c             ---------------------------------------------
c             eliminate 'mdnode' and perform quotient graph
c             transformation. reset tag value if necessary.
c             ---------------------------------------------
 
              tag = tag + 1
              if ( tag .ge. maxint ) then
                tag = 1
                do i = 1, neqns
                   if ( marker(i) .lt. maxint ) marker(i) = 0
                enddo
              endif
 
              call xislo2 ( mdnode, xadj,   adjncy, dhead,  invp,
     1                      perm,   xrowls, rowlst, cmpmap, qsize, 
     2                      llist,  marker, maxint, tag,    neqns,  
     3                      nadj,   neqp1,  xcliqu, clqsiz, sparnt )
 
              k1 = cmpmap(mdnode)
              num           = num + qsize(k1)
              llist(mdnode) = ehead
              ehead         = mdnode
              if ( delta .ge. 0 )  goto 500
 
c             -----------------------------------------------
c             update degrees of the nodes involved in the min
c             degree nodes elimination.
c             -----------------------------------------------
 
  700         continue
              if ( num .le. neqns ) then
                call xislo6 ( ehead,  neqns,  xadj,   adjncy, delta,
     1                        mdeg,   dhead,  invp,   perm,   xrowls,
     2                        rowlst, cmpmap, qsize,  llist,  marker,
     3                        maxint, tag,    nadj,   neqp1,  xcliqu,
     4                        clqsiz )
                goto 300
              endif
 
  800 continue
 
c     =======================================================
c     ---------------------------------
c     finish off preliminary version of
c     the parent vector
c     ---------------------------------
 
      if ( mdnode .eq. 0 ) go to 875
 
      istart = xadj(mdnode)
      istop  = xadj(mdnode+1)-1
 
      do i = istart, istop
          nabor = adjncy(i)
          if  ( nabor .eq. 0 )  go to 875
          if  ( invp(nabor) .lt. 0 )  then
              sparnt(nabor) = mdnode
          endif
      enddo
 
c     ------------------------------------
c     finish xsup(*) and compute lnzcol(*)
c     ------------------------------------
 
  875 continue
      xsup(nsuper+1) = neqns + 1
 
      do isuper = 1, nsuper
          fstcol = xsup(isuper)
          lstcol = xsup(isuper+1) - 1
          numnzs = lnzcol(fstcol)
          do icol = fstcol+1, lstcol
              numnzs       = numnzs - 1
              lnzcol(icol) = numnzs
          enddo
      enddo
 
c     =======================================================
 
c.debug
c     write(6,'("before xislo4")')
c     call xislp3 ( 'xrowls', ncomp+1, xrowls, 6 )
c     call xislp3 ( 'rowlst', neqns  , rowlst, 6 )
c.debug

      call xislo4 ( neqns, perm, invp, xrowls, rowlst, cmpmap, 
     1              qsize )
 
c     =====================================================
c     ----------------------------------------
c     compute parent vector for the supernodal
c     elimination tree
c     ----------------------------------------
      nsupp1 = nsuper + 1
 
      call xislo5 ( neqns , nsuper, nsupp1, xsup  , invp  ,
     1              sparnt, dhead , qsize  )
 
c     =====================================================
c.debug
c     write(6,'("at end of xislo1")')
c     call xislp3 ( 'dhead ', neqns, dhead , 6 )
c     call xislp3 ( 'invp  ', neqns, invp  , 6 )
c     call xislp3 ( 'qsize ', ncomp, qsize , 6 )
c     call xislp3 ( 'llist ', neqns, llist , 6 )
c     call xislp3 ( 'marker', neqns, marker, 6 )
c     call xislp3 ( 'xcliqu', neqns, xcliqu, 6 )
c     call xislp3 ( 'clqsiz', neqns, clqsiz, 6 )
c.debug
 
 
      return
      end
