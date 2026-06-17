      subroutine   xisls1   ( neqns , neqp1 , nadj  , nsuper, ncomp , 
     1                        xrowls, rowlst, cmpmap, xsup  , perm  ,
     2                        invp  , xadj  , adj   , mxtcol, mxanzf, 
     3                        temp )
 
c
c     created by john lewis
c     last modification   oct. 07, 1987 (bwp)
c                             dec.  3, 1992 (rgg) added computation of
c                                                 mxanzf, the max. no.
c                                                 of nonzeroes in a
c                                                 supernode (front)
c                             feb. 25, 1997 (rgg) added xrowls, rowlst,
c                                                 and cmpmap to compute 
c                                                 mxanzf for compressed 
c                                                 adjancies 
c                                                 (do 900 loop)
c
c     purpose - transform the full adjacency structure of  a
c               into the adjacency structure of the lower triangle
c               of  pap', where  p  is a permutation.
c               note that this is the adjacency of the original
c               matrix, not of the cholesky decomposition.
c
c     input parameters -
c         neqns       - number of equations
c         neqp1       - neqns+1, dimension of xadj
c         nadj        - number of entries in adjacency list,
c                       dimension of adj
c         nsuper      - number of supernodes
c         xrowls      - pointer to rowlst
c         rowlst      - list of permuted row numbers for each
c                       compressed node
c         cmpmap      - map of permuted row numbers to compressed 
c                       nodes.
c         xsup        - supernode partition array
c         (xadj, adj) - array pair containing the adjacency
c                       structure of the graph of the matrix.
c         perm, invp  - the permutation vectors
c
c     output parameters -
c         rowlst      - list of permuted row numbers for each
c                       compressed node
c         cmpmap      - map of permuted row numbers to compressed 
c                       nodes.
c         (xadj, adj) - array pair containing the column-wise adjacency
c                       structure of the graph of the lower triangle
c                       of the permuted matrix.
c         mxtcol      - maximum number of nonzeros in any single
c                       column of the permuted lower triangle
c         mxanzf      - maximum number of org. nonzeroes in any
c                       front
c
c     working vector -
c         temp        - temporary vector of length  neqns
c
c     note - this subroutine could be used to pack the adjacency
c            structure from long integers to short integers, by
c            providing a second copy of the adj parameter typed
c            as a short integer.
c
c     ==================================================================
 
c     --------------
c     ... parameters
c     --------------
 
      integer             neqns , neqp1 , nadj  , nsuper,
     1                    ncomp , mxtcol, mxanzf
 
      integer             xrowls (*), rowlst (*), cmpmap (*),
     1                    xsup   (*), adj    (*), invp   (*),
     2                    perm   (*), xadj   (*), temp   (*)
 
c     -------------------
c     ... local variables
c     -------------------
 
      integer             i     , ilen  , j     , jlen  , jstop ,
     1                    jstart, newi  , newj  , newlen, nextj ,
     2                    oldi  , rk    , top

      integer             cmpnod, count , isuper, k     , kbgn  , 
     1                    kend  , nbr   , nbrsiz, node  , nodbgn, 
     2                    nodend, nodsiz, jcol  , newcol, iold
 
c     ==================================================================

c     -----------------------------
c     ... permute entries in rowlst
c     -----------------------------

      top = xrowls(ncomp+1) - 1
c.debug
c     call xislp3 ( 'in xisls1 before 10 rowlst', top, rowlst, 6 )
c.debug

      do j = 1, top
          jcol      = rowlst(j)
          newcol    = invp(jcol)
          rowlst(j) = newcol
      enddo
c.debug
c     call xislp3 ( 'in xisls1 after 10 rowlst', top, rowlst, 6 )
c.debug

c  ---------------------------------------------------------
c  build new map of permuted row numbers to compressed nodes
c  ---------------------------------------------------------
c.debug
c     call xislp3 ( 'in xisls1 before 20 cmpmap', neqns, cmpmap, 6 )
c.debug
 
      temp(1:neqns) = cmpmap(1:neqns)
 
      do i = 1, neqns
          iold      = perm(i)
          cmpmap(i) = temp(iold)
      enddo
c.debug
c     call xislp3 ( 'in xisls1 after 20 cmpmap', neqns, cmpmap, 6 )
c.debug        

c     -------------------------------------------------------------
c     ... pack the adjacency structure of the lower triangle of the
c         permuted matrix into the back half of the full adjacency
c         structure of the original matrix.
c     -------------------------------------------------------------
 
      top    = xadj (neqns + 1)
      nextj  = xadj (neqns + 1)
      mxtcol = 0
 
      do i = neqns, 1, -1
          newi   = invp (i)
          ilen   = nextj - xadj (i)
          newlen = 0
          if  ( ilen .gt. 0 )  then
 
c             -------------------------------------------------
c             ... save only the entries in the 'newi'-th column
c                 of the lower triangle of the permuted matrix
c             -------------------------------------------------
 
              do rk = 1, ilen
                  nextj = nextj - 1
                  j     = adj (nextj)
                  newj  = invp (j)
                  if  ( newj .gt. newi )  then
                      top       = top - 1
                      adj (top) = newj
                      newlen    = newlen + 1
                  endif
              enddo
          endif
 
          xadj (i)    = top
          temp (newi) = newlen
          mxtcol      = max ( mxtcol, newlen )
      enddo
 
c     --------------------------------------------------------------
c     ... move lower triangle of permuted matrix to the front of the
c         adjacency list, permuting the row order as we go.
c     --------------------------------------------------------------
 
      nextj =  1
 
      do newi = 1, neqns
          jlen        = temp (newi)
          oldi        = perm (newi)
          jstart      = xadj (oldi)
          jstop       = jstart + jlen - 1
          temp (newi) = nextj
          if  ( jlen .gt. 0 )  then
 
              do j = jstart, jstop
                  adj (nextj) = adj (j)
                  nextj       = nextj + 1
              enddo
          endif
      enddo
 
c     ----------------------------------------------------------------
c     ... reset adjacency pointers  xadj  to moved and permuted values
c     ----------------------------------------------------------------
 
      do newi = 1, neqns
          xadj (newi) = temp (newi)
      enddo
 
      xadj (neqns+1) = nextj
 
c     -------------------------------------------------------
c     ... compute maximum number of original nonzeroes in any
c         supernode
c     -------------------------------------------------------
 
      mxanzf = 0
 
      do isuper = 1, nsuper

          count  = 0
 
          nodbgn = xsup(isuper)
          nodend = xsup(isuper+1) - 1

          do node = nodbgn, nodend

              cmpnod = cmpmap(node)
              k      = xrowls(cmpnod)
 
              if ( rowlst(k) .ne. node ) cycle

              nodsiz = xrowls(cmpnod+1) - xrowls(cmpnod)
              count  = count + ( nodsiz * ( nodsiz-1 ) ) / 2
         
              kbgn   = xadj(node)
              kend   = xadj(node+1) - 1

              do k = kbgn, kend
   
                  nbr    = adj(k)
                  cmpnod = cmpmap(nbr)
                  nbrsiz = xrowls(cmpnod+1) - xrowls(cmpnod)

                  count  = count + nodsiz * nbrsiz

              enddo

          enddo

          mxanzf = max ( mxanzf, count ) 
c.debug
c     write(6,'("in xisls1 at end of 800 - count, mxanzf = ", 2i8)')
c    1                                      count, mxanzf 
c.debug
 
      enddo
 
      return
      end
