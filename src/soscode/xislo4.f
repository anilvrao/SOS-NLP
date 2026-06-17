      subroutine   xislo4   ( neqns,  perm, invp, xrowls, rowlst,
     1                        cmpmap, qsize )
 
c
c     ==================================================================
c     ==================================================================
c     ====  xislo4 -- multiple minimum degree numbering             ====
c     ==================================================================
c     ==================================================================
c
c     created by joseph w. h. liu, york university
c     last modification   mar. 05, 1986
c
c     purpose: this routine performs the final step in
c              producing the permutation and inverse perm vectors
c              in the multiple elimination version of the minimum
c              degree ordering algorithm.
c     input parameters:
c         neqns - number of equations.
c         qsize - size of supernodes at elimination.
c
c     updated parameters:
c         invp - inverse perm vector.  on input,
c                if qsize(node)=0, then node has been merged
c                into the node '-invp(node)';
c                otherwise, -invp(node) is its inverse labelling.
c
c     output parameters:
c         perm - the permutation vector.
c
c     ==================================================================
 
c     --------------
c     ... parameters
c     --------------
 
      integer             neqns
 
      integer             perm (*),   invp (*),   qsize (*),
     1                    xrowls (*), rowlst (*), cmpmap (*)
 
c     -------------------
c     ... local variables
c     -------------------
 
      integer             father, nextf,  node,   nqsize, num,
     1                    root,   k1   ,  k2
 
c     ==================================================================
 
c.debug
c     write(6,'("in xislo4")')
c.debug
      do node = 1, neqns

c.old     nqsize = qsize(node)
c.old     if ( nqsize .le. 0 )  perm(node) = invp(node)
c.old     if ( nqsize .gt. 0 )  perm(node) = - invp(node)

          k1 = cmpmap(node)
          k2 = rowlst(xrowls(k1))

          if ( k2 .eq. node ) then

              nqsize = qsize(k1)

              if ( nqsize .le. 0 ) then
                  perm(node) =   invp(node)
              else
                  perm(node) = - invp(node)
              end if

          else

              perm(node) =   invp(node)

          end if

      enddo
 
c     --------------------------------------------
c     for each node which has been merged, do ....
c     --------------------------------------------
 
      do node = 1, neqns
          if ( perm(node) .gt. 0 )  cycle
 
c         -----------------------------------------
c         trace the merged tree until one which has
c         not been merged, call it root
c         -----------------------------------------
 
          father = node
  200     continue
          if ( perm(father) .le. 0 )  then
              father = - perm(father)
              goto 200
          endif
 
c         -----------------------
c         number node after root.
c         -----------------------
 
          root       = father
          num        = perm(root) + 1
          invp(node) = - num
          perm(root) = num
 
c         ------------------------
c         shorten the merged tree.
c         ------------------------
 
          father = node
 
  400     continue
          nextf = - perm(father)
          if ( nextf .gt. 0 )  then
              perm(father) = - root
              father       = nextf
              goto 400
          endif
 
      enddo
 
c         ----------------------
c         ready to compute perm.
c         ----------------------
 
      do node = 1, neqns
          num        = - invp(node)
          invp(node) = num
          perm(num)  = node
      enddo
 
      return
      end
