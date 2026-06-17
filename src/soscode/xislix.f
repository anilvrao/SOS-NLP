      subroutine xislix ( jrep  , irow  , loc   , neqns , nzcomp, 
     1                    colptr, rowlst )
c
c
c     ==================================================================
c     ====  xislix -- inserts new entry into compressed adjacency   ====
c     ====            structure                                     ====
c     ==================================================================
c     ==================================================================
c
c     purpose
c     -------
c
c     xislix inserts a new entry into the compressed adjacency 
c     structure.
c
c     created         24-jan-97   -- rgg --
c     last modified   
c
c     input arguments
c     ---------------
c
c     jrep        i   column of master compressed column
c     irow        i   row index to be added
c     loc         i   location to insert irow
c     neqns       i   order of uncompressed problem.
c     nzcomp      i   the number of entries in colptr and rowlst
c     colptr,
c       rowlst    i   compressed adjacency structure
c
c     output arguments
c     ----------------
c
c     nzcomp      i   the number of entries in colptr and rowlst
c     colptr,
c       rowlst    i   compressed adjacency structure
c
c     ==================================================================
 
c     --------------
c     ... parameters
c     --------------

      integer             jrep  , irow  , loc   , neqns , nzcomp

      integer             colptr (*)    , rowlst (*)
 
c     -------------------
c     ... local variables
c     -------------------

      integer             k

c     ==================================================================

c     -----------------------------------
c     ... slide everything down in collst
c     -----------------------------------

      do k = nzcomp, loc, - 1
          rowlst(k+1) = rowlst(k)
      enddo

      nzcomp = nzcomp + 1

c     ------------------------------
c     ... insert irow in freed space
c     ------------------------------

      rowlst(loc) = irow

c     -----------------
c     ... adjust colptr
c     -----------------

      do k = jrep+1, neqns+1
          colptr(k) = colptr(k) + 1
      enddo

c     ==================================================================

      return
      end
