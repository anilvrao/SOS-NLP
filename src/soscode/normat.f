
      subroutine NORMAT(rmat,irowr,jstrr,nonzr,nres,hmat,irowh,
     $    jstrh,nonzh,ndim,column,lnres)
c
      implicit double precision (a-h,o-z)
c
c         Purpose:  Construct the Hessian matrix for a sparse
c                   least squares problem, using the normal 
c                   matrix format, i.e. form
c
c            H  = V + (R)^T (R)
c
c                   where R is the residual Jacobian is sparse
c                   column format, and V is the residual Hessian
c                   in sparse column format.   The result,  H,
c                   must have the same sparsity pattern as V.
c
c         Arguments:
c
c           rmat    i    Residual Jacobian nonzeros (nonzr)
c           irowr   i    Row indices for R (nonzr)
c           jstrr   i    Column start locations (ndim+1)
c           nonzr   i    Number of nonzeros in R
c           nres    i    Number of residuals
c           hmat    i/o  Input:  Residual Hessian (nonzh)
c                        Output: Full Hessian (V destroyed)
c           irowh   i    Row indices for H and V (nonzh)
c           jstrh   i    Column start locations (ndim+1)
c           nonzh   i    Number of nonzeroes in H and V
c           ndim    i    Number of variables (columns in H and R)
c           column  i    Dense work array (lnres)
c           lnres   i    Length of column (lnres > nres)
c
      dimension rmat(nonzr), irowr(nonzr) ,jstrr(ndim+1),
     &          hmat(nonzh), irowh(nonzh) ,jstrh(ndim+1),
     &          column(lnres)
c
      parameter (zero=0.d0)
c
c         construct normal matrix elements (irown,jcoln)
c
      do jcoln = 1,ndim
c
c         save dense copy of rmat column jcoln
c
        column(1:nres) = zero
        do jj = jstrr(jcoln),jstrr(jcoln+1)-1
          irr = irowr(jj)
          column(irr) = rmat(jj)
        enddo
c
c         for column jcol loop over the rows in hmat and/or vmat.
c
        do kk = jstrh(jcoln),jstrh(jcoln+1)-1
c
          hsum = zero
c
          irown = irowh(kk)
c
c         sparse loop over column irown of rmat
c
          do ii = jstrr(irown),jstrr(irown+1)-1
            irclm = irowr(ii)
            hsum = hsum + rmat(ii)*column(irclm)
          enddo
c
          hmat(kk) = hmat(kk) + hsum
c
        enddo

      enddo
C
      RETURN
      END
