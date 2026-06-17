      subroutine xislfj ( nsuper, xsup  , xpanel, xlindx,
     1                    mxlfrt, mfront, mxtotf )
c
c
c     compute new maximum sizes for a front after the factorization
c
c---------------------------------------------------------------------
 
      integer    nsuper, mxlfrt, mfront, mxtotf
 
      integer    xsup (*), xpanel(*), xlindx (*)
 
      integer    isuper, ipanel, ipbgn , ipend , lfront,
     1           n1    , n2    , n3
 
c---------------------------------------------------------------------
 
      mxlfrt = 0
      mfront = 0
      mxtotf = 0
 
      do isuper = 1, nsuper

         ipbgn  = xsup(isuper)
         ipend  = xsup(isuper+1) - 1

         n2     = xpanel(ipend+1) - xpanel(ipbgn)
         n3     = xlindx ( isuper+1 ) - xlindx ( isuper )

         do ipanel = ipbgn, ipend
 
             n1     = xpanel(ipanel+1) - xpanel(ipanel)
             n2     = n2 - n1

             lfront = n1 + n2 + n3
 
             mxlfrt = max ( mxlfrt, lfront )
             mfront = max ( mfront, lfront*(lfront+1)/2 )
             mxtotf = max ( mxtotf, (n2+n3)*n1 + n1*(n1-1)/2 )
 
         enddo
 
      enddo
 
c---------------------------------------------------------------------
 
      return
      end
