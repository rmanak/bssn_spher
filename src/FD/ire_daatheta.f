      subroutine ire_daatheta(n_DDLthetatheta,n_DDLxx,n_JSthth,n_JSxx,n_
     &Pthetatheta,n_RRthetatheta,n_RRxx,n_Uthetatheta,n_alpha,n_em4phi,n
     &_metthetatheta,n_metxx,nm1_Athth,np1_Athth,Nx,ht,res)
      implicit none
      integer i
      integer Nx
      real*8 ht
      real*8 n_DDLthetatheta(Nx)
      real*8 n_DDLxx(Nx)
      real*8 n_JSthth(Nx)
      real*8 n_JSxx(Nx)
      real*8 n_Pthetatheta(Nx)
      real*8 n_RRthetatheta(Nx)
      real*8 n_RRxx(Nx)
      real*8 n_Uthetatheta(Nx)
      real*8 n_alpha(Nx)
      real*8 n_em4phi(Nx)
      real*8 n_metthetatheta(Nx)
      real*8 n_metxx(Nx)
      real*8 nm1_Athth(Nx)
      real*8 np1_Athth(Nx)
      real*8 res
      real*8 qb
      res = 0.0D0
      do i=1, Nx, 1
      qb = -0.5000000000000000D0 * (nm1_Athth(i) - 0.1D1 * np1_Athth(i))
     # / ht + 0.3333333333333333D0 / n_metxx(i) * n_alpha(i) * n_metthet
     #atheta(i) * n_RRxx(i) * n_em4phi(i) - 0.3333333333333333D0 / n_met
     #xx(i) * n_em4phi(i) * n_alpha(i) * n_JSxx(i) * n_metthetatheta(i) 
     #- 0.3333333333333333D0 * n_alpha(i) * n_RRthetatheta(i) * n_em4phi
     #(i) + 0.3333333333333333D0 * n_alpha(i) * n_JSthth(i) * n_em4phi(i
     #) - 0.1D1 * n_alpha(i) * n_Uthetatheta(i) - 0.3333333333333333D0 /
     # n_metxx(i) * n_DDLxx(i) * n_metthetatheta(i) * n_em4phi(i) + 0.33
     #33333333333333D0 * n_DDLthetatheta(i) * n_em4phi(i) - 0.1D1 * n_Pt
     #hetatheta(i)
      res = res + qb**2
      end do
      res = sqrt(res/(1*Nx))
      END
