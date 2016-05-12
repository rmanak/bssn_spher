      subroutine irev_daax(n_DDLthetatheta,n_DDLxx,n_JSthth,n_JSxx,n_Pxx
     &,n_RRthetatheta,n_RRxx,n_Uxx,n_alpha,n_em4phi,n_metthetatheta,n_me
     &txx,nm1_Axx,np1_Axx,Nx,ht,res)
      implicit none
      integer i
      integer Nx
      real*8 ht
      real*8 n_DDLthetatheta(Nx)
      real*8 n_DDLxx(Nx)
      real*8 n_JSthth(Nx)
      real*8 n_JSxx(Nx)
      real*8 n_Pxx(Nx)
      real*8 n_RRthetatheta(Nx)
      real*8 n_RRxx(Nx)
      real*8 n_Uxx(Nx)
      real*8 n_alpha(Nx)
      real*8 n_em4phi(Nx)
      real*8 n_metthetatheta(Nx)
      real*8 n_metxx(Nx)
      real*8 nm1_Axx(Nx)
      real*8 np1_Axx(Nx)
      real*8 res(Nx)
      real*8 qb
      do i=1, Nx, 1
      qb = 0.5000000000000000D0 * (-0.1D1 * nm1_Axx(i) + np1_Axx(i)) / h
     #t - 0.6666666666666667D0 * n_alpha(i) * n_RRxx(i) * n_em4phi(i) + 
     #0.6666666666666667D0 * n_em4phi(i) * n_alpha(i) * n_JSxx(i) + 0.66
     #66666666666667D0 / n_metthetatheta(i) * n_alpha(i) * n_metxx(i) * 
     #n_RRthetatheta(i) * n_em4phi(i) - 0.6666666666666667D0 / n_metthet
     #atheta(i) * n_alpha(i) * n_metxx(i) * n_JSthth(i) * n_em4phi(i) - 
     #0.1D1 * n_alpha(i) * n_Uxx(i) + 0.6666666666666667D0 * n_DDLxx(i) 
     #* n_em4phi(i) - 0.6666666666666667D0 / n_metthetatheta(i) * n_DDLt
     #hetatheta(i) * n_metxx(i) * n_em4phi(i) - 0.1D1 * n_Pxx(i)
      res(i)=qb
      end do
      END
