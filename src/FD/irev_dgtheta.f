      subroutine irev_dgtheta(ctfm,ctfmp,n_Athth,n_B,n_alpha,n_beta,n_di
     &vbeta,nm1_B,np1_B,Nx,ht,hx,vee,res)
      implicit none
      integer i
      integer Nx
      real*8 ht
      real*8 hx
      real*8 vee
      real*8 ctfm(Nx)
      real*8 ctfmp(Nx)
      real*8 n_Athth(Nx)
      real*8 n_B(Nx)
      real*8 n_alpha(Nx)
      real*8 n_beta(Nx)
      real*8 n_divbeta(Nx)
      real*8 nm1_B(Nx)
      real*8 np1_B(Nx)
      real*8 res(Nx)
      real*8 qb
      do i=2, Nx-1, 1
      qb = 0.5000000000000000D0 * (-0.1D1 * nm1_B(i) + np1_B(i)) / ht + 
     #0.6666666666666667D0 * vee * n_divbeta(i) * n_B(i) + 0.50000000000
     #00000D0 / ctfmp(i) * (n_B(i - 1) - 0.1D1 * n_B(i + 1)) / hx * n_be
     #ta(i) + 0.2D1 * n_Athth(i) * n_alpha(i) - 0.2D1 / ctfm(i) * n_B(i)
     # * n_beta(i)
      res(i)=qb
      end do
      END
