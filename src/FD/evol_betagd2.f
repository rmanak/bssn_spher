      subroutine evol_betagd2(n_Lamx,n_beta,np1_Lamx,np1_beta,x,Nx,ck,ep
     &sal,ht,myzero,phys_bdy,res)
      implicit none
      integer i
      integer Nx
      real*8 ck
      real*8 epsal
      real*8 ht
      real*8 myzero
      real*8 n_Lamx(Nx)
      real*8 n_beta(Nx)
      real*8 np1_Lamx(Nx)
      real*8 np1_beta(Nx)
      real*8 x(Nx)
      integer phys_bdy(2)
      real*8 res(Nx)
      real*8 qb
      if (phys_bdy(1) .eq. 1) then
      do i=1, 1, 1
      qb = np1_beta(i) - 0.1D1 * ht * ((-0.1D1 * n_beta(i) + np1_beta(i)
     #) / ht + epsal * (-0.1D1 * n_Lamx(i) + np1_Lamx(i)) / ht + 0.50000
     #00000000000D0 * epsal * ck * np1_Lamx(i) + 0.5000000000000000D0 * 
     #epsal * ck * n_Lamx(i))
      res(i)=qb
      end do
      endif
      do i=2, Nx-1, 1
      qb = np1_beta(i) - 0.1D1 * ht * ((-0.1D1 * n_beta(i) + np1_beta(i)
     #) / ht + epsal * (-0.1D1 * n_Lamx(i) + np1_Lamx(i)) / ht + 0.50000
     #00000000000D0 * epsal * ck * np1_Lamx(i) + 0.5000000000000000D0 * 
     #epsal * ck * n_Lamx(i))
      res(i)=qb
      end do
      if (phys_bdy(2) .eq. 1) then
      do i=Nx, Nx, 1
      qb = -0.1D1 * myzero * x(i)
      res(i)=qb
      end do
      endif
      END
