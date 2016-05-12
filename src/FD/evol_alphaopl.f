      subroutine evol_alphaopl(n_K,n_alpha,np1_K,np1_alpha,x,Nx,gma,gmb,
     &ht,myzero,phys_bdy,res)
      implicit none
      integer i
      integer Nx
      real*8 gma
      real*8 gmb
      real*8 ht
      real*8 myzero
      real*8 n_K(Nx)
      real*8 n_alpha(Nx)
      real*8 np1_K(Nx)
      real*8 np1_alpha(Nx)
      real*8 x(Nx)
      integer phys_bdy(2)
      real*8 res(Nx)
      real*8 qb
      if (phys_bdy(1) .eq. 1) then
      do i=1, 1, 1
      qb = np1_alpha(i) - 0.2D1 / (gma * np1_alpha(i) ** gmb * gmb * np1
     #_K(i) * ht + 0.2D1 * np1_alpha(i)) * ht * np1_alpha(i) * ((-0.1D1 
     #* n_alpha(i) + np1_alpha(i)) / ht + 0.5000000000000000D0 * gma * n
     #p1_alpha(i) ** gmb * np1_K(i) + 0.5000000000000000D0 * gma * n_alp
     #ha(i) ** gmb * n_K(i))
      res(i)=qb
      end do
      endif
      do i=2, Nx-1, 1
      qb = np1_alpha(i) - 0.2D1 / (gma * np1_alpha(i) ** gmb * gmb * np1
     #_K(i) * ht + 0.2D1 * np1_alpha(i)) * ht * np1_alpha(i) * ((-0.1D1 
     #* n_alpha(i) + np1_alpha(i)) / ht + 0.5000000000000000D0 * gma * n
     #p1_alpha(i) ** gmb * np1_K(i) + 0.5000000000000000D0 * gma * n_alp
     #ha(i) ** gmb * n_K(i))
      res(i)=qb
      end do
      if (phys_bdy(2) .eq. 1) then
      do i=Nx, Nx, 1
      qb = -0.1D1 * myzero * x(i) + 0.1D1
      res(i)=qb
      end do
      endif
      END
