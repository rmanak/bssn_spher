      subroutine resid_dgx(ctfmp,n_A,n_Axx,n_alpha,n_beta,n_divbeta,np1_
     &A,np1_Axx,np1_alpha,np1_beta,np1_divbeta,x,Nx,ht,hx,myzero,vee,zep
     &sdis,phys_bdy,res)
      implicit none
      integer i
      integer Nx
      real*8 ht
      real*8 hx
      real*8 myzero
      real*8 vee
      real*8 zepsdis
      real*8 ctfmp(Nx)
      real*8 n_A(Nx)
      real*8 n_Axx(Nx)
      real*8 n_alpha(Nx)
      real*8 n_beta(Nx)
      real*8 n_divbeta(Nx)
      real*8 np1_A(Nx)
      real*8 np1_Axx(Nx)
      real*8 np1_alpha(Nx)
      real*8 np1_beta(Nx)
      real*8 np1_divbeta(Nx)
      real*8 x(Nx)
      integer phys_bdy(2)
      real*8 res
      real*8 qb
      res = 0.0D0
      if (phys_bdy(1) .eq. 1) then
      do i=1, 1, 1
      qb = np1_A(i) - 0.1D1 + myzero * x(i)
      res = res + qb**2
      end do
      endif
      if (phys_bdy(1) .eq. 1) then
      do i=2, 2, 1
      qb = -0.5000000000000000D0 * (0.3D1 * np1_A(i - 1) - 0.4D1 * np1_A
     #(i) + np1_A(i + 1)) / hx + 0.3125000000000000D-1 * zepsdis / ht * 
     #(0.7D1 * np1_A(i) + np1_A(i + 2) - 0.4D1 * np1_A(i + 1) - 0.4D1 * 
     #np1_A(i - 1)) + 0.3125000000000000D-1 * zepsdis / ht * (0.7D1 * n_
     #A(i) + n_A(i + 2) - 0.4D1 * n_A(i + 1) - 0.4D1 * n_A(i - 1))
      res = res + qb**2
      end do
      endif
      do i=3, Nx-2, 1
      qb = -0.1D1 * (n_A(i) - 0.1D1 * np1_A(i)) / ht + np1_Axx(i) * np1_
     #alpha(i) + 0.3333333333333333D0 * np1_A(i) * vee * np1_divbeta(i) 
     #- 0.5000000000000000D0 * ((-0.1D1 * np1_beta(i - 1) + np1_beta(i +
     # 1)) * np1_A(i) + 0.5000000000000000D0 * np1_beta(i) * (-0.1D1 * n
     #p1_A(i - 1) + np1_A(i + 1))) / ctfmp(i) / hx + n_Axx(i) * n_alpha(
     #i) + 0.3333333333333333D0 * n_A(i) * vee * n_divbeta(i) - 0.500000
     #0000000000D0 * ((-0.1D1 * n_beta(i - 1) + n_beta(i + 1)) * n_A(i) 
     #+ 0.5000000000000000D0 * n_beta(i) * (-0.1D1 * n_A(i - 1) + n_A(i 
     #+ 1))) / ctfmp(i) / hx + 0.3125000000000000D-1 * zepsdis / ht * (0
     #.6D1 * np1_A(i) + np1_A(i + 2) + np1_A(i - 2) - 0.4D1 * np1_A(i + 
     #1) - 0.4D1 * np1_A(i - 1)) + 0.3125000000000000D-1 * zepsdis / ht 
     #* (0.6D1 * n_A(i) + n_A(i + 2) + n_A(i - 2) - 0.4D1 * n_A(i + 1) -
     # 0.4D1 * n_A(i - 1))
      res = res + qb**2
      end do
      do i=Nx-1, Nx-1, 1
      qb = -0.1D1 * (n_A(i) - 0.1D1 * np1_A(i)) / ht + np1_Axx(i) * np1_
     #alpha(i) + 0.3333333333333333D0 * np1_A(i) * vee * np1_divbeta(i) 
     #- 0.5000000000000000D0 * ((-0.1D1 * np1_beta(i - 1) + np1_beta(i +
     # 1)) * np1_A(i) + 0.5000000000000000D0 * np1_beta(i) * (-0.1D1 * n
     #p1_A(i - 1) + np1_A(i + 1))) / ctfmp(i) / hx + n_Axx(i) * n_alpha(
     #i) + 0.3333333333333333D0 * n_A(i) * vee * n_divbeta(i) - 0.500000
     #0000000000D0 * ((-0.1D1 * n_beta(i - 1) + n_beta(i + 1)) * n_A(i) 
     #+ 0.5000000000000000D0 * n_beta(i) * (-0.1D1 * n_A(i - 1) + n_A(i 
     #+ 1))) / ctfmp(i) / hx
      res = res + qb**2
      end do
      if (phys_bdy(2) .eq. 1) then
      do i=Nx, Nx, 1
      qb = -0.1D1 * (n_A(i) - 0.1D1 * np1_A(i)) / ht - 0.1D1 * myzero * 
     #x(i)
      res = res + qb**2
      end do
      endif
      res = sqrt(res/(1*Nx))
      END
