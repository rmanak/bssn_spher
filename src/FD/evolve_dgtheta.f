      subroutine evolve_dgtheta(ctfm,ctfmp,n_Athth,n_B,n_alpha,n_beta,n_
     &divbeta,np1_Athth,np1_B,np1_alpha,np1_beta,np1_divbeta,x,Nx,ht,hx,
     &myzero,vee,zepsdis,phys_bdy,res)
      implicit none
      integer i
      integer Nx
      real*8 ht
      real*8 hx
      real*8 myzero
      real*8 vee
      real*8 zepsdis
      real*8 ctfm(Nx)
      real*8 ctfmp(Nx)
      real*8 n_Athth(Nx)
      real*8 n_B(Nx)
      real*8 n_alpha(Nx)
      real*8 n_beta(Nx)
      real*8 n_divbeta(Nx)
      real*8 np1_Athth(Nx)
      real*8 np1_B(Nx)
      real*8 np1_alpha(Nx)
      real*8 np1_beta(Nx)
      real*8 np1_divbeta(Nx)
      real*8 x(Nx)
      integer phys_bdy(2)
      real*8 res(Nx)
      real*8 qb
      if (phys_bdy(1) .eq. 1) then
      do i=1, 1, 1
      qb = 0.1D1 - 0.1D1 * myzero * x(i)
      res(i)=qb
      end do
      endif
      if (phys_bdy(1) .eq. 1) then
      do i=2, 2, 1
      qb = np1_B(i) - 0.32D2 / (0.7D1 * zepsdis * hx + 0.64D2 * ht) * hx
     # * ht * (-0.5000000000000000D0 * (0.3D1 * np1_B(i - 1) - 0.4D1 * n
     #p1_B(i) + np1_B(i + 1)) / hx + 0.3125000000000000D-1 * zepsdis / h
     #t * (0.7D1 * np1_B(i) + np1_B(i + 2) - 0.4D1 * np1_B(i + 1) - 0.4D
     #1 * np1_B(i - 1)) + 0.3125000000000000D-1 * zepsdis / ht * (0.7D1 
     #* n_B(i) + n_B(i + 2) - 0.4D1 * n_B(i + 1) - 0.4D1 * n_B(i - 1)))
      res(i)=qb
      end do
      endif
      do i=3, Nx-2, 1
      qb = np1_B(i) - 0.48D2 / (0.16D2 * vee * np1_divbeta(i) * ht * ctf
     #m(i) + 0.9D1 * zepsdis * ctfm(i) - 0.48D2 * np1_beta(i) * ht + 0.4
     #8D2 * ctfm(i)) * ht * ctfm(i) * ((-0.1D1 * n_B(i) + np1_B(i)) / ht
     # + 0.3333333333333333D0 * vee * np1_divbeta(i) * np1_B(i) - 0.2500
     #000000000000D0 / ctfmp(i) * (-0.1D1 * np1_B(i - 1) + np1_B(i + 1))
     # / hx * np1_beta(i) + np1_Athth(i) * np1_alpha(i) - 0.1D1 / ctfm(i
     #) * np1_B(i) * np1_beta(i) + 0.3333333333333333D0 * vee * n_divbet
     #a(i) * n_B(i) - 0.2500000000000000D0 / ctfmp(i) * (-0.1D1 * n_B(i 
     #- 1) + n_B(i + 1)) / hx * n_beta(i) + n_Athth(i) * n_alpha(i) - 0.
     #1D1 / ctfm(i) * n_B(i) * n_beta(i) + 0.3125000000000000D-1 * zepsd
     #is / ht * (0.6D1 * np1_B(i) + np1_B(i + 2) + np1_B(i - 2) - 0.4D1 
     #* np1_B(i + 1) - 0.4D1 * np1_B(i - 1)) + 0.3125000000000000D-1 * z
     #epsdis / ht * (0.6D1 * n_B(i) + n_B(i + 2) + n_B(i - 2) - 0.4D1 * 
     #n_B(i + 1) - 0.4D1 * n_B(i - 1)))
      res(i)=qb
      end do
      do i=Nx-1, Nx-1, 1
      qb = np1_B(i) - 0.3D1 / (vee * np1_divbeta(i) * ht * ctfm(i) - 0.3
     #D1 * np1_beta(i) * ht + 0.3D1 * ctfm(i)) * ht * ctfm(i) * ((-0.1D1
     # * n_B(i) + np1_B(i)) / ht + 0.3333333333333333D0 * vee * np1_divb
     #eta(i) * np1_B(i) - 0.2500000000000000D0 / ctfmp(i) * (-0.1D1 * np
     #1_B(i - 1) + np1_B(i + 1)) / hx * np1_beta(i) + np1_Athth(i) * np1
     #_alpha(i) - 0.1D1 / ctfm(i) * np1_B(i) * np1_beta(i) + 0.333333333
     #3333333D0 * vee * n_divbeta(i) * n_B(i) - 0.2500000000000000D0 / c
     #tfmp(i) * (-0.1D1 * n_B(i - 1) + n_B(i + 1)) / hx * n_beta(i) + n_
     #Athth(i) * n_alpha(i) - 0.1D1 / ctfm(i) * n_B(i) * n_beta(i))
      res(i)=qb
      end do
      if (phys_bdy(2) .eq. 1) then
      do i=Nx, Nx, 1
      qb = np1_B(i) - 0.1D1 * ht * ((-0.1D1 * n_B(i) + np1_B(i)) / ht - 
     #0.1D1 * myzero * x(i))
      res(i)=qb
      end do
      endif
      END
