      subroutine evolve_eqrpsi(n_a1,n_alpha,n_beta,n_rPHI,n_rPI,n_rpsi,n
     &p1_a1,np1_alpha,np1_beta,np1_rPHI,np1_rPI,np1_rpsi,x,Nx,ht,myzero,
     &zepsdis,phys_bdy,res)
      implicit none
      integer i
      integer Nx
      real*8 ht
      real*8 myzero
      real*8 zepsdis
      real*8 n_a1(Nx)
      real*8 n_alpha(Nx)
      real*8 n_beta(Nx)
      real*8 n_rPHI(Nx)
      real*8 n_rPI(Nx)
      real*8 n_rpsi(Nx)
      real*8 np1_a1(Nx)
      real*8 np1_alpha(Nx)
      real*8 np1_beta(Nx)
      real*8 np1_rPHI(Nx)
      real*8 np1_rPI(Nx)
      real*8 np1_rpsi(Nx)
      real*8 x(Nx)
      integer phys_bdy(2)
      real*8 res(Nx)
      real*8 qb
      if (phys_bdy(1) .eq. 1) then
      do i=1, 1, 1
      qb = 0.1333333333333333D1 * np1_rpsi(i + 1) - 0.3333333333333333D0
     # * np1_rpsi(i + 2)
      res(i)=qb
      end do
      endif
      if (phys_bdy(1) .eq. 1) then
      do i=2, 2, 1
      qb = np1_rpsi(i) - 0.32D2 / (0.32D2 + 0.7D1 * zepsdis) * ht * ((-0
     #.1D1 * n_rpsi(i) + np1_rpsi(i)) / ht - 0.5000000000000000D0 * np1_
     #alpha(i) / np1_a1(i) * np1_rPI(i) - 0.5000000000000000D0 * np1_bet
     #a(i) * np1_rPHI(i) - 0.5000000000000000D0 * n_alpha(i) / n_a1(i) *
     # n_rPI(i) - 0.5000000000000000D0 * n_beta(i) * n_rPHI(i) + 0.31250
     #00000000000D-1 * zepsdis / ht * (0.7D1 * np1_rpsi(i) + np1_rpsi(i 
     #+ 2) - 0.4D1 * np1_rpsi(i + 1) - 0.4D1 * np1_rpsi(i - 1)) + 0.3125
     #000000000000D-1 * zepsdis / ht * (0.7D1 * n_rpsi(i) + n_rpsi(i + 2
     #) - 0.4D1 * n_rpsi(i + 1) - 0.4D1 * n_rpsi(i - 1)))
      res(i)=qb
      end do
      endif
      do i=3, Nx-2, 1
      qb = np1_rpsi(i) - 0.16D2 / (0.16D2 + 0.3D1 * zepsdis) * ht * ((-0
     #.1D1 * n_rpsi(i) + np1_rpsi(i)) / ht - 0.5000000000000000D0 * np1_
     #alpha(i) / np1_a1(i) * np1_rPI(i) - 0.5000000000000000D0 * np1_bet
     #a(i) * np1_rPHI(i) - 0.5000000000000000D0 * n_alpha(i) / n_a1(i) *
     # n_rPI(i) - 0.5000000000000000D0 * n_beta(i) * n_rPHI(i) + 0.31250
     #00000000000D-1 * zepsdis / ht * (0.6D1 * np1_rpsi(i) + np1_rpsi(i 
     #+ 2) + np1_rpsi(i - 2) - 0.4D1 * np1_rpsi(i + 1) - 0.4D1 * np1_rps
     #i(i - 1)) + 0.3125000000000000D-1 * zepsdis / ht * (0.6D1 * n_rpsi
     #(i) + n_rpsi(i + 2) + n_rpsi(i - 2) - 0.4D1 * n_rpsi(i + 1) - 0.4D
     #1 * n_rpsi(i - 1)))
      res(i)=qb
      end do
      do i=Nx-1, Nx-1, 1
      qb = np1_rpsi(i) - 0.1D1 * ht * ((-0.1D1 * n_rpsi(i) + np1_rpsi(i)
     #) / ht - 0.5000000000000000D0 * np1_alpha(i) / np1_a1(i) * np1_rPI
     #(i) - 0.5000000000000000D0 * np1_beta(i) * np1_rPHI(i) - 0.5000000
     #000000000D0 * n_alpha(i) / n_a1(i) * n_rPI(i) - 0.5000000000000000
     #D0 * n_beta(i) * n_rPHI(i))
      res(i)=qb
      end do
      if (phys_bdy(2) .eq. 1) then
      do i=Nx, Nx, 1
      qb = np1_rpsi(i) - 0.1D1 * ht * ((-0.1D1 * n_rpsi(i) + np1_rpsi(i)
     #) / ht - 0.1D1 * myzero * x(i))
      res(i)=qb
      end do
      endif
      END
