      subroutine resid_eqrpib2(ctfm,ctfmp,n_a1,n_alpha,n_b1,n_beta,n_rPI
     &b2,n_rpsi,np1_a1,np1_alpha,np1_b1,np1_beta,np1_rPIb2,np1_rpsi,octf
     &mp,x,Nx,ht,hx,mass,myzero,zepsdis,phys_bdy,res)
      implicit none
      integer i
      integer Nx
      real*8 ht
      real*8 hx
      real*8 mass
      real*8 myzero
      real*8 zepsdis
      real*8 ctfm(Nx)
      real*8 ctfmp(Nx)
      real*8 n_a1(Nx)
      real*8 n_alpha(Nx)
      real*8 n_b1(Nx)
      real*8 n_beta(Nx)
      real*8 n_rPIb2(Nx)
      real*8 n_rpsi(Nx)
      real*8 np1_a1(Nx)
      real*8 np1_alpha(Nx)
      real*8 np1_b1(Nx)
      real*8 np1_beta(Nx)
      real*8 np1_rPIb2(Nx)
      real*8 np1_rpsi(Nx)
      real*8 octfmp(Nx)
      real*8 x(Nx)
      integer phys_bdy(2)
      real*8 res
      real*8 qb
      res = 0.0D0
      if (phys_bdy(1) .eq. 1) then
      do i=1, 1, 1
      qb = np1_rPIb2(i) - 0.1333333333333333D1 * np1_rPIb2(i + 1) + 0.33
     #33333333333333D0 * np1_rPIb2(i + 2)
      res = res + qb**2
      end do
      endif
      if (phys_bdy(1) .eq. 1) then
      do i=2, 2, 1
      qb = (-0.1D1 * n_rPIb2(i) + np1_rPIb2(i)) / ht + 0.500000000000000
     #0D0 * np1_a1(i) * mass * np1_rpsi(i) * np1_alpha(i) * np1_b1(i) **
     # 2 - 0.1D1 * np1_rPIb2(i) * np1_beta(i) / ctfm(i) - 0.500000000000
     #0000D0 * ((0.5000000000000000D0 * np1_rPIb2(i) * (-0.1D1 * np1_bet
     #a(i - 1) + np1_beta(i + 1)) + 0.5000000000000000D0 * np1_beta(i) *
     # (-0.1D1 * np1_rPIb2(i - 1) + np1_rPIb2(i + 1)) - 0.50000000000000
     #00D0 / np1_a1(i) * (-0.1D1 * np1_rpsi(i - 1) + np1_rpsi(i + 1)) * 
     #octfmp(i) * np1_alpha(i) * np1_b1(i) ** 2) / ctfmp(i) + 0.1D1 / np
     #1_a1(i) * np1_alpha(i) * np1_b1(i) ** 2 * (-0.1D1 * np1_rpsi(i - 1
     #) + np1_rpsi(i + 1)) / ctfmp(i) / ctfm(i)) / hx - 0.50000000000000
     #00D0 * (-0.2500000000000000D0 / np1_a1(i) ** 2 * (-0.1D1 * np1_rps
     #i(i - 1) + np1_rpsi(i + 1)) * (-0.1D1 * np1_a1(i - 1) + np1_a1(i +
     # 1)) * np1_alpha(i) * np1_b1(i) ** 2 - 0.2500000000000000D0 / np1_
     #a1(i) * (-0.1D1 * np1_rpsi(i - 1) + np1_rpsi(i + 1)) * (np1_alpha(
     #i - 1) - 0.1D1 * np1_alpha(i + 1)) * np1_b1(i) ** 2 + 0.5000000000
     #000000D0 / np1_a1(i) * (-0.1D1 * np1_rpsi(i - 1) + np1_rpsi(i + 1)
     #) * (-0.1D1 * np1_b1(i - 1) + np1_b1(i + 1)) * np1_alpha(i) * np1_
     #b1(i) + 0.1D1 / np1_a1(i) * (np1_rpsi(i - 1) - 0.2D1 * np1_rpsi(i)
     # + np1_rpsi(i + 1)) * np1_alpha(i) * np1_b1(i) ** 2) / ctfmp(i) **
     # 2 / hx ** 2 + 0.5000000000000000D0 * n_a1(i) * mass * n_rpsi(i) *
     # n_alpha(i) * n_b1(i) ** 2 - 0.1D1 * n_rPIb2(i) * n_beta(i) / ctfm
     #(i) - 0.5000000000000000D0 * ((0.5000000000000000D0 * n_rPIb2(i) *
     # (-0.1D1 * n_beta(i - 1) + n_beta(i + 1)) + 0.5000000000000000D0 *
     # n_beta(i) * (-0.1D1 * n_rPIb2(i - 1) + n_rPIb2(i + 1)) - 0.500000
     #0000000000D0 / n_a1(i) * (-0.1D1 * n_rpsi(i - 1) + n_rpsi(i + 1)) 
     #* octfmp(i) * n_alpha(i) * n_b1(i) ** 2) / ctfmp(i) + 0.1D1 / n_a1
     #(i) * n_alpha(i) * n_b1(i) ** 2 * (-0.1D1 * n_rpsi(i - 1) + n_rpsi
     #(i + 1)) / ctfmp(i) / ctfm(i)) / hx - 0.5000000000000000D0 * (-0.2
     #500000000000000D0 / n_a1(i) ** 2 * (-0.1D1 * n_rpsi(i - 1) + n_rps
     #i(i + 1)) * (-0.1D1 * n_a1(i - 1) + n_a1(i + 1)) * n_alpha(i) * n_
     #b1(i) ** 2 - 0.2500000000000000D0 / n_a1(i) * (-0.1D1 * n_rpsi(i -
     # 1) + n_rpsi(i + 1)) * (n_alpha(i - 1) - 0.1D1 * n_alpha(i + 1)) *
     # n_b1(i) ** 2 + 0.5000000000000000D0 / n_a1(i) * (-0.1D1 * n_rpsi(
     #i - 1) + n_rpsi(i + 1)) * (-0.1D1 * n_b1(i - 1) + n_b1(i + 1)) * n
     #_alpha(i) * n_b1(i) + 0.1D1 / n_a1(i) * (n_rpsi(i - 1) - 0.2D1 * n
     #_rpsi(i) + n_rpsi(i + 1)) * n_alpha(i) * n_b1(i) ** 2) / ctfmp(i) 
     #** 2 / hx ** 2 + 0.3125000000000000D-1 * zepsdis / ht * (0.7D1 * n
     #p1_rPIb2(i) + np1_rPIb2(i + 2) - 0.4D1 * np1_rPIb2(i + 1) - 0.4D1 
     #* np1_rPIb2(i - 1)) + 0.3125000000000000D-1 * zepsdis / ht * (0.7D
     #1 * n_rPIb2(i) + n_rPIb2(i + 2) - 0.4D1 * n_rPIb2(i + 1) - 0.4D1 *
     # n_rPIb2(i - 1))
      res = res + qb**2
      end do
      endif
      do i=3, Nx-2, 1
      qb = (-0.1D1 * n_rPIb2(i) + np1_rPIb2(i)) / ht + 0.500000000000000
     #0D0 * np1_a1(i) * mass * np1_rpsi(i) * np1_alpha(i) * np1_b1(i) **
     # 2 - 0.1D1 * np1_rPIb2(i) * np1_beta(i) / ctfm(i) - 0.500000000000
     #0000D0 * ((0.5000000000000000D0 * np1_rPIb2(i) * (-0.1D1 * np1_bet
     #a(i - 1) + np1_beta(i + 1)) + 0.5000000000000000D0 * np1_beta(i) *
     # (-0.1D1 * np1_rPIb2(i - 1) + np1_rPIb2(i + 1)) - 0.50000000000000
     #00D0 / np1_a1(i) * (-0.1D1 * np1_rpsi(i - 1) + np1_rpsi(i + 1)) * 
     #octfmp(i) * np1_alpha(i) * np1_b1(i) ** 2) / ctfmp(i) + 0.1D1 / np
     #1_a1(i) * np1_alpha(i) * np1_b1(i) ** 2 * (-0.1D1 * np1_rpsi(i - 1
     #) + np1_rpsi(i + 1)) / ctfmp(i) / ctfm(i)) / hx - 0.50000000000000
     #00D0 * (-0.2500000000000000D0 / np1_a1(i) ** 2 * (-0.1D1 * np1_rps
     #i(i - 1) + np1_rpsi(i + 1)) * (-0.1D1 * np1_a1(i - 1) + np1_a1(i +
     # 1)) * np1_alpha(i) * np1_b1(i) ** 2 - 0.2500000000000000D0 / np1_
     #a1(i) * (-0.1D1 * np1_rpsi(i - 1) + np1_rpsi(i + 1)) * (np1_alpha(
     #i - 1) - 0.1D1 * np1_alpha(i + 1)) * np1_b1(i) ** 2 + 0.5000000000
     #000000D0 / np1_a1(i) * (-0.1D1 * np1_rpsi(i - 1) + np1_rpsi(i + 1)
     #) * (-0.1D1 * np1_b1(i - 1) + np1_b1(i + 1)) * np1_alpha(i) * np1_
     #b1(i) + 0.1D1 / np1_a1(i) * (np1_rpsi(i - 1) - 0.2D1 * np1_rpsi(i)
     # + np1_rpsi(i + 1)) * np1_alpha(i) * np1_b1(i) ** 2) / ctfmp(i) **
     # 2 / hx ** 2 + 0.5000000000000000D0 * n_a1(i) * mass * n_rpsi(i) *
     # n_alpha(i) * n_b1(i) ** 2 - 0.1D1 * n_rPIb2(i) * n_beta(i) / ctfm
     #(i) - 0.5000000000000000D0 * ((0.5000000000000000D0 * n_rPIb2(i) *
     # (-0.1D1 * n_beta(i - 1) + n_beta(i + 1)) + 0.5000000000000000D0 *
     # n_beta(i) * (-0.1D1 * n_rPIb2(i - 1) + n_rPIb2(i + 1)) - 0.500000
     #0000000000D0 / n_a1(i) * (-0.1D1 * n_rpsi(i - 1) + n_rpsi(i + 1)) 
     #* octfmp(i) * n_alpha(i) * n_b1(i) ** 2) / ctfmp(i) + 0.1D1 / n_a1
     #(i) * n_alpha(i) * n_b1(i) ** 2 * (-0.1D1 * n_rpsi(i - 1) + n_rpsi
     #(i + 1)) / ctfmp(i) / ctfm(i)) / hx - 0.5000000000000000D0 * (-0.2
     #500000000000000D0 / n_a1(i) ** 2 * (-0.1D1 * n_rpsi(i - 1) + n_rps
     #i(i + 1)) * (-0.1D1 * n_a1(i - 1) + n_a1(i + 1)) * n_alpha(i) * n_
     #b1(i) ** 2 - 0.2500000000000000D0 / n_a1(i) * (-0.1D1 * n_rpsi(i -
     # 1) + n_rpsi(i + 1)) * (n_alpha(i - 1) - 0.1D1 * n_alpha(i + 1)) *
     # n_b1(i) ** 2 + 0.5000000000000000D0 / n_a1(i) * (-0.1D1 * n_rpsi(
     #i - 1) + n_rpsi(i + 1)) * (-0.1D1 * n_b1(i - 1) + n_b1(i + 1)) * n
     #_alpha(i) * n_b1(i) + 0.1D1 / n_a1(i) * (n_rpsi(i - 1) - 0.2D1 * n
     #_rpsi(i) + n_rpsi(i + 1)) * n_alpha(i) * n_b1(i) ** 2) / ctfmp(i) 
     #** 2 / hx ** 2 + 0.3125000000000000D-1 * zepsdis / ht * (0.6D1 * n
     #p1_rPIb2(i) + np1_rPIb2(i + 2) + np1_rPIb2(i - 2) - 0.4D1 * np1_rP
     #Ib2(i + 1) - 0.4D1 * np1_rPIb2(i - 1)) + 0.3125000000000000D-1 * z
     #epsdis / ht * (0.6D1 * n_rPIb2(i) + n_rPIb2(i + 2) + n_rPIb2(i - 2
     #) - 0.4D1 * n_rPIb2(i + 1) - 0.4D1 * n_rPIb2(i - 1))
      res = res + qb**2
      end do
      do i=Nx-1, Nx-1, 1
      qb = (-0.1D1 * n_rPIb2(i) + np1_rPIb2(i)) / ht + 0.500000000000000
     #0D0 * np1_a1(i) * mass * np1_rpsi(i) * np1_alpha(i) * np1_b1(i) **
     # 2 - 0.1D1 * np1_rPIb2(i) * np1_beta(i) / ctfm(i) - 0.500000000000
     #0000D0 * ((0.5000000000000000D0 * np1_rPIb2(i) * (-0.1D1 * np1_bet
     #a(i - 1) + np1_beta(i + 1)) + 0.5000000000000000D0 * np1_beta(i) *
     # (-0.1D1 * np1_rPIb2(i - 1) + np1_rPIb2(i + 1)) - 0.50000000000000
     #00D0 / np1_a1(i) * (-0.1D1 * np1_rpsi(i - 1) + np1_rpsi(i + 1)) * 
     #octfmp(i) * np1_alpha(i) * np1_b1(i) ** 2) / ctfmp(i) + 0.1D1 / np
     #1_a1(i) * np1_alpha(i) * np1_b1(i) ** 2 * (-0.1D1 * np1_rpsi(i - 1
     #) + np1_rpsi(i + 1)) / ctfmp(i) / ctfm(i)) / hx - 0.50000000000000
     #00D0 * (-0.2500000000000000D0 / np1_a1(i) ** 2 * (-0.1D1 * np1_rps
     #i(i - 1) + np1_rpsi(i + 1)) * (-0.1D1 * np1_a1(i - 1) + np1_a1(i +
     # 1)) * np1_alpha(i) * np1_b1(i) ** 2 - 0.2500000000000000D0 / np1_
     #a1(i) * (-0.1D1 * np1_rpsi(i - 1) + np1_rpsi(i + 1)) * (np1_alpha(
     #i - 1) - 0.1D1 * np1_alpha(i + 1)) * np1_b1(i) ** 2 + 0.5000000000
     #000000D0 / np1_a1(i) * (-0.1D1 * np1_rpsi(i - 1) + np1_rpsi(i + 1)
     #) * (-0.1D1 * np1_b1(i - 1) + np1_b1(i + 1)) * np1_alpha(i) * np1_
     #b1(i) + 0.1D1 / np1_a1(i) * (np1_rpsi(i - 1) - 0.2D1 * np1_rpsi(i)
     # + np1_rpsi(i + 1)) * np1_alpha(i) * np1_b1(i) ** 2) / ctfmp(i) **
     # 2 / hx ** 2 + 0.5000000000000000D0 * n_a1(i) * mass * n_rpsi(i) *
     # n_alpha(i) * n_b1(i) ** 2 - 0.1D1 * n_rPIb2(i) * n_beta(i) / ctfm
     #(i) - 0.5000000000000000D0 * ((0.5000000000000000D0 * n_rPIb2(i) *
     # (-0.1D1 * n_beta(i - 1) + n_beta(i + 1)) + 0.5000000000000000D0 *
     # n_beta(i) * (-0.1D1 * n_rPIb2(i - 1) + n_rPIb2(i + 1)) - 0.500000
     #0000000000D0 / n_a1(i) * (-0.1D1 * n_rpsi(i - 1) + n_rpsi(i + 1)) 
     #* octfmp(i) * n_alpha(i) * n_b1(i) ** 2) / ctfmp(i) + 0.1D1 / n_a1
     #(i) * n_alpha(i) * n_b1(i) ** 2 * (-0.1D1 * n_rpsi(i - 1) + n_rpsi
     #(i + 1)) / ctfmp(i) / ctfm(i)) / hx - 0.5000000000000000D0 * (-0.2
     #500000000000000D0 / n_a1(i) ** 2 * (-0.1D1 * n_rpsi(i - 1) + n_rps
     #i(i + 1)) * (-0.1D1 * n_a1(i - 1) + n_a1(i + 1)) * n_alpha(i) * n_
     #b1(i) ** 2 - 0.2500000000000000D0 / n_a1(i) * (-0.1D1 * n_rpsi(i -
     # 1) + n_rpsi(i + 1)) * (n_alpha(i - 1) - 0.1D1 * n_alpha(i + 1)) *
     # n_b1(i) ** 2 + 0.5000000000000000D0 / n_a1(i) * (-0.1D1 * n_rpsi(
     #i - 1) + n_rpsi(i + 1)) * (-0.1D1 * n_b1(i - 1) + n_b1(i + 1)) * n
     #_alpha(i) * n_b1(i) + 0.1D1 / n_a1(i) * (n_rpsi(i - 1) - 0.2D1 * n
     #_rpsi(i) + n_rpsi(i + 1)) * n_alpha(i) * n_b1(i) ** 2) / ctfmp(i) 
     #** 2 / hx ** 2
      res = res + qb**2
      end do
      if (phys_bdy(2) .eq. 1) then
      do i=Nx, Nx, 1
      qb = (-0.1D1 * n_rPIb2(i) + np1_rPIb2(i)) / ht - 0.1D1 * myzero * 
     #x(i)
      res = res + qb**2
      end do
      endif
      res = sqrt(res/(1*Nx))
      END
