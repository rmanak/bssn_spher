      subroutine evolve_eqipib2(ctfm,ctfmp,n_a1,n_alpha,n_b1,n_beta,n_iP
     &Ib2,n_ipsi,np1_a1,np1_alpha,np1_b1,np1_beta,np1_iPIb2,np1_ipsi,oct
     &fmp,x,Nx,ht,hx,mass,myzero,zepsdis,phys_bdy,res)
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
      real*8 n_iPIb2(Nx)
      real*8 n_ipsi(Nx)
      real*8 np1_a1(Nx)
      real*8 np1_alpha(Nx)
      real*8 np1_b1(Nx)
      real*8 np1_beta(Nx)
      real*8 np1_iPIb2(Nx)
      real*8 np1_ipsi(Nx)
      real*8 octfmp(Nx)
      real*8 x(Nx)
      integer phys_bdy(2)
      real*8 res(Nx)
      real*8 qb
      if (phys_bdy(1) .eq. 1) then
      do i=1, 1, 1
      qb = 0.1333333333333333D1 * np1_iPIb2(i + 1) - 0.3333333333333333D
     #0 * np1_iPIb2(i + 2)
      res(i)=qb
      end do
      endif
      if (phys_bdy(1) .eq. 1) then
      do i=2, 2, 1
      qb = np1_iPIb2(i) + 0.32D2 / (0.32D2 * np1_beta(i) * ht * ctfmp(i)
     # * hx - 0.7D1 * zepsdis * ctfm(i) * ctfmp(i) * hx - 0.8D1 * ht * c
     #tfm(i) * np1_beta(i - 1) + 0.8D1 * ht * ctfm(i) * np1_beta(i + 1) 
     #- 0.32D2 * ctfm(i) * ctfmp(i) * hx) * ht * ctfm(i) * ctfmp(i) * hx
     # * (-0.1D1 * (n_iPIb2(i) - 0.1D1 * np1_iPIb2(i)) / ht + 0.50000000
     #00000000D0 * np1_a1(i) * mass * np1_ipsi(i) * np1_alpha(i) * np1_b
     #1(i) ** 2 - 0.1D1 * np1_iPIb2(i) * np1_beta(i) / ctfm(i) - 0.50000
     #00000000000D0 * ((0.5000000000000000D0 * np1_iPIb2(i) * (-0.1D1 * 
     #np1_beta(i - 1) + np1_beta(i + 1)) - 0.5000000000000000D0 * np1_be
     #ta(i) * (np1_iPIb2(i - 1) - 0.1D1 * np1_iPIb2(i + 1)) + 0.50000000
     #00000000D0 / np1_a1(i) * (np1_ipsi(i - 1) - 0.1D1 * np1_ipsi(i + 1
     #)) * octfmp(i) * np1_alpha(i) * np1_b1(i) ** 2) / ctfmp(i) - 0.1D1
     # / np1_a1(i) * np1_alpha(i) * np1_b1(i) ** 2 * (np1_ipsi(i - 1) - 
     #0.1D1 * np1_ipsi(i + 1)) / ctfmp(i) / ctfm(i)) / hx - 0.5000000000
     #000000D0 * (0.2500000000000000D0 / np1_a1(i) ** 2 * (np1_ipsi(i - 
     #1) - 0.1D1 * np1_ipsi(i + 1)) * (-0.1D1 * np1_a1(i - 1) + np1_a1(i
     # + 1)) * np1_alpha(i) * np1_b1(i) ** 2 - 0.2500000000000000D0 / np
     #1_a1(i) * (np1_ipsi(i - 1) - 0.1D1 * np1_ipsi(i + 1)) * (-0.1D1 * 
     #np1_alpha(i - 1) + np1_alpha(i + 1)) * np1_b1(i) ** 2 - 0.50000000
     #00000000D0 / np1_a1(i) * (np1_ipsi(i - 1) - 0.1D1 * np1_ipsi(i + 1
     #)) * (-0.1D1 * np1_b1(i - 1) + np1_b1(i + 1)) * np1_alpha(i) * np1
     #_b1(i) - 0.1D1 / np1_a1(i) * (-0.1D1 * np1_ipsi(i - 1) + 0.2D1 * n
     #p1_ipsi(i) - 0.1D1 * np1_ipsi(i + 1)) * np1_alpha(i) * np1_b1(i) *
     #* 2) / ctfmp(i) ** 2 / hx ** 2 + 0.5000000000000000D0 * n_a1(i) * 
     #mass * n_ipsi(i) * n_alpha(i) * n_b1(i) ** 2 - 0.1D1 * n_iPIb2(i) 
     #* n_beta(i) / ctfm(i) - 0.5000000000000000D0 * ((0.500000000000000
     #0D0 * n_iPIb2(i) * (-0.1D1 * n_beta(i - 1) + n_beta(i + 1)) - 0.50
     #00000000000000D0 * n_beta(i) * (n_iPIb2(i - 1) - 0.1D1 * n_iPIb2(i
     # + 1)) + 0.5000000000000000D0 / n_a1(i) * (n_ipsi(i - 1) - 0.1D1 *
     # n_ipsi(i + 1)) * octfmp(i) * n_alpha(i) * n_b1(i) ** 2) / ctfmp(i
     #) - 0.1D1 / n_a1(i) * n_alpha(i) * n_b1(i) ** 2 * (n_ipsi(i - 1) -
     # 0.1D1 * n_ipsi(i + 1)) / ctfmp(i) / ctfm(i)) / hx - 0.50000000000
     #00000D0 * (0.2500000000000000D0 / n_a1(i) ** 2 * (n_ipsi(i - 1) - 
     #0.1D1 * n_ipsi(i + 1)) * (-0.1D1 * n_a1(i - 1) + n_a1(i + 1)) * n_
     #alpha(i) * n_b1(i) ** 2 - 0.2500000000000000D0 / n_a1(i) * (n_ipsi
     #(i - 1) - 0.1D1 * n_ipsi(i + 1)) * (-0.1D1 * n_alpha(i - 1) + n_al
     #pha(i + 1)) * n_b1(i) ** 2 - 0.5000000000000000D0 / n_a1(i) * (n_i
     #psi(i - 1) - 0.1D1 * n_ipsi(i + 1)) * (-0.1D1 * n_b1(i - 1) + n_b1
     #(i + 1)) * n_alpha(i) * n_b1(i) - 0.1D1 / n_a1(i) * (-0.1D1 * n_ip
     #si(i - 1) + 0.2D1 * n_ipsi(i) - 0.1D1 * n_ipsi(i + 1)) * n_alpha(i
     #) * n_b1(i) ** 2) / ctfmp(i) ** 2 / hx ** 2 + 0.3125000000000000D-
     #1 * zepsdis / ht * (0.7D1 * np1_iPIb2(i) + np1_iPIb2(i + 2) - 0.4D
     #1 * np1_iPIb2(i + 1) - 0.4D1 * np1_iPIb2(i - 1)) + 0.3125000000000
     #000D-1 * zepsdis / ht * (0.7D1 * n_iPIb2(i) + n_iPIb2(i + 2) - 0.4
     #D1 * n_iPIb2(i + 1) - 0.4D1 * n_iPIb2(i - 1)))
      res(i)=qb
      end do
      endif
      do i=3, Nx-2, 1
      qb = np1_iPIb2(i) + 0.16D2 / (0.16D2 * np1_beta(i) * ht * ctfmp(i)
     # * hx - 0.3D1 * zepsdis * ctfm(i) * ctfmp(i) * hx - 0.4D1 * ht * c
     #tfm(i) * np1_beta(i - 1) + 0.4D1 * ht * ctfm(i) * np1_beta(i + 1) 
     #- 0.16D2 * ctfm(i) * ctfmp(i) * hx) * ht * ctfm(i) * ctfmp(i) * hx
     # * (-0.1D1 * (n_iPIb2(i) - 0.1D1 * np1_iPIb2(i)) / ht + 0.50000000
     #00000000D0 * np1_a1(i) * mass * np1_ipsi(i) * np1_alpha(i) * np1_b
     #1(i) ** 2 - 0.1D1 * np1_iPIb2(i) * np1_beta(i) / ctfm(i) - 0.50000
     #00000000000D0 * ((0.5000000000000000D0 * np1_iPIb2(i) * (-0.1D1 * 
     #np1_beta(i - 1) + np1_beta(i + 1)) - 0.5000000000000000D0 * np1_be
     #ta(i) * (np1_iPIb2(i - 1) - 0.1D1 * np1_iPIb2(i + 1)) + 0.50000000
     #00000000D0 / np1_a1(i) * (np1_ipsi(i - 1) - 0.1D1 * np1_ipsi(i + 1
     #)) * octfmp(i) * np1_alpha(i) * np1_b1(i) ** 2) / ctfmp(i) - 0.1D1
     # / np1_a1(i) * np1_alpha(i) * np1_b1(i) ** 2 * (np1_ipsi(i - 1) - 
     #0.1D1 * np1_ipsi(i + 1)) / ctfmp(i) / ctfm(i)) / hx - 0.5000000000
     #000000D0 * (0.2500000000000000D0 / np1_a1(i) ** 2 * (np1_ipsi(i - 
     #1) - 0.1D1 * np1_ipsi(i + 1)) * (-0.1D1 * np1_a1(i - 1) + np1_a1(i
     # + 1)) * np1_alpha(i) * np1_b1(i) ** 2 - 0.2500000000000000D0 / np
     #1_a1(i) * (np1_ipsi(i - 1) - 0.1D1 * np1_ipsi(i + 1)) * (-0.1D1 * 
     #np1_alpha(i - 1) + np1_alpha(i + 1)) * np1_b1(i) ** 2 - 0.50000000
     #00000000D0 / np1_a1(i) * (np1_ipsi(i - 1) - 0.1D1 * np1_ipsi(i + 1
     #)) * (-0.1D1 * np1_b1(i - 1) + np1_b1(i + 1)) * np1_alpha(i) * np1
     #_b1(i) - 0.1D1 / np1_a1(i) * (-0.1D1 * np1_ipsi(i - 1) + 0.2D1 * n
     #p1_ipsi(i) - 0.1D1 * np1_ipsi(i + 1)) * np1_alpha(i) * np1_b1(i) *
     #* 2) / ctfmp(i) ** 2 / hx ** 2 + 0.5000000000000000D0 * n_a1(i) * 
     #mass * n_ipsi(i) * n_alpha(i) * n_b1(i) ** 2 - 0.1D1 * n_iPIb2(i) 
     #* n_beta(i) / ctfm(i) - 0.5000000000000000D0 * ((0.500000000000000
     #0D0 * n_iPIb2(i) * (-0.1D1 * n_beta(i - 1) + n_beta(i + 1)) - 0.50
     #00000000000000D0 * n_beta(i) * (n_iPIb2(i - 1) - 0.1D1 * n_iPIb2(i
     # + 1)) + 0.5000000000000000D0 / n_a1(i) * (n_ipsi(i - 1) - 0.1D1 *
     # n_ipsi(i + 1)) * octfmp(i) * n_alpha(i) * n_b1(i) ** 2) / ctfmp(i
     #) - 0.1D1 / n_a1(i) * n_alpha(i) * n_b1(i) ** 2 * (n_ipsi(i - 1) -
     # 0.1D1 * n_ipsi(i + 1)) / ctfmp(i) / ctfm(i)) / hx - 0.50000000000
     #00000D0 * (0.2500000000000000D0 / n_a1(i) ** 2 * (n_ipsi(i - 1) - 
     #0.1D1 * n_ipsi(i + 1)) * (-0.1D1 * n_a1(i - 1) + n_a1(i + 1)) * n_
     #alpha(i) * n_b1(i) ** 2 - 0.2500000000000000D0 / n_a1(i) * (n_ipsi
     #(i - 1) - 0.1D1 * n_ipsi(i + 1)) * (-0.1D1 * n_alpha(i - 1) + n_al
     #pha(i + 1)) * n_b1(i) ** 2 - 0.5000000000000000D0 / n_a1(i) * (n_i
     #psi(i - 1) - 0.1D1 * n_ipsi(i + 1)) * (-0.1D1 * n_b1(i - 1) + n_b1
     #(i + 1)) * n_alpha(i) * n_b1(i) - 0.1D1 / n_a1(i) * (-0.1D1 * n_ip
     #si(i - 1) + 0.2D1 * n_ipsi(i) - 0.1D1 * n_ipsi(i + 1)) * n_alpha(i
     #) * n_b1(i) ** 2) / ctfmp(i) ** 2 / hx ** 2 + 0.3125000000000000D-
     #1 * zepsdis / ht * (0.6D1 * np1_iPIb2(i) + np1_iPIb2(i + 2) + np1_
     #iPIb2(i - 2) - 0.4D1 * np1_iPIb2(i + 1) - 0.4D1 * np1_iPIb2(i - 1)
     #) + 0.3125000000000000D-1 * zepsdis / ht * (0.6D1 * n_iPIb2(i) + n
     #_iPIb2(i + 2) + n_iPIb2(i - 2) - 0.4D1 * n_iPIb2(i + 1) - 0.4D1 * 
     #n_iPIb2(i - 1)))
      res(i)=qb
      end do
      do i=Nx-1, Nx-1, 1
      qb = np1_iPIb2(i) + 0.4D1 / (0.4D1 * np1_beta(i) * ht * ctfmp(i) *
     # hx - 0.1D1 * ht * ctfm(i) * np1_beta(i - 1) + ht * ctfm(i) * np1_
     #beta(i + 1) - 0.4D1 * ctfm(i) * ctfmp(i) * hx) * ht * ctfm(i) * ct
     #fmp(i) * hx * (-0.1D1 * (n_iPIb2(i) - 0.1D1 * np1_iPIb2(i)) / ht +
     # 0.5000000000000000D0 * np1_a1(i) * mass * np1_ipsi(i) * np1_alpha
     #(i) * np1_b1(i) ** 2 - 0.1D1 * np1_iPIb2(i) * np1_beta(i) / ctfm(i
     #) - 0.5000000000000000D0 * ((0.5000000000000000D0 * np1_iPIb2(i) *
     # (-0.1D1 * np1_beta(i - 1) + np1_beta(i + 1)) - 0.5000000000000000
     #D0 * np1_beta(i) * (np1_iPIb2(i - 1) - 0.1D1 * np1_iPIb2(i + 1)) +
     # 0.5000000000000000D0 / np1_a1(i) * (np1_ipsi(i - 1) - 0.1D1 * np1
     #_ipsi(i + 1)) * octfmp(i) * np1_alpha(i) * np1_b1(i) ** 2) / ctfmp
     #(i) - 0.1D1 / np1_a1(i) * np1_alpha(i) * np1_b1(i) ** 2 * (np1_ips
     #i(i - 1) - 0.1D1 * np1_ipsi(i + 1)) / ctfmp(i) / ctfm(i)) / hx - 0
     #.5000000000000000D0 * (0.2500000000000000D0 / np1_a1(i) ** 2 * (np
     #1_ipsi(i - 1) - 0.1D1 * np1_ipsi(i + 1)) * (-0.1D1 * np1_a1(i - 1)
     # + np1_a1(i + 1)) * np1_alpha(i) * np1_b1(i) ** 2 - 0.250000000000
     #0000D0 / np1_a1(i) * (np1_ipsi(i - 1) - 0.1D1 * np1_ipsi(i + 1)) *
     # (-0.1D1 * np1_alpha(i - 1) + np1_alpha(i + 1)) * np1_b1(i) ** 2 -
     # 0.5000000000000000D0 / np1_a1(i) * (np1_ipsi(i - 1) - 0.1D1 * np1
     #_ipsi(i + 1)) * (-0.1D1 * np1_b1(i - 1) + np1_b1(i + 1)) * np1_alp
     #ha(i) * np1_b1(i) - 0.1D1 / np1_a1(i) * (-0.1D1 * np1_ipsi(i - 1) 
     #+ 0.2D1 * np1_ipsi(i) - 0.1D1 * np1_ipsi(i + 1)) * np1_alpha(i) * 
     #np1_b1(i) ** 2) / ctfmp(i) ** 2 / hx ** 2 + 0.5000000000000000D0 *
     # n_a1(i) * mass * n_ipsi(i) * n_alpha(i) * n_b1(i) ** 2 - 0.1D1 * 
     #n_iPIb2(i) * n_beta(i) / ctfm(i) - 0.5000000000000000D0 * ((0.5000
     #000000000000D0 * n_iPIb2(i) * (-0.1D1 * n_beta(i - 1) + n_beta(i +
     # 1)) - 0.5000000000000000D0 * n_beta(i) * (n_iPIb2(i - 1) - 0.1D1 
     #* n_iPIb2(i + 1)) + 0.5000000000000000D0 / n_a1(i) * (n_ipsi(i - 1
     #) - 0.1D1 * n_ipsi(i + 1)) * octfmp(i) * n_alpha(i) * n_b1(i) ** 2
     #) / ctfmp(i) - 0.1D1 / n_a1(i) * n_alpha(i) * n_b1(i) ** 2 * (n_ip
     #si(i - 1) - 0.1D1 * n_ipsi(i + 1)) / ctfmp(i) / ctfm(i)) / hx - 0.
     #5000000000000000D0 * (0.2500000000000000D0 / n_a1(i) ** 2 * (n_ips
     #i(i - 1) - 0.1D1 * n_ipsi(i + 1)) * (-0.1D1 * n_a1(i - 1) + n_a1(i
     # + 1)) * n_alpha(i) * n_b1(i) ** 2 - 0.2500000000000000D0 / n_a1(i
     #) * (n_ipsi(i - 1) - 0.1D1 * n_ipsi(i + 1)) * (-0.1D1 * n_alpha(i 
     #- 1) + n_alpha(i + 1)) * n_b1(i) ** 2 - 0.5000000000000000D0 / n_a
     #1(i) * (n_ipsi(i - 1) - 0.1D1 * n_ipsi(i + 1)) * (-0.1D1 * n_b1(i 
     #- 1) + n_b1(i + 1)) * n_alpha(i) * n_b1(i) - 0.1D1 / n_a1(i) * (-0
     #.1D1 * n_ipsi(i - 1) + 0.2D1 * n_ipsi(i) - 0.1D1 * n_ipsi(i + 1)) 
     #* n_alpha(i) * n_b1(i) ** 2) / ctfmp(i) ** 2 / hx ** 2)
      res(i)=qb
      end do
      if (phys_bdy(2) .eq. 1) then
      do i=Nx, Nx, 1
      qb = np1_iPIb2(i) - 0.1D1 * ht * (-0.1D1 * (n_iPIb2(i) - 0.1D1 * n
     #p1_iPIb2(i)) / ht - 0.1D1 * myzero * x(i))
      res(i)=qb
      end do
      endif
      END
