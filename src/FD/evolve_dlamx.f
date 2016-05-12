      subroutine evolve_dlamx(ctfm,ctfmp,n_A,n_Athth,n_Axx,n_B,n_K,n_Lam
     &x,n_Sx,n_alpha,n_beta,n_divbeta,n_phi,np1_A,np1_Athth,np1_Axx,np1_
     &B,np1_K,np1_Lamx,np1_Sx,np1_alpha,np1_beta,np1_divbeta,np1_phi,oct
     &fmp,x,Nx,ht,hx,myzero,vee,zepsdis,phys_bdy,res)
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
      real*8 n_A(Nx)
      real*8 n_Athth(Nx)
      real*8 n_Axx(Nx)
      real*8 n_B(Nx)
      real*8 n_K(Nx)
      real*8 n_Lamx(Nx)
      real*8 n_Sx(Nx)
      real*8 n_alpha(Nx)
      real*8 n_beta(Nx)
      real*8 n_divbeta(Nx)
      real*8 n_phi(Nx)
      real*8 np1_A(Nx)
      real*8 np1_Athth(Nx)
      real*8 np1_Axx(Nx)
      real*8 np1_B(Nx)
      real*8 np1_K(Nx)
      real*8 np1_Lamx(Nx)
      real*8 np1_Sx(Nx)
      real*8 np1_alpha(Nx)
      real*8 np1_beta(Nx)
      real*8 np1_divbeta(Nx)
      real*8 np1_phi(Nx)
      real*8 octfmp(Nx)
      real*8 x(Nx)
      integer phys_bdy(2)
      real*8 res(Nx)
      real*8 qb
      if (phys_bdy(1) .eq. 1) then
      do i=1, 1, 1
      qb = -0.1D1 * myzero * x(i)
      res(i)=qb
      end do
      endif
      if (phys_bdy(1) .eq. 1) then
      do i=2, 2, 1
      qb = np1_Lamx(i) + 0.96D2 / (-0.15D2 * zepsdis * hx ** 2 + 0.832D3
     # * ht) * hx ** 2 * ht * (-0.8333333333333333D-1 * (-0.35D2 * np1_L
     #amx(i - 1) + 0.104D3 * np1_Lamx(i) - 0.114D3 * np1_Lamx(i + 1) + 0
     #.56D2 * np1_Lamx(i + 2) - 0.11D2 * np1_Lamx(i + 3)) / hx ** 2 + 0.
     #3125000000000000D-1 * zepsdis / ht * (0.5D1 * np1_Lamx(i) + np1_La
     #mx(i + 2) - 0.4D1 * np1_Lamx(i + 1) - 0.4D1 * np1_Lamx(i - 1)) + 0
     #.3125000000000000D-1 * zepsdis / ht * (0.5D1 * n_Lamx(i) + n_Lamx(
     #i + 2) - 0.4D1 * n_Lamx(i + 1) - 0.4D1 * n_Lamx(i - 1)))
      res(i)=qb
      end do
      endif
      do i=3, Nx-2, 1
      qb = np1_Lamx(i) + 0.48D2 / (0.16D2 * vee * np1_divbeta(i) * ht * 
     #ctfmp(i) * hx - 0.9D1 * zepsdis * ctfmp(i) * hx + 0.12D2 * ht * np
     #1_beta(i - 1) - 0.12D2 * ht * np1_beta(i + 1) - 0.48D2 * hx * ctfm
     #p(i)) * ht * ctfmp(i) * hx * (-0.1D1 * (n_Lamx(i) - 0.1D1 * np1_La
     #mx(i)) / ht - 0.3333333333333333D0 * np1_divbeta(i) * np1_Lamx(i) 
     #* vee + 0.1D1 / np1_A(i) * np1_alpha(i) * np1_Sx(i) - 0.5000000000
     #000000D0 * (-0.4D1 / np1_B(i) / np1_A(i) * np1_Athth(i) * np1_alph
     #a(i) + 0.4D1 / np1_B(i) ** 2 * np1_alpha(i) * np1_Athth(i)) / ctfm
     #(i) + 0.1D1 / ctfm(i) ** 2 / np1_B(i) * np1_beta(i) - 0.5000000000
     #000000D0 * ((0.5000000000000000D0 / np1_A(i) ** 3 * np1_Axx(i) * n
     #p1_alpha(i) * (-0.1D1 * np1_A(i - 1) + np1_A(i + 1)) - 0.1D1 / np1
     #_B(i) ** 2 / np1_A(i) * (-0.1D1 * np1_B(i - 1) + np1_B(i + 1)) * n
     #p1_Athth(i) * np1_alpha(i) + 0.6666666666666667D0 / np1_A(i) * np1
     #_alpha(i) * (np1_K(i - 1) - 0.1D1 * np1_K(i + 1)) + 0.500000000000
     #0000D0 * np1_beta(i) * (-0.1D1 * np1_Lamx(i - 1) + np1_Lamx(i + 1)
     #) + 0.1D1 / np1_A(i) ** 2 * np1_Axx(i) * (np1_alpha(i - 1) - 0.1D1
     # * np1_alpha(i + 1)) - 0.5000000000000000D0 * np1_Lamx(i) * (-0.1D
     #1 * np1_beta(i - 1) + np1_beta(i + 1)) - 0.5000000000000000D0 * (-
     #0.1D1 * np1_beta(i - 1) + np1_beta(i + 1)) / np1_A(i) * octfmp(i) 
     #- 0.1666666666666667D0 / np1_A(i) * vee * (np1_divbeta(i - 1) - 0.
     #1D1 * np1_divbeta(i + 1)) - 0.6D1 / np1_A(i) ** 2 * np1_alpha(i) *
     # np1_Axx(i) * (np1_phi(i - 1) - 0.1D1 * np1_phi(i + 1))) / ctfmp(i
     #) + 0.1D1 / np1_B(i) * (-0.1D1 * np1_beta(i - 1) + np1_beta(i + 1)
     #) / ctfmp(i) / ctfm(i)) / hx + 0.5000000000000000D0 / np1_A(i) / c
     #tfmp(i) ** 2 * (-0.1D1 * np1_beta(i - 1) + 0.2D1 * np1_beta(i) - 0
     #.1D1 * np1_beta(i + 1)) / hx ** 2 - 0.3333333333333333D0 * n_divbe
     #ta(i) * n_Lamx(i) * vee + 0.1D1 / n_A(i) * n_alpha(i) * n_Sx(i) - 
     #0.5000000000000000D0 * (-0.4D1 / n_B(i) / n_A(i) * n_Athth(i) * n_
     #alpha(i) + 0.4D1 / n_B(i) ** 2 * n_alpha(i) * n_Athth(i)) / ctfm(i
     #) + 0.1D1 / ctfm(i) ** 2 / n_B(i) * n_beta(i) - 0.5000000000000000
     #D0 * ((0.5000000000000000D0 / n_A(i) ** 3 * n_Axx(i) * n_alpha(i) 
     #* (-0.1D1 * n_A(i - 1) + n_A(i + 1)) - 0.1D1 / n_B(i) ** 2 / n_A(i
     #) * (-0.1D1 * n_B(i - 1) + n_B(i + 1)) * n_Athth(i) * n_alpha(i) +
     # 0.6666666666666667D0 / n_A(i) * n_alpha(i) * (n_K(i - 1) - 0.1D1 
     #* n_K(i + 1)) + 0.5000000000000000D0 * n_beta(i) * (-0.1D1 * n_Lam
     #x(i - 1) + n_Lamx(i + 1)) + 0.1D1 / n_A(i) ** 2 * n_Axx(i) * (n_al
     #pha(i - 1) - 0.1D1 * n_alpha(i + 1)) - 0.5000000000000000D0 * n_La
     #mx(i) * (-0.1D1 * n_beta(i - 1) + n_beta(i + 1)) - 0.5000000000000
     #000D0 * (-0.1D1 * n_beta(i - 1) + n_beta(i + 1)) / n_A(i) * octfmp
     #(i) - 0.1666666666666667D0 / n_A(i) * vee * (n_divbeta(i - 1) - 0.
     #1D1 * n_divbeta(i + 1)) - 0.6D1 / n_A(i) ** 2 * n_alpha(i) * n_Axx
     #(i) * (n_phi(i - 1) - 0.1D1 * n_phi(i + 1))) / ctfmp(i) + 0.1D1 / 
     #n_B(i) * (-0.1D1 * n_beta(i - 1) + n_beta(i + 1)) / ctfmp(i) / ctf
     #m(i)) / hx + 0.5000000000000000D0 / n_A(i) / ctfmp(i) ** 2 * (-0.1
     #D1 * n_beta(i - 1) + 0.2D1 * n_beta(i) - 0.1D1 * n_beta(i + 1)) / 
     #hx ** 2 + 0.3125000000000000D-1 * zepsdis / ht * (0.6D1 * np1_Lamx
     #(i) + np1_Lamx(i + 2) + np1_Lamx(i - 2) - 0.4D1 * np1_Lamx(i + 1) 
     #- 0.4D1 * np1_Lamx(i - 1)) + 0.3125000000000000D-1 * zepsdis / ht 
     #* (0.6D1 * n_Lamx(i) + n_Lamx(i + 2) + n_Lamx(i - 2) - 0.4D1 * n_L
     #amx(i + 1) - 0.4D1 * n_Lamx(i - 1)))
      res(i)=qb
      end do
      do i=Nx-1, Nx-1, 1
      qb = np1_Lamx(i) + 0.12D2 / (0.4D1 * vee * np1_divbeta(i) * ht * c
     #tfmp(i) * hx + 0.3D1 * ht * np1_beta(i - 1) - 0.3D1 * ht * np1_bet
     #a(i + 1) - 0.12D2 * hx * ctfmp(i)) * ht * ctfmp(i) * hx * (-0.1D1 
     #* (n_Lamx(i) - 0.1D1 * np1_Lamx(i)) / ht - 0.3333333333333333D0 * 
     #np1_divbeta(i) * np1_Lamx(i) * vee + 0.1D1 / np1_A(i) * np1_alpha(
     #i) * np1_Sx(i) - 0.5000000000000000D0 * (-0.4D1 / np1_B(i) / np1_A
     #(i) * np1_Athth(i) * np1_alpha(i) + 0.4D1 / np1_B(i) ** 2 * np1_al
     #pha(i) * np1_Athth(i)) / ctfm(i) + 0.1D1 / ctfm(i) ** 2 / np1_B(i)
     # * np1_beta(i) - 0.5000000000000000D0 * ((0.5000000000000000D0 / n
     #p1_A(i) ** 3 * np1_Axx(i) * np1_alpha(i) * (-0.1D1 * np1_A(i - 1) 
     #+ np1_A(i + 1)) - 0.1D1 / np1_B(i) ** 2 / np1_A(i) * (-0.1D1 * np1
     #_B(i - 1) + np1_B(i + 1)) * np1_Athth(i) * np1_alpha(i) + 0.666666
     #6666666667D0 / np1_A(i) * np1_alpha(i) * (np1_K(i - 1) - 0.1D1 * n
     #p1_K(i + 1)) + 0.5000000000000000D0 * np1_beta(i) * (-0.1D1 * np1_
     #Lamx(i - 1) + np1_Lamx(i + 1)) + 0.1D1 / np1_A(i) ** 2 * np1_Axx(i
     #) * (np1_alpha(i - 1) - 0.1D1 * np1_alpha(i + 1)) - 0.500000000000
     #0000D0 * np1_Lamx(i) * (-0.1D1 * np1_beta(i - 1) + np1_beta(i + 1)
     #) - 0.5000000000000000D0 * (-0.1D1 * np1_beta(i - 1) + np1_beta(i 
     #+ 1)) / np1_A(i) * octfmp(i) - 0.1666666666666667D0 / np1_A(i) * v
     #ee * (np1_divbeta(i - 1) - 0.1D1 * np1_divbeta(i + 1)) - 0.6D1 / n
     #p1_A(i) ** 2 * np1_alpha(i) * np1_Axx(i) * (np1_phi(i - 1) - 0.1D1
     # * np1_phi(i + 1))) / ctfmp(i) + 0.1D1 / np1_B(i) * (-0.1D1 * np1_
     #beta(i - 1) + np1_beta(i + 1)) / ctfmp(i) / ctfm(i)) / hx + 0.5000
     #000000000000D0 / np1_A(i) / ctfmp(i) ** 2 * (-0.1D1 * np1_beta(i -
     # 1) + 0.2D1 * np1_beta(i) - 0.1D1 * np1_beta(i + 1)) / hx ** 2 - 0
     #.3333333333333333D0 * n_divbeta(i) * n_Lamx(i) * vee + 0.1D1 / n_A
     #(i) * n_alpha(i) * n_Sx(i) - 0.5000000000000000D0 * (-0.4D1 / n_B(
     #i) / n_A(i) * n_Athth(i) * n_alpha(i) + 0.4D1 / n_B(i) ** 2 * n_al
     #pha(i) * n_Athth(i)) / ctfm(i) + 0.1D1 / ctfm(i) ** 2 / n_B(i) * n
     #_beta(i) - 0.5000000000000000D0 * ((0.5000000000000000D0 / n_A(i) 
     #** 3 * n_Axx(i) * n_alpha(i) * (-0.1D1 * n_A(i - 1) + n_A(i + 1)) 
     #- 0.1D1 / n_B(i) ** 2 / n_A(i) * (-0.1D1 * n_B(i - 1) + n_B(i + 1)
     #) * n_Athth(i) * n_alpha(i) + 0.6666666666666667D0 / n_A(i) * n_al
     #pha(i) * (n_K(i - 1) - 0.1D1 * n_K(i + 1)) + 0.5000000000000000D0 
     #* n_beta(i) * (-0.1D1 * n_Lamx(i - 1) + n_Lamx(i + 1)) + 0.1D1 / n
     #_A(i) ** 2 * n_Axx(i) * (n_alpha(i - 1) - 0.1D1 * n_alpha(i + 1)) 
     #- 0.5000000000000000D0 * n_Lamx(i) * (-0.1D1 * n_beta(i - 1) + n_b
     #eta(i + 1)) - 0.5000000000000000D0 * (-0.1D1 * n_beta(i - 1) + n_b
     #eta(i + 1)) / n_A(i) * octfmp(i) - 0.1666666666666667D0 / n_A(i) *
     # vee * (n_divbeta(i - 1) - 0.1D1 * n_divbeta(i + 1)) - 0.6D1 / n_A
     #(i) ** 2 * n_alpha(i) * n_Axx(i) * (n_phi(i - 1) - 0.1D1 * n_phi(i
     # + 1))) / ctfmp(i) + 0.1D1 / n_B(i) * (-0.1D1 * n_beta(i - 1) + n_
     #beta(i + 1)) / ctfmp(i) / ctfm(i)) / hx + 0.5000000000000000D0 / n
     #_A(i) / ctfmp(i) ** 2 * (-0.1D1 * n_beta(i - 1) + 0.2D1 * n_beta(i
     #) - 0.1D1 * n_beta(i + 1)) / hx ** 2)
      res(i)=qb
      end do
      if (phys_bdy(2) .eq. 1) then
      do i=Nx, Nx, 1
      qb = np1_Lamx(i) - 0.1D1 * ht * (-0.1D1 * (n_Lamx(i) - 0.1D1 * np1
     #_Lamx(i)) / ht - 0.1D1 * myzero * x(i))
      res(i)=qb
      end do
      endif
      END
