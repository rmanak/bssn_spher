      subroutine ire_hs(ctfm,ctfmp,n_A,n_Athth,n_Axx,n_B,n_K,n_Lamx,n_ps
     &i,n_rho,octfmp,Nx,hx,res)
      implicit none
      integer i
      integer Nx
      real*8 hx
      real*8 ctfm(Nx)
      real*8 ctfmp(Nx)
      real*8 n_A(Nx)
      real*8 n_Athth(Nx)
      real*8 n_Axx(Nx)
      real*8 n_B(Nx)
      real*8 n_K(Nx)
      real*8 n_Lamx(Nx)
      real*8 n_psi(Nx)
      real*8 n_rho(Nx)
      real*8 octfmp(Nx)
      real*8 res(Nx)
      real*8 qb
      do i=5, Nx, 1
      qb = -0.8333333333333333D-1 * n_psi(i) ** 5 * n_K(i) ** 2 + 0.2500
     #000000000000D0 * n_rho(i) * n_psi(i) ** 5 + 0.1250000000000000D0 /
     # n_A(i) ** 2 * n_Axx(i) ** 2 * n_psi(i) ** 5 + 0.2500000000000000D
     #0 / n_B(i) ** 2 * n_Athth(i) ** 2 * n_psi(i) ** 5 + (-0.3125000000
     #000000D-1 / n_A(i) * n_Lamx(i) * (n_A(i - 2) - 0.4D1 * n_A(i - 1) 
     #+ 0.3D1 * n_A(i)) / hx * n_psi(i) - 0.3125000000000000D-1 * (n_A(i
     # - 2) - 0.4D1 * n_A(i - 1) + 0.3D1 * n_A(i)) / hx / n_A(i) ** 2 * 
     #n_psi(i) * octfmp(i) - 0.6250000000000000D-1 / n_B(i) * (n_B(i - 2
     #) - 0.4D1 * n_B(i - 1) + 0.3D1 * n_B(i)) / hx * n_Lamx(i) * n_psi(
     #i) - 0.6250000000000000D-1 * (n_B(i - 2) - 0.4D1 * n_B(i - 1) + 0.
     #3D1 * n_B(i)) / hx / n_B(i) / n_A(i) * n_psi(i) * octfmp(i) + 0.62
     #50000000000000D-1 * (-0.1D1 * n_Lamx(i - 2) + 0.4D1 * n_Lamx(i - 1
     #) - 0.3D1 * n_Lamx(i)) / hx * n_psi(i) - 0.5000000000000000D0 / n_
     #A(i) * (n_psi(i - 2) - 0.4D1 * n_psi(i - 1) + 0.3D1 * n_psi(i)) / 
     #hx * octfmp(i)) / ctfmp(i) + (-0.2343750000000000D-1 / n_A(i) ** 3
     # * (n_A(i - 2) - 0.4D1 * n_A(i - 1) + 0.3D1 * n_A(i)) ** 2 / hx **
     # 2 * n_psi(i) - 0.1250000000000000D0 / n_A(i) ** 2 * (n_A(i - 2) -
     # 0.4D1 * n_A(i - 1) + 0.3D1 * n_A(i)) / hx ** 2 * (n_psi(i - 2) - 
     #0.4D1 * n_psi(i - 1) + 0.3D1 * n_psi(i)) - 0.1562500000000000D-1 /
     # n_B(i) ** 2 / n_A(i) * (n_B(i - 2) - 0.4D1 * n_B(i - 1) + 0.3D1 *
     # n_B(i)) ** 2 / hx ** 2 * n_psi(i) + 0.2500000000000000D0 / n_B(i)
     # / n_A(i) * (n_B(i - 2) - 0.4D1 * n_B(i - 1) + 0.3D1 * n_B(i)) / h
     #x ** 2 * (n_psi(i - 2) - 0.4D1 * n_psi(i - 1) + 0.3D1 * n_psi(i)) 
     #+ 0.8333333333333333D-1 / n_A(i) * (0.11D2 * n_psi(i - 4) - 0.56D2
     # * n_psi(i - 3) + 0.114D3 * n_psi(i - 2) - 0.104D3 * n_psi(i - 1) 
     #+ 0.35D2 * n_psi(i)) / hx ** 2 + 0.5208333333333333D-2 / n_A(i) **
     # 2 * (0.11D2 * n_A(i - 4) - 0.56D2 * n_A(i - 3) + 0.114D3 * n_A(i 
     #- 2) - 0.104D3 * n_A(i - 1) + 0.35D2 * n_A(i)) / hx ** 2 * n_psi(i
     #) + 0.1041666666666667D-1 / n_B(i) / n_A(i) * n_psi(i) * (0.11D2 *
     # n_B(i - 4) - 0.56D2 * n_B(i - 3) + 0.114D3 * n_B(i - 2) - 0.104D3
     # * n_B(i - 1) + 0.35D2 * n_B(i)) / hx ** 2) / ctfmp(i) ** 2 + (-0.
     #2500000000000000D0 * n_Lamx(i) * n_psi(i) + (0.6250000000000000D-1
     # * (n_A(i - 2) - 0.4D1 * n_A(i - 1) + 0.3D1 * n_A(i)) / hx / n_B(i
     #) / n_A(i) * n_psi(i) + 0.1250000000000000D0 * (n_B(i - 2) - 0.4D1
     # * n_B(i - 1) + 0.3D1 * n_B(i)) / hx / n_B(i) / n_A(i) * n_psi(i) 
     #+ 0.1D1 / n_A(i) * (n_psi(i - 2) - 0.4D1 * n_psi(i - 1) + 0.3D1 * 
     #n_psi(i)) / hx) / ctfmp(i)) / ctfm(i)
      res(i)=qb
      end do
      END
