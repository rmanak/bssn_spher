      subroutine solve_hs3(ctfm,ctfmp,n_A,n_B,n_C1s,n_C5psi,n_psi,octfmp
     &,Nx,hx,myzero,phys_bdy,res)
      implicit none
      integer i
      integer Nx
      real*8 hx
      real*8 myzero
      real*8 ctfm(Nx)
      real*8 ctfmp(Nx)
      real*8 n_A(Nx)
      real*8 n_B(Nx)
      real*8 n_C1s(Nx)
      real*8 n_C5psi(Nx)
      real*8 n_psi(Nx)
      real*8 octfmp(Nx)
      integer phys_bdy(2)
      real*8 res(Nx)
      real*8 qb
      if (phys_bdy(1) .eq. 1) then
      do i=1, 1, 1
      qb = 0.1333333333333333D1 * n_psi(i + 1) - 0.3333333333333333D0 * 
     #n_psi(i + 2)
      res(i)=qb
      end do
      endif
      do i=2, Nx-1, 1
      qb = n_psi(i) - 0.1D1 / (n_C1s(i) * n_A(i) * ctfmp(i) ** 2 * hx **
     # 2 - 0.2D1) * n_A(i) * ctfmp(i) ** 2 * hx ** 2 * (n_C1s(i) * n_psi
     #(i) + n_C5psi(i) + (-0.5000000000000000D0 / n_A(i) / ctfmp(i) * (-
     #0.1D1 * n_psi(i - 1) + n_psi(i + 1)) * octfmp(i) + 0.1D1 / ctfm(i)
     # / n_A(i) / ctfmp(i) * (-0.1D1 * n_psi(i - 1) + n_psi(i + 1))) / h
     #x + (0.1D1 / n_A(i) * (n_psi(i - 1) - 0.2D1 * n_psi(i) + n_psi(i +
     # 1)) - 0.1250000000000000D0 / n_A(i) ** 2 * (-0.1D1 * n_A(i - 1) +
     # n_A(i + 1)) * (-0.1D1 * n_psi(i - 1) + n_psi(i + 1)) + 0.25000000
     #00000000D0 / n_B(i) / n_A(i) * (-0.1D1 * n_B(i - 1) + n_B(i + 1)) 
     #* (-0.1D1 * n_psi(i - 1) + n_psi(i + 1))) / ctfmp(i) ** 2 / hx ** 
     #2)
      res(i)=qb
      end do
      if (phys_bdy(2) .eq. 1) then
      do i=Nx, Nx, 1
      qb = 0.1D1 - 0.1D1 * myzero * ctfm(i)
      res(i)=qb
      end do
      endif
      END
