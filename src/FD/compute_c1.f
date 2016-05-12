      subroutine compute_c1(ctfm,ctfmp,n_A,n_B,n_Lamx,octfmp,x,Nx,hx,myz
     &ero,phys_bdy,res)
      implicit none
      integer i
      integer Nx
      real*8 hx
      real*8 myzero
      real*8 ctfm(Nx)
      real*8 ctfmp(Nx)
      real*8 n_A(Nx)
      real*8 n_B(Nx)
      real*8 n_Lamx(Nx)
      real*8 octfmp(Nx)
      real*8 x(Nx)
      integer phys_bdy(2)
      real*8 res(Nx)
      real*8 qb
      if (phys_bdy(1) .eq. 1) then
      do i=1, 1, 1
      qb = myzero * x(i)
      res(i)=qb
      end do
      endif
      do i=2, Nx-1, 1
      qb = -0.2500000000000000D0 * n_Lamx(i) / ctfm(i) + ((-0.3125000000
     #000000D-1 / n_A(i) * n_Lamx(i) * (-0.1D1 * n_A(i - 1) + n_A(i + 1)
     #) - 0.3125000000000000D-1 * (-0.1D1 * n_A(i - 1) + n_A(i + 1)) / n
     #_A(i) ** 2 * octfmp(i) - 0.6250000000000000D-1 / n_B(i) * (-0.1D1 
     #* n_B(i - 1) + n_B(i + 1)) * n_Lamx(i) - 0.6250000000000000D-1 * (
     #-0.1D1 * n_B(i - 1) + n_B(i + 1)) / n_B(i) / n_A(i) * octfmp(i) + 
     #0.6250000000000000D-1 * n_Lamx(i - 1) - 0.6250000000000000D-1 * n_
     #Lamx(i + 1)) / ctfmp(i) + (0.6250000000000000D-1 / n_B(i) / n_A(i)
     # * (-0.1D1 * n_A(i - 1) + n_A(i + 1)) + 0.1250000000000000D0 / n_A
     #(i) / n_B(i) * (-0.1D1 * n_B(i - 1) + n_B(i + 1))) / ctfmp(i) / ct
     #fm(i)) / hx + (-0.2343750000000000D-1 / n_A(i) ** 3 * (-0.1D1 * n_
     #A(i - 1) + n_A(i + 1)) ** 2 - 0.1562500000000000D-1 / n_B(i) ** 2 
     #/ n_A(i) * (-0.1D1 * n_B(i - 1) + n_B(i + 1)) ** 2 + 0.62500000000
     #00000D-1 / n_A(i) ** 2 * (n_A(i - 1) - 0.2D1 * n_A(i) + n_A(i + 1)
     #) - 0.1250000000000000D0 / n_B(i) / n_A(i) * (-0.1D1 * n_B(i - 1) 
     #+ 0.2D1 * n_B(i) - 0.1D1 * n_B(i + 1))) / ctfmp(i) ** 2 / hx ** 2
      res(i)=qb
      end do
      if (phys_bdy(2) .eq. 1) then
      do i=Nx, Nx, 1
      qb = myzero * x(i)
      res(i)=qb
      end do
      endif
      END
