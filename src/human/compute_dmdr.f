      subroutine compute_dmdr(ctfm, rho, Nx, dmdr)
         implicit none
         integer Nx
         real*8 rho(Nx), ctfm(Nx), dmdr(Nx)
         integer i
         do i=1, Nx, 1
           dmdr(i) = rho(i)*ctfm(i)**2
         end do
         END
