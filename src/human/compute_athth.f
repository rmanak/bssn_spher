      subroutine compute_athth(Axx,A,B,Nx,Athth)
         implicit none
         integer Nx
         real*8 Axx(Nx)
         real*8 A(Nx)
         real*8 B(Nx)
         real*8 Athth(Nx)
         integer i
         ! Axx/A + 2*Athth/B = 0
         do i=1,Nx,1
           Athth(i) = -B(i)*Axx(i)/(2.0D0*A(i))
         end do
         END
