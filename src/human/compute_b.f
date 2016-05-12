      subroutine compute_b(A,Nx,B)
         implicit none
         integer Nx
         real*8 A(Nx)
         real*8 B(Nx)
         integer i
         do i=1,Nx,1
           B(i) = 1.0D0/dsqrt(A(i))
         end do
      END
