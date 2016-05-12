      subroutine rel_diff(f1,f2,N,res)
         implicit none
         integer N
         real*8 f1(N)
         real*8 f2(N)
         real*8 res
         integer i
         res =0.0D0
         do i=1, N, 1
           if ( dabs(f1(i)-f2(i))/abs(f1(i)) > res )  then
              res = dabs(f1(i)-f2(i))/abs(f1(i))
           end if
         end do
         END
