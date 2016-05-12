c-------------------------------
c interpolator for periodic boundry condition
c input:
c    x_c,f_c,n_c
c    (values on corse grid)
c output:
c    f, rc
c    (value of function on finer grid (must be 2x finer)
c    rc is the return code for un/succesfull interpolation
c  
c  method:
c     1: linear interpolation
c     2: C^2 smooth spiline interpolation
c---------------------------------------------
        subroutine intplt(x_c,f_c,n_c,x,f,n,method,rc)
        
        implicit none
        integer method
        integer n_c, n
        real*8 x_c(n_c), f_c(n_c)
        real*8 x(n), f(n)
        real*8 dx
        real*8 xmid
        integer rc 
        real*8 fmid
        integer i,k,l,j
        real*8 f_c_dd(n_c)
        if (method .ne. 1 .and. method .ne. 2) then
           write(0,*) 'interpolator: not supported method'
           rc=0
           return
        endif
c-------
c Checking for Compatibilities
c-------
            if ( (n-1)/(n_c-1) .ne. 2) then
                write(0,*) 'interpolater: not compatible grids'
                rc=0
                return
            endif 
            if ( x(n) .ne. x_c(n_c) ) then
               write(0,*) 'interpolator: incompatible grid length'
               write(0,*) ' x_max = ',x(n)
               write(0,*) ' xc_max = ', x_c(n_c)
               !rc=0
               !return
            endif
c--------------------------------------    
        if ( method .eq. 1 ) then
           do i=2, n_c
               k = 2*i - 1
               l= k - 1
               f(k) = f_c(i)
               f(l) = (f_c(i) + f_c(i-1))/2 
            enddo
            f(1) = f_c(1)
            rc=1
            return
        endif 

        if ( method .eq. 2 ) then
          dx = ( x_c(n_c)-x_c(1) )/(n_c-1)
          call Divxx(f_c, dx , n_c, f_c_dd)
          do i=3, n_c
             k=2*i-1
             l=k-1
             f(k) = f_c(i)
             xmid = x(l)
             j=i-1
             fmid = f_c_dd(j)*(xmid-x_c(j-1))**3/(6*dx) 
     &              + f_c_dd(j-1)*(x_c(j)-xmid)**3/(6*dx)
     &              + (f_c(j)/dx - f_c_dd(j)*dx/6)*(xmid-x_c(j-1))
     &              + (f_c(j-1)/dx - f_c_dd(j-1)*dx/6)*(x_c(j)-xmid)
             f(l) = fmid
          enddo

          i=2
             k=2*i-1
             l=k-1
             f(k) = f_c(i)
             xmid = x(l)
             j=i-1
             fmid = f_c_dd(j)*(xmid-(x_c(j)-dx))**3/(6*dx) 
     &              + f_c_dd(n_c-1)*(x_c(j)-xmid)**3/(6*dx)
     &              + (f_c(j)/dx - f_c_dd(j)*dx/6)*(xmid-(x_c(j)-dx))
     &              + (f_c(n_c-1)/dx - f_c_dd(n_c-1)*dx/6)*(x_c(j)-xmid)
             f(2) = fmid

          f(1)=f_c(1)  
          f(n)=f_c(n_c)
          rc=1
          return
        endif

        end
