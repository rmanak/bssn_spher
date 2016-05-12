c------------------------------------------------------------------
c   Calculates all first/second x/t derivative of metric functions
c------------------------------------------------------------------
          subroutine calcall ( 
     &              Gxx_a  , Gxx  , Gxx_r,
     &              Gyy_a  , Gyy  , Gyy_r,
     &              Gzz_a  , Gzz  , Gzz_r,
     &              alpha_a, alpha, alpha_r,
     &              beta_a , beta , beta_r,
     &              Gxx_x, Gxx_t, Gxx_xt, Gxx_tt, Gxx_xx,
     &              Gyy_x, Gyy_t, Gyy_xt, Gyy_tt, Gyy_xx,
     &              Gzz_x, Gzz_t, Gzz_xt, Gzz_tt, Gzz_xx,
     &              alpha_x, alpha_t, alpha_xt, alpha_tt, alpha_xx,
     &              beta_x , beta_t , beta_xt , beta_tt , beta_xx,
     &              n, dx, dt
     &                     )
          implicit none
          integer n
          real*8 dx, dt
          real*8 Gxx_a(1:n), Gxx(1:n), Gxx_r(1:n)
          real*8 Gyy_a(1:n), Gyy(1:n), Gyy_r(1:n)
          real*8 Gzz_a(1:n), Gzz(1:n), Gzz_r(1:n)
 
          real*8 alpha_a(1:n), alpha(1:n), alpha_r(1:n)
          
          real*8 beta_a(1:n), beta(1:n), beta_r(1:n)
          
          real*8 Gxx_x (1:n), Gyy_x (1:n), Gzz_x (1:n)
          real*8 Gxx_t (1:n), Gyy_t (1:n), Gzz_t (1:n)
          real*8 Gxx_xt(1:n), Gyy_xt(1:n), Gzz_xt(1:n)
          real*8 Gxx_tt(1:n), Gyy_tt(1:n), Gzz_tt(1:n) 
          real*8 Gxx_xx(1:n), Gyy_xx(1:n), Gzz_xx(1:n)
 
          real*8 alpha_x (1:n), beta_x (1:n)
          real*8 alpha_t (1:n), beta_t (1:n)
          real*8 alpha_xx(1:n), beta_xx(1:n)
          real*8 alpha_xt(1:n), beta_xt(1:n)
          real*8 alpha_tt(1:n), beta_tt(1:n)   
          
          call calcallDiv(Gxx_a,   Gxx ,    Gxx_r , 
     &                        dx ,   dt,    n   ,
     &                        Gxx_x,   Gxx_t,   Gxx_xx, Gxx_tt , Gxx_xt)
          
          call calcallDiv(Gyy_a,   Gyy ,    Gyy_r , 
     &                        dx ,   dt,    n   ,
     &                        Gyy_x,   Gyy_t,   Gyy_xx, Gyy_tt , Gyy_xt)

          call calcallDiv(Gzz_a,   Gzz ,    Gzz_r , 
     &                        dx ,   dt,    n   ,
     &                        Gzz_x,   Gzz_t,   Gzz_xx, Gzz_tt , Gzz_xt)

          call calcallDiv(alpha_a,   alpha ,    alpha_r , 
     &                        dx ,   dt,    n   ,
     &                        alpha_x,   alpha_t,   alpha_xx, alpha_tt 
     &                          , alpha_xt)
          call calcallDiv(beta_a,   beta ,    beta_r , 
     &                        dx ,   dt,    n   ,
     &                        beta_x,   beta_t,   beta_xx, beta_tt 
     &                          , beta_xt)
            
          return
          end

c-----------------------------------------------------------------
c calculates all kind of derivative required in our problem for 
c a given grid function
c all O(h^2) approximation
c----------------------------------------------------------------
        subroutine calcallDiv(f_a,   f ,    f_r , 
     &                        dx ,   dt,    n   ,
     &                        f_x,   f_t,   f_xx, f_tt , f_xt)

        implicit none
        
        integer n
        real*8 f_a(1:n), f(1:n), f_r(1:n)
        real*8 dx, dt
        real*8 f_x(1:n), f_t(1:n), f_xx(1:n), f_tt(1:n), f_xt(1:n)
        
        call Divx(f, dx, n, f_x)
        call Divt(f_a, f_r, dt, n, f_t)  
        call Divtt(f_a, f, f_r, dt, n, f_tt)
        call Divxx(f, dx, n, f_xx)
        call Divxt(f_a, f_r, dx, dt, n, f_xt)

        return

        end
c-----------------------------------------------------------------
c calculates x derivative of f using left and right values O(h^2)
c assuming periodic boundry condition
c-----------------------------------------------------------------

         subroutine Divx(f, dx, n, f_x)
         
         implicit none
         
         integer n, i
         real*8 f(1:n), f_x(1:n)
         real*8 dx 

         do i=2, n-1
            f_x(i) = ( f(i+1) - f(i-1) ) / (2*dx)
         enddo 

         i = 1
            f_x(i) = ( f(i+1) - f(n-1) ) / (2*dx)

         i = n
            f_x(i) = ( f(1+1) - f(i-1) ) / (2*dx)

         return
         
         end 
c-------------------------------------------------------------------
c calculates time derivative of f using advanced and retarded values
c O(h^2) approximation
c-------------------------------------------------------------------

         subroutine Divt(f_a, f_r, dt, n, f_t)

         implicit none
         integer n, i
         real*8 dt
         real*8 f_a(1:n), f_r(1:n), f_t(1:n)

         do i=1, n
             f_t(i) = ( f_a(i) - f_r(i) ) / (2*dt)
         enddo
         
         return

         end
c---------------------------------------------------------------------
c  calculates second time derivative of f using advanced and retarded
c  and current values. O(h^2) approximation
c---------------------------------------------------------------------


        subroutine Divtt(f_a, f, f_r, dt, n, f_tt)

        implicit none

        real*8 dt
        integer n, i
        real*8 f_a(1:n), f(1:n), f_r(1:n), f_tt(1:n)
        
        do i=1, n
           f_tt(i) = ( f_a(i) - 2*f(i) + f_r(i) ) / (dt*dt)
        enddo
        
        return

        end
c----------------------------------------------------------------
c calculates second x derivative of f using left and right values
c O(h^2) approximation 
c----------------------------------------------------------------

        subroutine Divxx(f, dx, n, f_xx)

        implicit none
        integer n, i
        real*8 dx
        real*8 f(1:n), f_xx(1:n)
        
        do i=2, n-1
            f_xx(i) = ( f(i+1) - 2*f(i) + f(i-1) ) / (dx*dx)
        enddo

        i = 1
            f_xx(i) = ( f(i+1) - 2*f(i) + f(n-1) ) / (dx*dx)

        i = n
            f_xx(i) = ( f(1+1) - 2*f(i) + f(i-1) ) / (dx*dx)

        return

        end
c--------------------------------------------------------------
c calculates mixed x-t derivative of f using a/r values
c O(h^2) approximation
c--------------------------------------------------------------

        subroutine Divxt(f_a, f_r, dx, dt, n, f_xt)

        implicit none
        integer n, i
        real*8 dx, dt
        real*8 f_a(1:n), f_r(1:n), f_xt(1:n)

        do i=2, n-1
           f_xt(i) = ( f_a(i+1) - f_a(i-1) - f_r(i+1) + f_r(i-1) )
     &                      / (4*dt*dx)
        enddo

        i = 1
           f_xt(i) = ( f_a(i+1) - f_a(n-1) - f_r(i+1) + f_r(n-1) )
     &                      / (4*dt*dx)

        i = n 
           f_xt(i) = ( f_a(1+1) - f_a(i-1) - f_r(1+1) + f_r(i-1) )
     &                      / (4*dt*dx)

        return
        
        end
