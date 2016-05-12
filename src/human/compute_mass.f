      subroutine compute_mass(A,B,psi,ctfm,ctfmp,rho,b1,a1,mass,
     &  mass2,tmr
     &  ,maxtmr,bhmass,dx,Nx)
       implicit none
       real*8 pi
       parameter ( pi = 3.141592653589793D0 )
       integer i
       integer Nx
       real*8 ctfm(Nx)
       real*8 rho(Nx)
       real*8 ctfmp(Nx)
       real*8 b1(Nx)
       real*8 avg
       real*8 dx
       real*8 mass(Nx)
       real*8 mass2(Nx)
       real*8 dR
       real*8 a1(Nx)
       real*8 grr
       real*8 r,b1p,rp
       real*8 tmr(Nx)
       real*8 maxtmr
       real*8 A(Nx)
       real*8 B(Nx)
       real*8 psi(Nx)
       real*8 bhmass


       mass2(1) = 0.0D0
       do i=2, Nx-1, 1
        dR = ctfmp(i)*dx
        avg = 0.5*( ctfm(i)**2*rho(i) +
     &              ctfm(i-1)**2*rho(i-1)
     &            )
        mass2(i) = mass2(i-1) + 4*pi*avg*dR
       end do
       mass2(Nx) = mass2(Nx-1)

        mass(1) = 0.0D0
        tmr(1) = 0.0D0
        maxtmr = 0.0D0
        bhmass = 0.0D0
        do i=2, Nx-1, 1
          rp = b1(i)*ctfm(i)
          mass(i) = ctfm(i)*psi(i)**2*sqrt(B(i))/2.0D0 * (
     &     1.0D0 - B(i)/A(i)* ( 
     &     1.0D0 + ctfm(i)*(B(i)-B(i-1))/dx/ctfmp(i)/(2.0D0*B(i))
     &     + 2*ctfm(i) * (psi(i)-psi(i-1))/dx/ctfmp(i)/psi(i)
     &                        )**2
     &     )

          tmr(i) = 2.0D0*mass(i)/rp
          if ( maxtmr .lt. tmr(i) ) then
             maxtmr = tmr(i)
             bhmass = mass(i)
          end if
        end do
        mass(Nx) = mass(Nx-1)
        tmr(Nx) = tmr(Nx-1)
       END
