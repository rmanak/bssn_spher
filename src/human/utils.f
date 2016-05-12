
        subroutine bailout(Max2MovR,time,flag)

        implicit none
        integer ubh
        integer getu
        integer flag
        real*8 time, Max2MovR

        if ( flag .eq. 1 ) then
            ubh = 5 ! Dangerous
            open(ubh, file='DISPERSAL', form='formatted')
            write(ubh,*) 1
            close(ubh)
            write(0,*) 'DISPERSAL'
            write(0,*) 't=',time
            write(0,*) 'Max 2M/r =',Max2MovR
        endif

        if ( flag .eq. 2 ) then
              ubh = 5 !Dangerous
              open(ubh, file='BH', form='formatted')
              write(ubh,*) 1
              close(ubh)
              write(0,*) 'Black Hole detected at time:'
              write(0,*) 't=',time
        endif

        return
        end

        subroutine check_BH(a,Nx,maxtmovr,time,BHT,rc)
        implicit none
        integer Nx
        real*8 a(Nx)
        real*8 maxtmovr
        real*8 time, BHT
        real*8 tmp1
        integer getu
        integer i
        integer ubh
        integer rc

        maxtmovr = 0.0D0
        rc = 0

        do i=1, Nx

           tmp1 = (1.0D0-1.0D0/a(i)**2)

           if (abs(tmp1) .gt. maxtmovr) then
              maxtmovr = abs(tmp1)
           endif

        enddo

        if (maxtmovr .gt. BHT) then
              ubh = 4 !Dangerous
              open(ubh, file='BH', form='formatted')
              write(ubh,*) 1
              close(ubh)
              write(0,*) 'Black Hole detected at time:'
              write(0,*) 't=',time
              write(0,*) 'Max 2M/r=', maxtmovr
              rc = 1
              return
        endif

        return
        END

