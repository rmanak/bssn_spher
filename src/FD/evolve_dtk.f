      subroutine evolve_dtk(ctfmp,n_A,n_Athth,n_Axx,n_B,n_DDLthetatheta,
     &n_DDLxx,n_K,n_TS,n_alpha,n_beta,n_metthetatheta,n_metxx,n_rho,np1_
     &A,np1_Athth,np1_Axx,np1_B,np1_DDLthetatheta,np1_DDLxx,np1_K,np1_TS
     &,np1_alpha,np1_beta,np1_metthetatheta,np1_metxx,np1_rho,x,Nx,ht,hx
     &,myzero,zepsdis,phys_bdy,res)
      implicit none
      integer i
      integer Nx
      real*8 ht
      real*8 hx
      real*8 myzero
      real*8 zepsdis
      real*8 ctfmp(Nx)
      real*8 n_A(Nx)
      real*8 n_Athth(Nx)
      real*8 n_Axx(Nx)
      real*8 n_B(Nx)
      real*8 n_DDLthetatheta(Nx)
      real*8 n_DDLxx(Nx)
      real*8 n_K(Nx)
      real*8 n_TS(Nx)
      real*8 n_alpha(Nx)
      real*8 n_beta(Nx)
      real*8 n_metthetatheta(Nx)
      real*8 n_metxx(Nx)
      real*8 n_rho(Nx)
      real*8 np1_A(Nx)
      real*8 np1_Athth(Nx)
      real*8 np1_Axx(Nx)
      real*8 np1_B(Nx)
      real*8 np1_DDLthetatheta(Nx)
      real*8 np1_DDLxx(Nx)
      real*8 np1_K(Nx)
      real*8 np1_TS(Nx)
      real*8 np1_alpha(Nx)
      real*8 np1_beta(Nx)
      real*8 np1_metthetatheta(Nx)
      real*8 np1_metxx(Nx)
      real*8 np1_rho(Nx)
      real*8 x(Nx)
      integer phys_bdy(2)
      real*8 res(Nx)
      real*8 qb
      if (phys_bdy(1) .eq. 1) then
      do i=1, 1, 1
      qb = 0.1333333333333333D1 * np1_K(i + 1) - 0.3333333333333333D0 * 
     #np1_K(i + 2)
      res(i)=qb
      end do
      endif
      if (phys_bdy(1) .eq. 1) then
      do i=2, 2, 1
      qb = np1_K(i) + 0.96D2 / (0.32D2 * np1_alpha(i) * np1_K(i) * ht - 
     #0.21D2 * zepsdis - 0.96D2) * ht * (-0.1D1 * (n_K(i) - 0.1D1 * np1_
     #K(i)) / ht - 0.1666666666666667D0 * np1_alpha(i) * np1_K(i) ** 2 +
     # 0.2500000000000000D0 / ctfmp(i) * (np1_K(i - 1) - 0.1D1 * np1_K(i
     # + 1)) / hx * np1_beta(i) - 0.2500000000000000D0 * np1_alpha(i) * 
     #np1_rho(i) - 0.2500000000000000D0 * np1_alpha(i) * np1_TS(i) - 0.5
     #000000000000000D0 / np1_A(i) ** 2 * np1_Axx(i) ** 2 * np1_alpha(i)
     # - 0.1D1 / np1_B(i) ** 2 * np1_Athth(i) ** 2 * np1_alpha(i) + 0.1D
     #1 / np1_metthetatheta(i) * np1_DDLthetatheta(i) + 0.50000000000000
     #00D0 / np1_metxx(i) * np1_DDLxx(i) - 0.1666666666666667D0 * n_alph
     #a(i) * n_K(i) ** 2 + 0.2500000000000000D0 / ctfmp(i) * (n_K(i - 1)
     # - 0.1D1 * n_K(i + 1)) / hx * n_beta(i) - 0.2500000000000000D0 * n
     #_alpha(i) * n_rho(i) - 0.2500000000000000D0 * n_alpha(i) * n_TS(i)
     # - 0.5000000000000000D0 / n_A(i) ** 2 * n_Axx(i) ** 2 * n_alpha(i)
     # - 0.1D1 / n_B(i) ** 2 * n_Athth(i) ** 2 * n_alpha(i) + 0.1D1 / n_
     #metthetatheta(i) * n_DDLthetatheta(i) + 0.5000000000000000D0 / n_m
     #etxx(i) * n_DDLxx(i) + 0.3125000000000000D-1 * zepsdis / ht * (0.7
     #D1 * np1_K(i) + np1_K(i + 2) - 0.4D1 * np1_K(i + 1) - 0.4D1 * np1_
     #K(i - 1)) + 0.3125000000000000D-1 * zepsdis / ht * (0.7D1 * n_K(i)
     # + n_K(i + 2) - 0.4D1 * n_K(i + 1) - 0.4D1 * n_K(i - 1)))
      res(i)=qb
      end do
      endif
      do i=3, Nx-2, 1
      qb = np1_K(i) + 0.48D2 / (0.16D2 * np1_alpha(i) * np1_K(i) * ht - 
     #0.9D1 * zepsdis - 0.48D2) * ht * (-0.1D1 * (n_K(i) - 0.1D1 * np1_K
     #(i)) / ht - 0.1666666666666667D0 * np1_alpha(i) * np1_K(i) ** 2 + 
     #0.2500000000000000D0 / ctfmp(i) * (np1_K(i - 1) - 0.1D1 * np1_K(i 
     #+ 1)) / hx * np1_beta(i) - 0.2500000000000000D0 * np1_alpha(i) * n
     #p1_rho(i) - 0.2500000000000000D0 * np1_alpha(i) * np1_TS(i) - 0.50
     #00000000000000D0 / np1_A(i) ** 2 * np1_Axx(i) ** 2 * np1_alpha(i) 
     #- 0.1D1 / np1_B(i) ** 2 * np1_Athth(i) ** 2 * np1_alpha(i) + 0.1D1
     # / np1_metthetatheta(i) * np1_DDLthetatheta(i) + 0.500000000000000
     #0D0 / np1_metxx(i) * np1_DDLxx(i) - 0.1666666666666667D0 * n_alpha
     #(i) * n_K(i) ** 2 + 0.2500000000000000D0 / ctfmp(i) * (n_K(i - 1) 
     #- 0.1D1 * n_K(i + 1)) / hx * n_beta(i) - 0.2500000000000000D0 * n_
     #alpha(i) * n_rho(i) - 0.2500000000000000D0 * n_alpha(i) * n_TS(i) 
     #- 0.5000000000000000D0 / n_A(i) ** 2 * n_Axx(i) ** 2 * n_alpha(i) 
     #- 0.1D1 / n_B(i) ** 2 * n_Athth(i) ** 2 * n_alpha(i) + 0.1D1 / n_m
     #etthetatheta(i) * n_DDLthetatheta(i) + 0.5000000000000000D0 / n_me
     #txx(i) * n_DDLxx(i) + 0.3125000000000000D-1 * zepsdis / ht * (0.6D
     #1 * np1_K(i) + np1_K(i + 2) + np1_K(i - 2) - 0.4D1 * np1_K(i + 1) 
     #- 0.4D1 * np1_K(i - 1)) + 0.3125000000000000D-1 * zepsdis / ht * (
     #0.6D1 * n_K(i) + n_K(i + 2) + n_K(i - 2) - 0.4D1 * n_K(i + 1) - 0.
     #4D1 * n_K(i - 1)))
      res(i)=qb
      end do
      do i=Nx-1, Nx-1, 1
      qb = np1_K(i) + 0.3D1 / (np1_alpha(i) * np1_K(i) * ht - 0.3D1) * h
     #t * (-0.1D1 * (n_K(i) - 0.1D1 * np1_K(i)) / ht - 0.166666666666666
     #7D0 * np1_alpha(i) * np1_K(i) ** 2 + 0.2500000000000000D0 / ctfmp(
     #i) * (np1_K(i - 1) - 0.1D1 * np1_K(i + 1)) / hx * np1_beta(i) - 0.
     #2500000000000000D0 * np1_alpha(i) * np1_rho(i) - 0.250000000000000
     #0D0 * np1_alpha(i) * np1_TS(i) - 0.5000000000000000D0 / np1_A(i) *
     #* 2 * np1_Axx(i) ** 2 * np1_alpha(i) - 0.1D1 / np1_B(i) ** 2 * np1
     #_Athth(i) ** 2 * np1_alpha(i) + 0.1D1 / np1_metthetatheta(i) * np1
     #_DDLthetatheta(i) + 0.5000000000000000D0 / np1_metxx(i) * np1_DDLx
     #x(i) - 0.1666666666666667D0 * n_alpha(i) * n_K(i) ** 2 + 0.2500000
     #000000000D0 / ctfmp(i) * (n_K(i - 1) - 0.1D1 * n_K(i + 1)) / hx * 
     #n_beta(i) - 0.2500000000000000D0 * n_alpha(i) * n_rho(i) - 0.25000
     #00000000000D0 * n_alpha(i) * n_TS(i) - 0.5000000000000000D0 / n_A(
     #i) ** 2 * n_Axx(i) ** 2 * n_alpha(i) - 0.1D1 / n_B(i) ** 2 * n_Ath
     #th(i) ** 2 * n_alpha(i) + 0.1D1 / n_metthetatheta(i) * n_DDLthetat
     #heta(i) + 0.5000000000000000D0 / n_metxx(i) * n_DDLxx(i))
      res(i)=qb
      end do
      if (phys_bdy(2) .eq. 1) then
      do i=Nx, Nx, 1
      qb = np1_K(i) - 0.1D1 * ht * (-0.1D1 * (n_K(i) - 0.1D1 * np1_K(i))
     # / ht - 0.1D1 * myzero * x(i))
      res(i)=qb
      end do
      endif
      END
