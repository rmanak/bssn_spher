      subroutine evolve_daax(n_Axx,n_DDLthetatheta,n_DDLxx,n_JSthth,n_JS
     &xx,n_Pxx,n_RRthetatheta,n_RRxx,n_Uxx,n_alpha,n_em4phi,n_metthetath
     &eta,n_metxx,np1_Axx,np1_DDLthetatheta,np1_DDLxx,np1_JSthth,np1_JSx
     &x,np1_Pxx,np1_RRthetatheta,np1_RRxx,np1_Uxx,np1_alpha,np1_em4phi,n
     &p1_metthetatheta,np1_metxx,x,Nx,ht,hx,myzero,zepsdis,phys_bdy,res)
      implicit none
      integer i
      integer Nx
      real*8 ht
      real*8 hx
      real*8 myzero
      real*8 zepsdis
      real*8 n_Axx(Nx)
      real*8 n_DDLthetatheta(Nx)
      real*8 n_DDLxx(Nx)
      real*8 n_JSthth(Nx)
      real*8 n_JSxx(Nx)
      real*8 n_Pxx(Nx)
      real*8 n_RRthetatheta(Nx)
      real*8 n_RRxx(Nx)
      real*8 n_Uxx(Nx)
      real*8 n_alpha(Nx)
      real*8 n_em4phi(Nx)
      real*8 n_metthetatheta(Nx)
      real*8 n_metxx(Nx)
      real*8 np1_Axx(Nx)
      real*8 np1_DDLthetatheta(Nx)
      real*8 np1_DDLxx(Nx)
      real*8 np1_JSthth(Nx)
      real*8 np1_JSxx(Nx)
      real*8 np1_Pxx(Nx)
      real*8 np1_RRthetatheta(Nx)
      real*8 np1_RRxx(Nx)
      real*8 np1_Uxx(Nx)
      real*8 np1_alpha(Nx)
      real*8 np1_em4phi(Nx)
      real*8 np1_metthetatheta(Nx)
      real*8 np1_metxx(Nx)
      real*8 x(Nx)
      integer phys_bdy(2)
      real*8 res(Nx)
      real*8 qb
      if (phys_bdy(1) .eq. 1) then
      do i=1, 1, 1
      qb = -0.1D1 * myzero * x(i)
      res(i)=qb
      end do
      endif
      if (phys_bdy(1) .eq. 1) then
      do i=2, 2, 1
      qb = np1_Axx(i) - 0.32D2 / (0.7D1 * zepsdis * hx + 0.64D2 * ht) * 
     #hx * ht * (-0.5000000000000000D0 * (0.3D1 * np1_Axx(i - 1) - 0.4D1
     # * np1_Axx(i) + np1_Axx(i + 1)) / hx + hx * myzero + 0.31250000000
     #00000D-1 * zepsdis / ht * (0.7D1 * np1_Axx(i) + np1_Axx(i + 2) - 0
     #.4D1 * np1_Axx(i + 1) - 0.4D1 * np1_Axx(i - 1)) + 0.31250000000000
     #00D-1 * zepsdis / ht * (0.7D1 * n_Axx(i) + n_Axx(i + 2) - 0.4D1 * 
     #n_Axx(i + 1) - 0.4D1 * n_Axx(i - 1)))
      res(i)=qb
      end do
      endif
      do i=3, Nx-2, 1
      qb = np1_Axx(i) - 0.16D2 / (0.16D2 + 0.3D1 * zepsdis) * ht * (-0.1
     #D1 * (n_Axx(i) - 0.1D1 * np1_Axx(i)) / ht - 0.3333333333333333D0 *
     # np1_alpha(i) * np1_RRxx(i) * np1_em4phi(i) + 0.3333333333333333D0
     # * np1_em4phi(i) * np1_alpha(i) * np1_JSxx(i) + 0.3333333333333333
     #D0 / np1_metthetatheta(i) * np1_alpha(i) * np1_metxx(i) * np1_RRth
     #etatheta(i) * np1_em4phi(i) - 0.3333333333333333D0 / np1_metthetat
     #heta(i) * np1_alpha(i) * np1_metxx(i) * np1_JSthth(i) * np1_em4phi
     #(i) - 0.5000000000000000D0 * np1_alpha(i) * np1_Uxx(i) + 0.3333333
     #333333333D0 * np1_DDLxx(i) * np1_em4phi(i) - 0.3333333333333333D0 
     #/ np1_metthetatheta(i) * np1_DDLthetatheta(i) * np1_metxx(i) * np1
     #_em4phi(i) - 0.5000000000000000D0 * np1_Pxx(i) - 0.333333333333333
     #3D0 * n_alpha(i) * n_RRxx(i) * n_em4phi(i) + 0.3333333333333333D0 
     #* n_em4phi(i) * n_alpha(i) * n_JSxx(i) + 0.3333333333333333D0 / n_
     #metthetatheta(i) * n_alpha(i) * n_metxx(i) * n_RRthetatheta(i) * n
     #_em4phi(i) - 0.3333333333333333D0 / n_metthetatheta(i) * n_alpha(i
     #) * n_metxx(i) * n_JSthth(i) * n_em4phi(i) - 0.5000000000000000D0 
     #* n_alpha(i) * n_Uxx(i) + 0.3333333333333333D0 * n_DDLxx(i) * n_em
     #4phi(i) - 0.3333333333333333D0 / n_metthetatheta(i) * n_DDLthetath
     #eta(i) * n_metxx(i) * n_em4phi(i) - 0.5000000000000000D0 * n_Pxx(i
     #) + 0.3125000000000000D-1 * zepsdis / ht * (0.6D1 * np1_Axx(i) + n
     #p1_Axx(i + 2) + np1_Axx(i - 2) - 0.4D1 * np1_Axx(i + 1) - 0.4D1 * 
     #np1_Axx(i - 1)) + 0.3125000000000000D-1 * zepsdis / ht * (0.6D1 * 
     #n_Axx(i) + n_Axx(i + 2) + n_Axx(i - 2) - 0.4D1 * n_Axx(i + 1) - 0.
     #4D1 * n_Axx(i - 1)))
      res(i)=qb
      end do
      do i=Nx-1, Nx-1, 1
      qb = np1_Axx(i) - 0.1D1 * ht * (-0.1D1 * (n_Axx(i) - 0.1D1 * np1_A
     #xx(i)) / ht - 0.3333333333333333D0 * np1_alpha(i) * np1_RRxx(i) * 
     #np1_em4phi(i) + 0.3333333333333333D0 * np1_em4phi(i) * np1_alpha(i
     #) * np1_JSxx(i) + 0.3333333333333333D0 / np1_metthetatheta(i) * np
     #1_alpha(i) * np1_metxx(i) * np1_RRthetatheta(i) * np1_em4phi(i) - 
     #0.3333333333333333D0 / np1_metthetatheta(i) * np1_alpha(i) * np1_m
     #etxx(i) * np1_JSthth(i) * np1_em4phi(i) - 0.5000000000000000D0 * n
     #p1_alpha(i) * np1_Uxx(i) + 0.3333333333333333D0 * np1_DDLxx(i) * n
     #p1_em4phi(i) - 0.3333333333333333D0 / np1_metthetatheta(i) * np1_D
     #DLthetatheta(i) * np1_metxx(i) * np1_em4phi(i) - 0.500000000000000
     #0D0 * np1_Pxx(i) - 0.3333333333333333D0 * n_alpha(i) * n_RRxx(i) *
     # n_em4phi(i) + 0.3333333333333333D0 * n_em4phi(i) * n_alpha(i) * n
     #_JSxx(i) + 0.3333333333333333D0 / n_metthetatheta(i) * n_alpha(i) 
     #* n_metxx(i) * n_RRthetatheta(i) * n_em4phi(i) - 0.333333333333333
     #3D0 / n_metthetatheta(i) * n_alpha(i) * n_metxx(i) * n_JSthth(i) *
     # n_em4phi(i) - 0.5000000000000000D0 * n_alpha(i) * n_Uxx(i) + 0.33
     #33333333333333D0 * n_DDLxx(i) * n_em4phi(i) - 0.3333333333333333D0
     # / n_metthetatheta(i) * n_DDLthetatheta(i) * n_metxx(i) * n_em4phi
     #(i) - 0.5000000000000000D0 * n_Pxx(i))
      res(i)=qb
      end do
      if (phys_bdy(2) .eq. 1) then
      do i=Nx, Nx, 1
      qb = np1_Axx(i) - 0.1D1 * ht * (-0.1D1 * (n_Axx(i) - 0.1D1 * np1_A
     #xx(i)) / ht - 0.1D1 * myzero * x(i))
      res(i)=qb
      end do
      endif
      END
