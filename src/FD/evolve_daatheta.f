      subroutine evolve_daatheta(n_Athth,n_DDLthetatheta,n_DDLxx,n_JStht
     &h,n_JSxx,n_Pthetatheta,n_RRthetatheta,n_RRxx,n_Uthetatheta,n_alpha
     &,n_em4phi,n_metthetatheta,n_metxx,np1_Athth,np1_DDLthetatheta,np1_
     &DDLxx,np1_JSthth,np1_JSxx,np1_Pthetatheta,np1_RRthetatheta,np1_RRx
     &x,np1_Uthetatheta,np1_alpha,np1_em4phi,np1_metthetatheta,np1_metxx
     &,x,Nx,ht,hx,myzero,zepsdis,phys_bdy,res)
      implicit none
      integer i
      integer Nx
      real*8 ht
      real*8 hx
      real*8 myzero
      real*8 zepsdis
      real*8 n_Athth(Nx)
      real*8 n_DDLthetatheta(Nx)
      real*8 n_DDLxx(Nx)
      real*8 n_JSthth(Nx)
      real*8 n_JSxx(Nx)
      real*8 n_Pthetatheta(Nx)
      real*8 n_RRthetatheta(Nx)
      real*8 n_RRxx(Nx)
      real*8 n_Uthetatheta(Nx)
      real*8 n_alpha(Nx)
      real*8 n_em4phi(Nx)
      real*8 n_metthetatheta(Nx)
      real*8 n_metxx(Nx)
      real*8 np1_Athth(Nx)
      real*8 np1_DDLthetatheta(Nx)
      real*8 np1_DDLxx(Nx)
      real*8 np1_JSthth(Nx)
      real*8 np1_JSxx(Nx)
      real*8 np1_Pthetatheta(Nx)
      real*8 np1_RRthetatheta(Nx)
      real*8 np1_RRxx(Nx)
      real*8 np1_Uthetatheta(Nx)
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
      qb = np1_Athth(i) - 0.32D2 / (0.7D1 * zepsdis * hx + 0.64D2 * ht) 
     #* hx * ht * (hx * myzero + 0.5000000000000000D0 * (-0.3D1 * np1_At
     #hth(i - 1) + 0.4D1 * np1_Athth(i) - 0.1D1 * np1_Athth(i + 1)) / hx
     # + 0.3125000000000000D-1 * zepsdis / ht * (0.7D1 * np1_Athth(i) + 
     #np1_Athth(i + 2) - 0.4D1 * np1_Athth(i + 1) - 0.4D1 * np1_Athth(i 
     #- 1)) + 0.3125000000000000D-1 * zepsdis / ht * (0.7D1 * n_Athth(i)
     # + n_Athth(i + 2) - 0.4D1 * n_Athth(i + 1) - 0.4D1 * n_Athth(i - 1
     #)))
      res(i)=qb
      end do
      endif
      do i=3, Nx-2, 1
      qb = np1_Athth(i) - 0.16D2 / (0.16D2 + 0.3D1 * zepsdis) * ht * (-0
     #.1D1 * (n_Athth(i) - 0.1D1 * np1_Athth(i)) / ht + 0.16666666666666
     #67D0 / np1_metxx(i) * np1_alpha(i) * np1_metthetatheta(i) * np1_RR
     #xx(i) * np1_em4phi(i) - 0.1666666666666667D0 / np1_metxx(i) * np1_
     #em4phi(i) * np1_alpha(i) * np1_JSxx(i) * np1_metthetatheta(i) - 0.
     #1666666666666667D0 * np1_alpha(i) * np1_RRthetatheta(i) * np1_em4p
     #hi(i) + 0.1666666666666667D0 * np1_alpha(i) * np1_JSthth(i) * np1_
     #em4phi(i) - 0.5000000000000000D0 * np1_alpha(i) * np1_Uthetatheta(
     #i) - 0.1666666666666667D0 / np1_metxx(i) * np1_DDLxx(i) * np1_mett
     #hetatheta(i) * np1_em4phi(i) + 0.1666666666666667D0 * np1_DDLtheta
     #theta(i) * np1_em4phi(i) - 0.5000000000000000D0 * np1_Pthetatheta(
     #i) + 0.1666666666666667D0 / n_metxx(i) * n_alpha(i) * n_metthetath
     #eta(i) * n_RRxx(i) * n_em4phi(i) - 0.1666666666666667D0 / n_metxx(
     #i) * n_em4phi(i) * n_alpha(i) * n_JSxx(i) * n_metthetatheta(i) - 0
     #.1666666666666667D0 * n_alpha(i) * n_RRthetatheta(i) * n_em4phi(i)
     # + 0.1666666666666667D0 * n_alpha(i) * n_JSthth(i) * n_em4phi(i) -
     # 0.5000000000000000D0 * n_alpha(i) * n_Uthetatheta(i) - 0.16666666
     #66666667D0 / n_metxx(i) * n_DDLxx(i) * n_metthetatheta(i) * n_em4p
     #hi(i) + 0.1666666666666667D0 * n_DDLthetatheta(i) * n_em4phi(i) - 
     #0.5000000000000000D0 * n_Pthetatheta(i) + 0.3125000000000000D-1 * 
     #zepsdis / ht * (0.6D1 * np1_Athth(i) + np1_Athth(i + 2) + np1_Atht
     #h(i - 2) - 0.4D1 * np1_Athth(i + 1) - 0.4D1 * np1_Athth(i - 1)) + 
     #0.3125000000000000D-1 * zepsdis / ht * (0.6D1 * n_Athth(i) + n_Ath
     #th(i + 2) + n_Athth(i - 2) - 0.4D1 * n_Athth(i + 1) - 0.4D1 * n_At
     #hth(i - 1)))
      res(i)=qb
      end do
      do i=Nx-1, Nx-1, 1
      qb = np1_Athth(i) - 0.1D1 * ht * (-0.1D1 * (n_Athth(i) - 0.1D1 * n
     #p1_Athth(i)) / ht + 0.1666666666666667D0 / np1_metxx(i) * np1_alph
     #a(i) * np1_metthetatheta(i) * np1_RRxx(i) * np1_em4phi(i) - 0.1666
     #666666666667D0 / np1_metxx(i) * np1_em4phi(i) * np1_alpha(i) * np1
     #_JSxx(i) * np1_metthetatheta(i) - 0.1666666666666667D0 * np1_alpha
     #(i) * np1_RRthetatheta(i) * np1_em4phi(i) + 0.1666666666666667D0 *
     # np1_alpha(i) * np1_JSthth(i) * np1_em4phi(i) - 0.5000000000000000
     #D0 * np1_alpha(i) * np1_Uthetatheta(i) - 0.1666666666666667D0 / np
     #1_metxx(i) * np1_DDLxx(i) * np1_metthetatheta(i) * np1_em4phi(i) +
     # 0.1666666666666667D0 * np1_DDLthetatheta(i) * np1_em4phi(i) - 0.5
     #000000000000000D0 * np1_Pthetatheta(i) + 0.1666666666666667D0 / n_
     #metxx(i) * n_alpha(i) * n_metthetatheta(i) * n_RRxx(i) * n_em4phi(
     #i) - 0.1666666666666667D0 / n_metxx(i) * n_em4phi(i) * n_alpha(i) 
     #* n_JSxx(i) * n_metthetatheta(i) - 0.1666666666666667D0 * n_alpha(
     #i) * n_RRthetatheta(i) * n_em4phi(i) + 0.1666666666666667D0 * n_al
     #pha(i) * n_JSthth(i) * n_em4phi(i) - 0.5000000000000000D0 * n_alph
     #a(i) * n_Uthetatheta(i) - 0.1666666666666667D0 / n_metxx(i) * n_DD
     #Lxx(i) * n_metthetatheta(i) * n_em4phi(i) + 0.1666666666666667D0 *
     # n_DDLthetatheta(i) * n_em4phi(i) - 0.5000000000000000D0 * n_Pthet
     #atheta(i))
      res(i)=qb
      end do
      if (phys_bdy(2) .eq. 1) then
      do i=Nx, Nx, 1
      qb = np1_Athth(i) - 0.1D1 * ht * (-0.1D1 * (n_Athth(i) - 0.1D1 * n
     #p1_Athth(i)) / ht - 0.1D1 * myzero * x(i))
      res(i)=qb
      end do
      endif
      END
