!
!------------------------------------------------------------------
! 11-20-2020
!
! Artificial Neural Network module for aladyn_pi.f code
! Open_MP version
!
! Vesselin Yamakov
! National Institute of Aerospace
! 100 Exploration Way,
! Hampton, VA 23666 
! phone: (757)-864-2850
! fax:   (757)-864-8911
! e-mail: yamakov@nianet.org
!------------------------------------------------------------------
! Use a trained NN to get pot. param for a specific potential form
!------------------------------------------------------------------
!
      MODULE ANN

      use sys_OMP
      use constants
      use sim_box
      use pot_module
      save

      integer, parameter :: Max_net_layers = 8

      integer :: n_set_ann, net_atom_types, iflag_ann
      integer :: net_layers, net_in, net_out, mG_dim, max_tri_index
      integer :: nBOP_params ! = Nodes_of_layer(net_layers) !

      real(kind=4) :: Rc_ann,d_ann,d4_ann,Gauss_ann,range_min_ann
      real(kind=4) :: ActFunc_shift

      integer :: Nodes_of_layer(Max_net_layers)

      real(kind=4), dimension(:), allocatable :: r0_value,r0G_value,
     1 base_pot_param

      real(kind=4), dimension(:,:,:), allocatable :: Gi_atom, dG_i
      !dir$ attributes align:64 :: Gi_atom
      !dir$ attributes align:64 :: dG_i
      real(kind=4), dimension(:,:), allocatable :: Gi_list, Gi_new,
     1  U3_of
      !dir$ attributes align:64 :: Gi_list
      !dir$ attributes align:64 :: Gi_new
      !dir$ attributes align:64 :: U3_of

      real(kind=4), dimension(:), allocatable :: U0,U1,U2,U3
      !dir$ attributes align:64 :: U1
      !dir$ attributes align:64 :: U2
      !dir$ attributes align:64 :: U3

      real(kind=4), dimension(:,:), allocatable :: W1_ann,W3_ann
      real(kind=4), dimension(:,:,:), allocatable :: W2_ann
      !dir$ attributes align:64 :: W1_ann
      !dir$ attributes align:64 :: W2_ann
      !dir$ attributes align:64 :: W3_ann

      real(kind=4), dimension(:), allocatable :: B1_ann, B3_ann
      real(kind=4), dimension(:,:), allocatable :: B2_ann
      !dir$ attributes align:64 :: B1_ann
      !dir$ attributes align:64 :: B2_ann
      !dir$ attributes align:64 :: B3_ann

      real(kind=4), dimension(:,:,:,:), allocatable::dBOP_param_dxij
      !dir$ attributes align:64 :: dBOP_param_dxij

      real(kind=4), dimension(:), allocatable :: buf_ann

      real(kind=4), dimension(:), allocatable :: r0Rc,r0pRc
      real(kind=4), dimension(:), allocatable :: U0f1,U1f1,U2f1
      real(kind=4), dimension(:), allocatable :: U0f2,U1f2,U2f2
      real(kind=4), dimension(:), allocatable :: U0f3,U1f3,U2f3

      real(kind=4), dimension(:,:), allocatable :: 
     1 Gi_dev,xr_ij0,xr_ij1,xr_ij2,xr_ij3
      real(kind=4), dimension(:,:,:), allocatable :: 
     1 xr_ij_dev,fsij_dev, dfuN_dev, dGi_dx,dGi_dy,dGi_dz
      real(kind=4), dimension(:,:,:), allocatable :: 
     1 dfs_rij_3D1,dfs_rij_3D2,dfs_rij_3D3
      real(kind=4), dimension(:,:,:,:), allocatable :: 
     1 dfs_rij_3D,Gi_3D_dev, dANN_dxij

      CONTAINS
!
! ---------------------------------------------------------------------
! Calculates ANN and its Analytical derivatives for the force calculation.
! ---------------------------------------------------------------------
!
      subroutine dANN_OMP(maxn_nbrs, Rcn)

      use atoms

      implicit real(kind=4) (a-h,o-z)

      integer, intent(in) :: maxn_nbrs
      real(kind=4), intent(in) :: Rcn
  
      real(kind=4) :: dgij1,dgij2,dgij3,U_x1,U_x2,U_x3
      real(kind=4) :: Gi_sum_01,Gi_sum_11,Gi_sum_21,Gi_sum_31,
     1 Gi_sum_41
      real(kind=4) :: Gi_sum_02,Gi_sum_12,Gi_sum_22,Gi_sum_32,
     1 Gi_sum_42
      real(kind=4) :: Gi_sum_03,Gi_sum_13,Gi_sum_23,Gi_sum_33,
     1 Gi_sum_43
      real(kind=4), dimension(3) :: dcos_ijk

      if(ihalt.ne.0) return

      Rc = Rcn; Rc2 = Rcn**2; d4 = d4_ann
      Sigma2 = 2.0*Gauss_ann**2
      r0Rc(:) = Rcn*r0_value(:)
      r0pRc(:) = Rcn + r0_value(:)

      h11 = h(1,1); h12 = h(1,2); h13 = h(1,3)
      h22 = h(2,2); h23 = h(2,3); h33 = h(3,3)
!
! --- Start I loop over ni atoms ---
!
!$OMP  PARALLEL DO DEFAULT(NONE)
!$OMP& PRIVATE(nj,ij_delta,sxn,syn,szn,sx0,sy0,sz0,xij,yij,zij,rij,r0)
!$OMP& PRIVATE(r1ij,r2ij,RcRij,RcRij3,RcRij4,RcRij5,fc_rij,dfsij)
!$OMP& PRIVATE(denominator,RijRo,expRijRo)
!$OMP& SHARED(natoms,maxn_nbrs, nbr_list, sx,sy,sz, fsij_dev)
!$OMP& SHARED(h11,h22,h33,h12,h13,h23, Rc,d4, Gauss_ann)
!$OMP& SHARED(n_set_ann, r0_value, Rc2,r0Rc,r0pRc, Sigma2)
!$OMP& SHARED(xr_ij0, xr_ij1, xr_ij2, xr_ij3)
!$OMP& SHARED(dfs_rij_3D1, dfs_rij_3D2, dfs_rij_3D3)
!$OMP& SCHEDULE(STATIC)

      do ni=1,natoms   ! loop over atoms !
       sxn = sx(ni); syn = sy(ni); szn = sz(ni)

       !DIR$ ASSUME (mod(maxn_nbrs,8) .eq. 0)
       do nb1 = 1,maxn_nbrs     ! Loop I: i - j bond !

        nj = nbr_list(nb1,ni)
        ij_delta = 1-min(abs(nj-ni),1) ! 1 if ni=nj; 0 if ni=/=nj !

        sx0 = sx(nj) - sxn
        sx0 = sx0 - nint(sx0)     ! make periodic along X !
        sy0 = sy(nj) - syn
        sy0 = sy0 - nint(sy0)     ! make periodic along Y !
        sz0 = sz(nj) - szn
        sz0 = sz0 - nint(sz0)     ! make periodic along Z !

        xij = h11*sx0 + h12*sy0 + h13*sz0
        yij = h22*sy0 + h23*sz0
        zij = h33*sz0
        r2ij = xij**2 + yij**2 + zij**2 + ij_delta*Rc2 ! rij+Rc when i=j
        rij = sqrt(r2ij)
        r1ij = 1.0/rij

        xr_ij0(nb1,ni) = rij
        xr_ij1(nb1,ni) = xij*r1ij
        xr_ij2(nb1,ni) = yij*r1ij
        xr_ij3(nb1,ni) = zij*r1ij

        RcRij = max(Rc-rij, 0.0)
        RcRij3 = RcRij**3
        RcRij4 = RcRij3*RcRij
        RcRij5 = RcRij4*RcRij
        fc_rij = RcRij4 / (d4+RcRij4)
        denominator = (Gauss_ann * (d4+RcRij4))**2

        do n_set = 1,n_set_ann
         r0 = r0_value(n_set)  ! do not use r0G_value() here !
         RijRo = rij-r0
         expRijRo = exp(-(RijRo/Gauss_ann)**2)
         fsij_dev(nb1,n_set,ni) = fc_rij * expRijRo / r0
         dfsij = 
     = (r0Rc(n_set) + rij*(rij-r0pRc(n_set)) - Sigma2)*d4 - RijRo*RcRij5
         dfsij = 2.0*dfsij*RcRij3*expRijRo / (r0*denominator)
         dfs_rij_3D1(n_set,nb1,ni) = dfsij*xr_ij1(nb1,ni) ! dfs_ij*xij/rij !
         dfs_rij_3D2(n_set,nb1,ni) = dfsij*xr_ij2(nb1,ni) ! dfs_ij*xij/rij !
         dfs_rij_3D3(n_set,nb1,ni) = dfsij*xr_ij3(nb1,ni) ! dfs_ij*xij/rij !
        enddo ! do n_set = 1,n_set_ann !

       enddo ! do nb1 = 1,maxn_nbrs !
      enddo ! do ni=1,natoms !

!$OMP END PARALLEL DO

!
! --- Start II loop over ni atoms ---
!
!$OMP  PARALLEL DO DEFAULT(NONE)
!$OMP& PRIVATE(Gi_sum0,Gi_sum1,Gi_sum2,Gi_sum3,Gi_sum4, fsij_fsik)
!$OMP& PRIVATE(icol, cos_ijk,cos2_ijk, Pl2,Pl4,Pl6)
!$OMP& SHARED(natoms, maxn_nbrs,n_set_ann, fsij_dev, Gi_dev)
!$OMP& SHARED(xr_ij0, xr_ij1, xr_ij2, xr_ij3)
!$OMP& SCHEDULE(STATIC)

      do ni=1,natoms   ! loop over atoms !
       do n_set = 1,n_set_ann

        Gi_sum0 = 0.0; Gi_sum1 = 0.0; Gi_sum2 = 0.0
        Gi_sum3 = 0.0; Gi_sum4 = 0.0

        !DIR$ VECTOR UNALIGNED
        !DIR$ ASSUME (mod(maxn_nbrs,8) .eq. 0)
        do nb1 = 1,maxn_nbrs     ! Loop I: i - j bond !
         !DIR$ VECTOR UNALIGNED
         !DIR$ ASSUME (mod(maxn_nbrs,8) .eq. 0)
         !$OMP simd simdlen(8)
         do nb2 = 1,maxn_nbrs     ! Loop I: i - j bond !

          fsij_fsik = fsij_dev(nb1,n_set,ni) * fsij_dev(nb2,n_set,ni)

          cos_ijk = xr_ij1(nb1,ni)*xr_ij1(nb2,ni) + 
     +    xr_ij2(nb1,ni)*xr_ij2(nb2,ni) + xr_ij3(nb1,ni)*xr_ij3(nb2,ni)
          cos2_ijk = cos_ijk**2
          Pl2 = 1.50*cos2_ijk - 0.50
          Pl4 = (4.3750*cos2_ijk - 3.750)*cos2_ijk + 0.3750
          pl6 = ((14.43750*cos2_ijk - 19.68750)*cos2_ijk + 6.5625)*
     1          cos2_ijk - 0.31250

          Gi_sum0 = Gi_sum0 + fsij_fsik
          Gi_sum1 = Gi_sum1 + cos_ijk*fsij_fsik
          Gi_sum2 = Gi_sum2 + Pl2*fsij_fsik
          Gi_sum3 = Gi_sum3 + Pl4*fsij_fsik
          Gi_sum4 = Gi_sum4 + Pl6*fsij_fsik

         enddo ! do nb2 = 1,maxn_nbrs !
        enddo ! do nb1 = 1,maxn_nbrs !

        ! Calc. global G-vector !
        icol = n_set*5 - 4
        Gi_dev(icol,ni) = Gi_sum0
        Gi_dev(icol+1,ni) = Gi_sum1
        Gi_dev(icol+2,ni) = Gi_sum2
        Gi_dev(icol+3,ni) = Gi_sum3
        Gi_dev(icol+4,ni) = Gi_sum4

       enddo  ! do n_set=1,n_set_ann !
      enddo ! do ni=1,natoms !

!$OMP END PARALLEL DO

!
! --- Start III loop over ni atoms ---
!

!$OMP  PARALLEL DO DEFAULT(NONE)
!$OMP& PRIVATE(rij1, dcos_ijk, dgij, Pl2,Pl4,Pl6, dPl2,dPl4,dPl6)
!$OMP& PRIVATE(fsij,fsik, cos_ijk,cos2_ijk)
!$OMP& PRIVATE(Gi_fact0, Gi_fact1, Gi_fact2, Gi_fact3, Gi_fact4, icol)
!$OMP& PRIVATE(dgij1,dgij2,dgij3)
!$OMP& PRIVATE(Gi_sum_01,Gi_sum_11,Gi_sum_21,Gi_sum_31,Gi_sum_41)
!$OMP& PRIVATE(Gi_sum_02,Gi_sum_12,Gi_sum_22,Gi_sum_32,Gi_sum_42)
!$OMP& PRIVATE(Gi_sum_03,Gi_sum_13,Gi_sum_23,Gi_sum_33,Gi_sum_43)
!$OMP& SHARED(natoms, maxn_nbrs,n_set_ann, fsij_dev, Gi_dev)
!$OMP& SHARED(iPOT_file_ver)
!$OMP& SHARED(xr_ij0, xr_ij1, xr_ij2, xr_ij3)
!$OMP& SHARED(dGi_dx, dGi_dy, dGi_dz)
!$OMP& SHARED(dfs_rij_3D1, dfs_rij_3D2, dfs_rij_3D3)
!$OMP& SCHEDULE(STATIC)

      do ni = 1, natoms   ! loop over atoms !
       !DIR$ ASSUME (mod(maxn_nbrs,8) .eq. 0)
       do nb1 = 1,maxn_nbrs     ! Loop I: i - j bond !
         rij1 = 1.0/xr_ij0(nb1,ni)
         do n_set = 1,n_set_ann

         Gi_sum_01 = 0.0; Gi_sum_02 = 0.0; Gi_sum_03 = 0.0
         Gi_sum_11 = 0.0; Gi_sum_12 = 0.0; Gi_sum_13 = 0.0
         Gi_sum_21 = 0.0; Gi_sum_22 = 0.0; Gi_sum_23 = 0.0
         Gi_sum_31 = 0.0; Gi_sum_32 = 0.0; Gi_sum_33 = 0.0
         Gi_sum_41 = 0.0; Gi_sum_42 = 0.0; Gi_sum_43 = 0.0

         fsij = fsij_dev(nb1,n_set,ni)

        !DIR$ VECTOR UNALIGNED
        !DIR$ ASSUME (mod(maxn_nbrs,8) .eq. 0)
!$OMP simd simdlen(8)
        do nb2 = 1,maxn_nbrs     ! Loop I: i - j bond !

          fsik = 2.0*fsij_dev(nb2,n_set,ni)

          cos_ijk = xr_ij1(nb1,ni)*xr_ij1(nb2,ni) + 
     +    xr_ij2(nb1,ni)*xr_ij2(nb2,ni) + xr_ij3(nb1,ni)*xr_ij3(nb2,ni)
          cos2_ijk = cos_ijk**2

          dcos_ijk(1) = (xr_ij1(nb2,ni) - xr_ij1(nb1,ni)*cos_ijk)*rij1
          dcos_ijk(2) = (xr_ij2(nb2,ni) - xr_ij2(nb1,ni)*cos_ijk)*rij1
          dcos_ijk(3) = (xr_ij3(nb2,ni) - xr_ij3(nb1,ni)*cos_ijk)*rij1

          Gi_sum_01 = Gi_sum_01 + fsik*dfs_rij_3D1(n_set,nb1,ni)
          Gi_sum_02 = Gi_sum_02 + fsik*dfs_rij_3D2(n_set,nb1,ni)
          Gi_sum_03 = Gi_sum_03 + fsik*dfs_rij_3D3(n_set,nb1,ni)

          dgij1 = cos_ijk*dfs_rij_3D1(n_set,nb1,ni) + fsij*dcos_ijk(1)
          dgij2 = cos_ijk*dfs_rij_3D2(n_set,nb1,ni) + fsij*dcos_ijk(2)
          dgij3 = cos_ijk*dfs_rij_3D3(n_set,nb1,ni) + fsij*dcos_ijk(3)
          Gi_sum_11 = Gi_sum_11 + fsik*dgij1
          Gi_sum_12 = Gi_sum_12 + fsik*dgij2
          Gi_sum_13 = Gi_sum_13 + fsik*dgij3

          Pl2 = (1.50*cos2_ijk - 0.50); dPl2 = 3.0*cos_ijk
          dgij1 = Pl2*dfs_rij_3D1(n_set,nb1,ni) + fsij*dPl2*dcos_ijk(1)
          dgij2 = Pl2*dfs_rij_3D2(n_set,nb1,ni) + fsij*dPl2*dcos_ijk(2)
          dgij3 = Pl2*dfs_rij_3D3(n_set,nb1,ni) + fsij*dPl2*dcos_ijk(3)
          Gi_sum_21 = Gi_sum_21 + fsik*dgij1
          Gi_sum_22 = Gi_sum_22 + fsik*dgij2
          Gi_sum_23 = Gi_sum_23 + fsik*dgij3

          Pl4 = (4.3750*cos2_ijk - 3.750)*cos2_ijk + 0.3750
          dPl4 = (17.50*cos2_ijk - 7.50)*cos_ijk
          dgij1 = Pl4*dfs_rij_3D1(n_set,nb1,ni) + fsij*dPl4*dcos_ijk(1)
          dgij2 = Pl4*dfs_rij_3D2(n_set,nb1,ni) + fsij*dPl4*dcos_ijk(2)
          dgij3 = Pl4*dfs_rij_3D3(n_set,nb1,ni) + fsij*dPl4*dcos_ijk(3)
          Gi_sum_31 = Gi_sum_31 + fsik*dgij1
          Gi_sum_32 = Gi_sum_32 + fsik*dgij2
          Gi_sum_33 = Gi_sum_33 + fsik*dgij3

          Pl6 =((14.43750*cos2_ijk - 19.68750)*cos2_ijk + 6.5625)*
     &         cos2_ijk - 0.31250
          dPl6=(86.6250*cos2_ijk**2 - 78.750*cos2_ijk + 13.1250)*
     &          cos_ijk
          dgij1 = Pl6*dfs_rij_3D1(n_set,nb1,ni) + fsij*dPl6*dcos_ijk(1)
          dgij2 = Pl6*dfs_rij_3D2(n_set,nb1,ni) + fsij*dPl6*dcos_ijk(2)
          dgij3 = Pl6*dfs_rij_3D3(n_set,nb1,ni) + fsij*dPl6*dcos_ijk(3)
          Gi_sum_41 = Gi_sum_41 + fsik*dgij1
          Gi_sum_42 = Gi_sum_42 + fsik*dgij2
          Gi_sum_43 = Gi_sum_43 + fsik*dgij3

        enddo ! do nb2 = 1,maxn_nbrs !
!$OMP end simd

        icol = n_set*5 - 4

        Gi_fact0 = 1.0/sqrt(Gi_dev(icol,ni)**2 + 1.0)
        Gi_fact1 = 1.0/sqrt(Gi_dev(icol+1,ni)**2 + 1.0)
        Gi_fact2 = 1.0/sqrt(Gi_dev(icol+2,ni)**2 + 1.0)
        Gi_fact3 = 1.0/sqrt(Gi_dev(icol+3,ni)**2 + 1.0)
        Gi_fact4 = 1.0/sqrt(Gi_dev(icol+4,ni)**2 + 1.0)

        dGi_dx(icol,nb1,ni) = Gi_sum_01*Gi_fact0
        dGi_dy(icol,nb1,ni) = Gi_sum_02*Gi_fact0
        dGi_dz(icol,nb1,ni) = Gi_sum_03*Gi_fact0

        dGi_dx(icol+1,nb1,ni) = Gi_sum_11*Gi_fact1
        dGi_dy(icol+1,nb1,ni) = Gi_sum_12*Gi_fact1
        dGi_dz(icol+1,nb1,ni) = Gi_sum_13*Gi_fact1

        dGi_dx(icol+2,nb1,ni) = Gi_sum_21*Gi_fact2
        dGi_dy(icol+2,nb1,ni) = Gi_sum_22*Gi_fact2
        dGi_dz(icol+2,nb1,ni) = Gi_sum_23*Gi_fact2

        dGi_dx(icol+3,nb1,ni) = Gi_sum_31*Gi_fact3
        dGi_dy(icol+3,nb1,ni) = Gi_sum_32*Gi_fact3
        dGi_dz(icol+3,nb1,ni) = Gi_sum_33*Gi_fact3

        dGi_dx(icol+4,nb1,ni) = Gi_sum_41*Gi_fact4
        dGi_dy(icol+4,nb1,ni) = Gi_sum_42*Gi_fact4
        dGi_dz(icol+4,nb1,ni) = Gi_sum_43*Gi_fact4

         enddo ! do n_set = 1,n_set_ann !
       enddo ! do nb1 = 1,maxn_nbrs !
      enddo ! do ni=1,natoms !

!$OMP END PARALLEL DO

! --- Gis are done here for atom ni ---
! --- Start NN on atom ni ---

C     write(300+mynod,22) Gi_dev(1:8,1), ActFunc_shift
C     write(300+mynod,22) Gi_dev(9:16,1)
C     write(300+mynod,22) Gi_dev(17:24,1)
C     write(300+mynod,22) Gi_dev(25:32,1)
C     write(300+mynod,22) Gi_dev(33:40,1)
C     write(300+mynod,24)B1_ann(1:8)
C 22  format('dANN: Gi_dev(1)=',8f12.8,' Af=',f12.8)
C 24  format('dANN: B1_ann(1)=',8f12.8)

!$OMP  PARALLEL DO DEFAULT(NONE)
!$OMP& PRIVATE(U0,U1,U2, U_vect, NLcurr,NLnext, expU, TrFunc)
!$OMP& SHARED(natoms, mG_dim, Gi_dev, Nodes_of_layer, net_layers)
!$OMP& SHARED(B1_ann,B2_ann,B3_ann, W1_ann,W2_ann,W3_ann, Ep_of, U3_of)
!$OMP& SHARED(ActFunc_shift, dfuN_dev, base_pot_param)
!$OMP& SCHEDULE(STATIC)

      do ni=1,natoms
      
! -- Input Layer ---
       do iraw=1,Nodes_of_layer(2)
        U_vect = B1_ann(iraw)
        do icol=1,Nodes_of_layer(1)
         U_vect = U_vect + 
     + log(Gi_dev(icol,ni)+sqrt(Gi_dev(icol,ni)**2+1.0))*
     & W1_ann(icol,iraw)
        enddo
        TrFunc = 1.0 / (1.0 + exp(-U_vect))
        U2(iraw) = TrFunc + ActFunc_shift
        dfuN_dev(iraw,1,ni) = TrFunc*(1.0-TrFunc)
       enddo

C     if(ident(ni).eq.1) write(300+mynod,23) U2(1:8)
C 23  format('dANN: U2(1)=',8f12.8)

! -- Hidden Layers ---
       do layer=2, net_layers-2
        NLcurr = Nodes_of_layer(layer)
        NLnext = Nodes_of_layer(layer+1)
!       U1(1:NLcurr)=U2(1:NLcurr)
        call move_alloc(from=U1, to=U0)  ! swap U1 and U2 !
        call move_alloc(from=U2, to=U1)
        call move_alloc(from=U0, to=U2)
        do iraw=1,NLnext
         U_vect = B2_ann(iraw,layer)
         do icol=1,NLcurr
          U_vect = U_vect + U1(icol)*W2_ann(icol,iraw,layer)
         enddo
         TrFunc = 1.0 / (1.0 + exp(-U_vect))
         U2(iraw) = TrFunc + ActFunc_shift
         dfuN_dev(iraw,layer,ni) = TrFunc*(1.0-TrFunc)
        enddo
       enddo ! do layer=2, net_layers-1 !

! -- Output Layer ---     
       do iraw=1,Nodes_of_layer(net_layers)  ! 1.. 8 !
        U_vect = B3_ann(iraw)
        do icol=1,Nodes_of_layer(net_layers-1)   ! 1.. 20 !
         U_vect = U_vect + U2(icol)*W3_ann(icol,iraw)
        enddo
        U3_of(iraw,ni) = U_vect + base_pot_param(iraw)
       enddo

      enddo ! do ni=1,natoms !

!$OMP END PARALLEL DO

C     write(300+mynod,21) U3_of(1:8,1)
C 21  format('dANN: U3_of(1)=',8f12.8)

!$OMP  PARALLEL DO DEFAULT(NONE)
!$OMP& PRIVATE(U0f1,U0f2,U0f3, U1f1,U1f2,U1f3, U2f1,U2f2,U2f3)
!$OMP& PRIVATE(U_x1,U_x2,U_x3, NLcurr,NLnext)
!$OMP& SHARED(natoms, maxn_nbrs, Nodes_of_layer, net_layers)
!$OMP& SHARED(B1_ann,B2_ann,B3_ann, W1_ann,W2_ann,W3_ann, Ep_of)
!$OMP& SHARED(dfuN_dev,dANN_dxij)
!$OMP& SHARED(dGi_dx, dGi_dy, dGi_dz)
!$OMP& SCHEDULE(STATIC)

      do ni=1,natoms
       !DIR$ ASSUME (mod(maxn_nbrs,8) .eq. 0)
       do nb1 = 1,maxn_nbrs

! --- DO ANN for each (i-j) pair using Gis as input vectors ---

! -- Input Layer ---     
       do iraw=1,Nodes_of_layer(2)
        U_x1 = 0.0; U_x2 = 0.0; U_x3 = 0.0
        !DIR$ VECTOR UNALIGNED
        do icol=1,Nodes_of_layer(1)
         U_x1 = U_x1 + dGi_dx(icol,nb1,ni)*W1_ann(icol,iraw)
         U_x2 = U_x2 + dGi_dy(icol,nb1,ni)*W1_ann(icol,iraw)
         U_x3 = U_x3 + dGi_dz(icol,nb1,ni)*W1_ann(icol,iraw)
        enddo
        U2f1(iraw) = U_x1*dfuN_dev(iraw,1,ni)
        U2f2(iraw) = U_x2*dfuN_dev(iraw,1,ni)
        U2f3(iraw) = U_x3*dfuN_dev(iraw,1,ni)
       enddo                  ! exp(-u)/(1+exp(-u))^2 !

! -- Hidden Layers ---     
       do layer=2, net_layers-2
        NLcurr = Nodes_of_layer(layer)
        NLnext = Nodes_of_layer(layer+1)

!       !DIR$ VECTOR UNALIGNED
!       do i=1,NLcurr
!        U1f1(i)=U2f1(i)
!        U1f2(i)=U2f2(i)
!        U1f3(i)=U2f3(i)
!       enddo

        call move_alloc(from=U1f1, to=U0f1)  ! swap U1f and U2f !
        call move_alloc(from=U2f1, to=U1f1) 
        call move_alloc(from=U0f1, to=U2f1)

        call move_alloc(from=U1f2, to=U0f2)  
        call move_alloc(from=U2f2, to=U1f2)  
        call move_alloc(from=U0f2, to=U2f2)  

        call move_alloc(from=U1f3, to=U0f3)  
        call move_alloc(from=U2f3, to=U1f3)  
        call move_alloc(from=U0f3, to=U2f3)  

        !DIR$ CODE_ALIGN :64
        do iraw=1,NLnext 
          U_x1 = 0.0; U_x2 = 0.0; U_x3 = 0.0
          !DIR$ VECTOR UNALIGNED
          do icol=1,NLcurr
            U_x1 = U_x1 + U1f1(icol)*W2_ann(icol,iraw,layer)
            U_x2 = U_x2 + U1f2(icol)*W2_ann(icol,iraw,layer)
            U_x3 = U_x3 + U1f3(icol)*W2_ann(icol,iraw,layer)
          enddo
          U2f1(iraw) = U_x1*dfuN_dev(iraw,layer,ni)
          U2f2(iraw) = U_x2*dfuN_dev(iraw,layer,ni)
          U2f3(iraw) = U_x3*dfuN_dev(iraw,layer,ni)
        enddo               ! exp(-u)/(1+exp(-u))^2
       enddo ! do layer=2, net_layers-1 !

! -- Output Layer ---

       do iraw=1,Nodes_of_layer(net_layers)  ! 1.. 1 !
        U_x1 = 0.0; U_x2 = 0.0; U_x3 = 0.0
        !DIR$ VECTOR UNALIGNED
        do icol=1,Nodes_of_layer(net_layers-1)   ! 1.. 20 !
         U_x1 = U_x1 + U2f1(icol)*W3_ann(icol,iraw)
         U_x2 = U_x2 + U2f2(icol)*W3_ann(icol,iraw)
         U_x3 = U_x3 + U2f3(icol)*W3_ann(icol,iraw)
        enddo
        dANN_dxij(1,iraw,nb1,ni) = U_x1
        dANN_dxij(2,iraw,nb1,ni) = U_x2
        dANN_dxij(3,iraw,nb1,ni) = U_x3
       enddo

! --- End of ANN for each (i-j) pair using gij as input vectors ---

       enddo ! do nb1 = 1,maxn_nbrs !
      enddo ! do ni = 1, natoms !

!$OMP END PARALLEL DO

      return
      end subroutine      ! dANN_OMP  !
!
! ---------------------------------------------------------------------
! Calculates Analytical derivatives and force calculation.
! ---------------------------------------------------------------------
!
      subroutine Frc_ANN_OMP(ecohe)

      use atoms

      implicit double precision (a-h,o-z)
  
      double precision, intent(out) :: ecohe

      double precision, dimension(3) :: fij,fr

      if(ihalt.ne.0) return

      call dANN_OMP(max_nbrs, Rc_ann)  ! Get dANN_dxij() !
 
! --- Calc Actual Force Vectors and Energy ---

      ecohe = 0.0

!$OMP  PARALLEL DO DEFAULT(NONE) REDUCTION(+:ecohe)
!$OMP& PRIVATE(Ep, fr, fij, nbi, nj,nj2, ij_delta, nbrs_of_j)
!$OMP& SHARED(natoms, max_nbrs, nbr_list)
!$OMP& SHARED(frr, Ep_of, U3_of, dANN_dxij)
!$OMP& SCHEDULE(STATIC)

      do ni = 1,natoms   ! loop over atoms !

       fr(1:3) = 0.0

       !DIR$ ASSUME (mod(max_nbrs,8) .eq. 0)
       do nb1 = 1,max_nbrs     ! Loop I: i - j bond !
         nj = nbr_list(nb1,ni)
         !call mm_prefetch(dANN_dxij(1,1,nb1,nj),1)
         !call mm_prefetch(dANN_dxij(1,1,nb2,nj),1)
         ij_delta = min(abs(nj-ni),1) ! 0 if ni=nj; 1 if ni=/=nj !
         nbrs_of_j = ij_delta*max_nbrs  ! rij < Rc !

         nbi=nb1
         !DIR$ VECTOR UNALIGNED
         !DIR$ ASSUME (nbrs_of_j.eq.0 .or. mod(nbrs_of_j,8).eq.0)
         do nb2 = 1,nbrs_of_j  ! search for i as a neighbor of j !
          nj2 = nbr_list(nb2,nj)
          if(nj2.eq.ni) nbi=nb2       ! i is the nbi-th nbr of j !
         enddo
                                   ! Y3ij - Y3ji !
         fij(:)=dANN_dxij(:,1,nb1,ni)-dANN_dxij(:,1,nbi,nj)
         fr(1:3) = fr(1:3) + fij(1:3)   ! x,y,z forces !
       enddo ! do nb1 = 1,max_nbrs !

       frr(1:3,ni) = fr(1:3)
   
       Ep = U3_of(1,ni)
       Ep_of(ni) = 2.0*Ep  ! Twice the atomic energy  !
       ecohe = ecohe + Ep   ! devided by 2 later in MSR for !
                            ! compatibility with other potentials !
      enddo ! do ni = 1,natoms !

!$OMP END PARALLEL DO

      return
      end subroutine      ! Frc_ANN_OMP  !
!
! ---------------------------------------------------------------------
!
      subroutine alloc_types_ANN(ierror)

      integer, intent(out) :: ierror
      integer :: ialloc(8)

      ialloc(:) = 0

      max_cols = 1; max_raws = 1
      do i=2,net_layers-2
       if(Nodes_of_layer(i).gt.max_cols) max_cols=Nodes_of_layer(i)
       if(Nodes_of_layer(i+1).gt.max_raws) max_raws=Nodes_of_layer(i+1)
      enddo

      Ncols1 = Nodes_of_layer(1)
      Nraws1 = Nodes_of_layer(2)
      Ncols2 = max_cols
      Nraws2 = max_raws
      Ncols3 = Nodes_of_layer(net_layers-1)
      Nraws3 = Nodes_of_layer(net_layers)

      nbuf_dim = Ncols1*Nraws1 + (net_layers-2)*max_cols*max_raws + 
     +           Ncols3*Nraws3

      allocate(W1_ann(Ncols1,Nraws1), stat=ialloc(1))
      allocate(W2_ann(Ncols2, Nraws2, 2:net_layers-2), stat=ialloc(2))
      allocate(W3_ann(Ncols3,Nraws3), stat=ialloc(3))

      allocate(B1_ann(Nraws1), stat=ialloc(4))
      allocate(B2_ann(Nraws2, 2:net_layers-2), stat=ialloc(5))
      allocate(B3_ann(Nraws3), stat=ialloc(6))
      allocate(buf_ann(nbuf_dim), stat=ialloc(7))

      ierror = 0
      do i=1,7
       ierror = ierror + ialloc(i)
      enddo

      return
      end subroutine        ! alloc_types_ANN !
!
! ---------------------------------------------------------------------
!
      subroutine deall_types_ANN(ierror)

      implicit double precision (a-h,o-z)

      integer, intent(out) :: ierror
      integer ialloc(10)

      ialloc(:) = 0

      if (allocated(W1_ann)) deallocate(W1_ann, stat=ialloc(1))
      if (allocated(W2_ann)) deallocate(W2_ann, stat=ialloc(2))
      if (allocated(W3_ann)) deallocate(W3_ann, stat=ialloc(3))

      if (allocated(B1_ann)) deallocate(B1_ann, stat=ialloc(4))
      if (allocated(B2_ann)) deallocate(B2_ann, stat=ialloc(5))
      if (allocated(B3_ann)) deallocate(B3_ann, stat=ialloc(6))

      if (allocated(r0_value)) deallocate(r0_value, stat=ialloc(7))
      if (allocated(r0G_value)) deallocate(r0G_value, stat=ialloc(8))
      if (allocated(buf_ann)) deallocate(buf_ann, stat=ialloc(9))

      if (allocated(base_pot_param)) deallocate(base_pot_param,
     1 stat=ialloc(10))

      ierror = 0
      do i=1,10
       ierror = ierror + ialloc(i)
      enddo

      return
      end subroutine        ! deall_types_ANN !
!
! ---------------------------------------------------------------------
!
      subroutine alloc_atoms_ANN(ierror)

      integer, intent(out) :: ierror
      integer ialloc(17)

      ialloc(:) = 0

      max_ls = 1
      do i=1,net_layers
       if (Nodes_of_layer(i).gt.max_ls) max_ls=Nodes_of_layer(i)
      enddo

      if(I_have_GPU.eq.0) then
       allocate(dBOP_param_dxij(3, nBOP_params, nbrs_per_atom,
     & natoms_alloc), stat=ialloc(1))
      endif

      allocate(Gi_atom(0:4,n_set_ann,max_tri_index), stat=ialloc(3))
      allocate(dG_i(0:4,n_set_ann,max_tri_index), stat=ialloc(4))
      allocate(Gi_list(mG_dim,natoms_alloc), stat=ialloc(5))
      allocate(Gi_new(mG_dim,nbrs_per_atom), stat=ialloc(6))

      allocate(U1(max_ls), stat=ialloc(7))
      allocate(U2(max_ls), stat=ialloc(8))
      allocate(U3(nBOP_params), stat=ialloc(8))
      allocate(U3_of(nBOP_params, natoms_alloc), stat=ialloc(8))

      allocate(xr_ij0(nbrs_per_atom,natoms_alloc), stat=ialloc(9))
      allocate(xr_ij1(nbrs_per_atom,natoms_alloc), stat=ialloc(9))
      allocate(xr_ij2(nbrs_per_atom,natoms_alloc), stat=ialloc(9))
      allocate(xr_ij3(nbrs_per_atom,natoms_alloc), stat=ialloc(9))
      allocate(fsij_dev(nbrs_per_atom,n_set_ann,natoms_alloc), 
     & stat=ialloc(10))
      allocate(dfs_rij_3D1(n_set_ann,nbrs_per_atom,natoms_alloc), 
     & stat=ialloc(11))
      allocate(dfs_rij_3D2(n_set_ann,nbrs_per_atom,natoms_alloc), 
     & stat=ialloc(11))
      allocate(dfs_rij_3D3(n_set_ann,nbrs_per_atom,natoms_alloc),
     & stat=ialloc(11))
      allocate(Gi_dev(mG_dim,natoms_alloc), stat=ialloc(12))
      allocate(dGi_dx(mG_dim,nbrs_per_atom,natoms_alloc), 
     & stat=ialloc(13))
      allocate(dGi_dy(mG_dim,nbrs_per_atom,natoms_alloc), 
     & stat=ialloc(13))
      allocate(dGi_dz(mG_dim,nbrs_per_atom,natoms_alloc), 
     & stat=ialloc(13))
      allocate(dfuN_dev(max_ls,net_layers-2,natoms_alloc), 
     & stat=ialloc(14))
      allocate(U1f1(max_ls), U1f2(max_ls),U1f3(max_ls), stat=ialloc(15))
      allocate(U2f1(max_ls), U2f2(max_ls),U2f3(max_ls), stat=ialloc(15))
      allocate(dANN_dxij(3,nBOP_params,nbrs_per_atom,natoms_alloc), 
     & stat=ialloc(16))
      allocate(r0Rc(n_set_ann), r0pRc(n_set_ann), stat=ialloc(17))

      ierror = 0
      do i=1,17
       ierror = ierror + ialloc(i)
      enddo
      if(ierror.ne.0) then
       write(6,*)'ERROR allocating x in alloc_atoms_ANN'
       ihalt = 1
      endif

      return
      end subroutine        ! alloc_atoms_ANN !
!
! ---------------------------------------------------------------------
!
      subroutine deall_atoms_ANN(ierror)

      integer, intent(out) :: ierror
      integer ialloc(17)

      ialloc(:) = 0

      if (allocated(dBOP_param_dxij)) deallocate(dBOP_param_dxij, 
     & stat=ialloc(2))

      if (allocated(Gi_atom)) deallocate(Gi_atom, stat=ialloc(3))
      if (allocated(dG_i)) deallocate(dG_i, stat=ialloc(4))
      if (allocated(Gi_list)) deallocate(Gi_list, stat=ialloc(5))

      if (allocated(U1)) deallocate(U1, stat=ialloc(6))
      if (allocated(U2)) deallocate(U2, stat=ialloc(7))
      if (allocated(U3)) deallocate(U3, stat=ialloc(7))
      if (allocated(U3_of)) deallocate(U3_of, stat=ialloc(7))

      if (allocated(Gi_new)) deallocate(Gi_new, stat=ialloc(8))

      if (allocated(xr_ij0)) deallocate(xr_ij0, stat=ialloc(9))
      if (allocated(xr_ij1)) deallocate(xr_ij1, stat=ialloc(9))
      if (allocated(xr_ij2)) deallocate(xr_ij2, stat=ialloc(9))
      if (allocated(xr_ij3)) deallocate(xr_ij3, stat=ialloc(9))
      if (allocated(fsij_dev)) deallocate(fsij_dev, stat=ialloc(10))
      if (allocated(dfs_rij_3D1))deallocate(dfs_rij_3D1,stat=ialloc(11))
      if (allocated(dfs_rij_3D2))deallocate(dfs_rij_3D2,stat=ialloc(11))
      if (allocated(dfs_rij_3D3))deallocate(dfs_rij_3D3,stat=ialloc(11))
      if (allocated(Gi_dev)) deallocate(Gi_dev, stat=ialloc(12))
      if (allocated(dGi_dx)) deallocate(dGi_dx, stat=ialloc(13))
      if (allocated(dGi_dy)) deallocate(dGi_dy, stat=ialloc(13))
      if (allocated(dGi_dz)) deallocate(dGi_dz, stat=ialloc(13))
      if (allocated(dfuN_dev)) deallocate(dfuN_dev, stat=ialloc(14))
      if (allocated(U1f1)) deallocate(U1f1,U1f2,U1f3, stat=ialloc(15))
      if (allocated(dANN_dxij)) deallocate(dANN_dxij, stat=ialloc(16))
      if (allocated(r0Rc)) deallocate(r0Rc,r0pRc, stat=ialloc(17))

      ierror = 0
      do i=1,17
       ierror = ierror + ialloc(i)
      enddo

      return
      end subroutine        ! deall_atoms_ANN !
!
! ------------------------------------------------------------------
!
      END MODULE  ! BOP !
