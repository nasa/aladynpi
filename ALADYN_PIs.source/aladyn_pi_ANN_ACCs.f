!
!------------------------------------------------------------------
! 11-20-2020
!
! Artificial Neural Network module for aladyn_pi.f code
! Open_ACC version
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
      MODULE ANN_ACC

      use sys_OMP
      use constants
      use sim_box
      use pot_module
      use ANN
      save

      real(kind=4), dimension(:,:), allocatable :: r_ij_dev
      real(kind=4), dimension(:,:), allocatable :: x_ij_dev
      real(kind=4), dimension(:,:), allocatable :: y_ij_dev
      real(kind=4), dimension(:,:), allocatable :: z_ij_dev
      real(kind=4), dimension(:,:), allocatable :: U1f,U2f,Gi_acc
      real(kind=4), dimension(:,:,:), allocatable :: 
     1 fsij_acc, dfuN_acc
      real(kind=4), dimension(:,:,:), allocatable :: dfs_dxij
      real(kind=4), dimension(:,:,:), allocatable :: dfs_dyij
      real(kind=4), dimension(:,:,:), allocatable :: dfs_dzij
      real(kind=4), dimension(:,:,:,:), allocatable :: Gi_3D_acc


      CONTAINS
!
! ---------------------------------------------------------------------
! Calculates ANN and its Analytical derivatives for the 
! force calculation using OpenACC.
! ---------------------------------------------------------------------
!
      subroutine dANN_ACC(maxn_nbrs, Rcn)

      use sys_ACC
      use atoms

      implicit real(kind=4) (a-h,o-z)

      integer, intent(in) :: maxn_nbrs
      real(kind=4), intent(in) :: Rcn
  
      integer (kind=8) :: iMb = 1024**2
      integer (kind=8) :: memory_request, memory1
 
      real(kind=4), dimension(3) :: dgij
      real(kind=4), dimension(3) :: xr_ij, xr_ik, dcos_ijk
      real(kind=4), dimension(n_set_ann) :: r0Rc, r0pRc

C     write(6,*)'dANN_ACC: natoms=',natoms,' maxn_nbrs=',maxn_nbrs

      max_ls = 1
      do i=2,net_layers-1
       if (Nodes_of_layer(i).gt.max_ls) max_ls=Nodes_of_layer(i)
      enddo

      h11 = h(1,1); h12 = h(1,2); h13 = h(1,3)
      h22 = h(2,2); h23 = h(2,3); h33 = h(3,3)

      ! Those are replacements of ACC_* equivalents    !
      ! redefined in pgmc_sys_ACC.f and pgmc_sys_OMP.f !

      My_GPU_free_mem = GET_GPU_FREE_MEM(My_GPU) / iMb  ! in Mbs !
!     write(6,*)'1: My_GPU_free_mem=',My_GPU_free_mem

       memory1 = 4*(2*maxn_nbrs+n_set_ann*maxn_nbrs) +
     + mG_dim + 3*(mG_dim+nBOP_params)*maxn_nbrs + max_ls*net_layers

       memory_request = memory1*natoms*4  ! kind=4 !

       memory_request = memory_request / iMb ! in Mbs !
       memory_request = memory_request * 1.3 ! safe margin !

       if(nACC_devices.gt.0) then
       if(memory_request.gt.My_GPU_free_mem) then
        write(6,13) memory_request, My_GPU_free_mem, My_GPU, mynod
        ihalt = 1; return
       endif
       endif
  13   format(/,'ERROR: Requested memory in dANN_ACC is:',i6,
     1 ' Mb,',/,'which exceeds the available GPU memory of',i6,
     2 ' Mb on GPU:',i2,' on node:',i5,/,'Increase number of nodes.')

      Rc = Rcn; Rc2 = Rcn**2; d4 = d4_ann
      Sigma2 = 2.0*Gauss_ann**2
      r0Rc(:) = Rcn*r0_value(:)
      r0pRc(:) = Rcn + r0_value(:)

!
! --- Start I loop over ni atoms ---
!
!$ACC PARALLEL LOOP COLLAPSE (2)
!$ACC& PRIVATE(xr_ij)
      do ni=1,natoms   ! loop over atoms !
       do nb1 = 1,maxn_nbrs     ! Loop I: i - j bond !

        nj = nbr_list(nb1,ni)
        ij_delta = 1-min(abs(nj-ni),1) ! 1 if ni=nj; 0 if ni=/=nj !

        sxn = sx(ni); syn = sy(ni); szn = sz(ni)

        sx0 = sx(nj) - sxn
        sx0 = sx0 - nint(sx0)     ! make periodic along X !
        sy0 = sy(nj) - syn
        sy0 = sy0 - nint(sy0)     ! make periodic along Y !
        sz0 = sz(nj) - szn
        sz0 = sz0 - nint(sz0)     ! make periodic along Z !

        xij = h11*sx0 + h12*sy0 + h13*sz0
        yij = h22*sy0 + h23*sz0
        zij = h33*sz0
        r2ij = xij**2 + yij**2 + zij**2 + ij_delta*Rc2 ! rij+Rc if i=j
        rij = sqrt(r2ij)
        r1ij = 1.0/rij

        xr_ij(1) = xij*r1ij
        xr_ij(2) = yij*r1ij
        xr_ij(3) = zij*r1ij
        r_ij_dev(nb1,ni) = rij
        x_ij_dev(nb1,ni) = xr_ij(1) ! x/r !
        y_ij_dev(nb1,ni) = xr_ij(2) ! y/r !
        z_ij_dev(nb1,ni) = xr_ij(3) ! z/r !

        RcRij = max(Rc-rij, 0.0)
        RcRij3 = RcRij**3
        RcRij4 = RcRij3*RcRij
        RcRij5 = RcRij4*RcRij
        fc_rij = RcRij4 / (d4+RcRij4)
        denominator = (Gauss_ann * (d4+RcRij4))**2

!$ACC LOOP SEQ
        do n_set = 1,n_set_ann
         r0 = r0_value(n_set)  ! do not use r0G_value() here !
         RijRo = rij-r0
         expRijRo = exp(-(RijRo/Gauss_ann)**2)
         dfsij = (r0Rc(n_set) + rij*(rij-r0pRc(n_set)) - Sigma2)*d4 - 
     -            RijRo*RcRij5
         dfsij = 2.0*dfsij*RcRij3*expRijRo / (r0*denominator)
         dfs_dxij(nb1,ni,n_set) = dfsij*xr_ij(1) ! dfs_ij*xij/rij !
         dfs_dyij(nb1,ni,n_set) = dfsij*xr_ij(2) ! dfs_ij*xij/rij !
         dfs_dzij(nb1,ni,n_set) = dfsij*xr_ij(3) ! dfs_ij*xij/rij !
         fsij_acc(nb1,n_set,ni) = fc_rij * expRijRo / r0
        enddo ! do n_set = 1,n_set_ann !

       enddo ! do nb1 = 1,maxn_nbrs !
      enddo ! do ni=1,natoms !
!
! --- Start II loop over ni atoms ---
!
!$ACC PARALLEL LOOP COLLAPSE(2)
!$ACC& PRIVATE(Gi_sum0, Gi_sum1, Gi_sum2, Gi_sum3, Gi_sum4)
      do ni=1,natoms   ! loop over atoms !
      do n_set = 1,n_set_ann

       Gi_sum0 = 0.0; Gi_sum1 = 0.0; Gi_sum2 = 0.0
       Gi_sum3 = 0.0; Gi_sum4 = 0.0

!$ACC LOOP REDUCTION(+:Gi_sum0, Gi_sum1, Gi_sum2, Gi_sum3, Gi_sum4)
!  !$ACC& COLLAPSE (2)
!$ACC& PRIVATE(xr_ij, xr_ik, dcos_ijk, fsij_fsik)
       do nb1 = 1,maxn_nbrs     ! Loop I: i - j bond !
        do nb2 = 1,maxn_nbrs     ! Loop I: i - j bond !

         fsij_fsik = fsij_acc(nb1,n_set,ni) * fsij_acc(nb2,n_set,ni)

         xr_ij(1) = x_ij_dev(nb1,ni)
         xr_ij(2) = y_ij_dev(nb1,ni)
         xr_ij(3) = z_ij_dev(nb1,ni)

         xr_ik(1) = x_ij_dev(nb2,ni)
         xr_ik(2) = y_ij_dev(nb2,ni)
         xr_ik(3) = z_ij_dev(nb2,ni)

       cos_ijk = xr_ij(1)*xr_ik(1)+xr_ij(2)*xr_ik(2)+xr_ij(3)*xr_ik(3)
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
       Gi_acc(icol,ni) = Gi_sum0
       Gi_acc(icol+1,ni) = Gi_sum1
       Gi_acc(icol+2,ni) = Gi_sum2
       Gi_acc(icol+3,ni) = Gi_sum3
       Gi_acc(icol+4,ni) = Gi_sum4

       enddo  ! do n_set=1,n_set_ann !
      enddo ! do ni=1,natoms !
!
! --- Start III loop over ni atoms ---
!
!$ACC PARALLEL LOOP VECTOR_LENGTH(64)
!$ACC& PRIVATE(xr_ij, xr_ik, dcos_ijk, dgij, fsij, fsik,rij1,
!$ACC& Gi_sum0_x, Gi_sum0_y, Gi_sum0_z,
!$ACC& Gi_sum1_x, Gi_sum1_y, Gi_sum1_z,
!$ACC& Gi_sum2_x, Gi_sum2_y, Gi_sum2_z,
!$ACC& Gi_sum3_x, Gi_sum3_y, Gi_sum3_z,
!$ACC& Gi_sum4_x, Gi_sum4_y, Gi_sum4_z)
      do ni = 1, natoms   ! loop over atoms ! 
!$ACC LOOP COLLAPSE(2) VECTOR
       do nb1 = 1,maxn_nbrs     ! Loop I: i - j bond !
         do n_set = 1,n_set_ann

         Gi_sum0_x = 0.0; Gi_sum0_y = 0.0; Gi_sum0_z = 0.0
         Gi_sum1_x = 0.0; Gi_sum1_y = 0.0; Gi_sum1_z = 0.0
         Gi_sum2_x = 0.0; Gi_sum2_y = 0.0; Gi_sum2_z = 0.0
         Gi_sum3_x = 0.0; Gi_sum3_y = 0.0; Gi_sum3_z = 0.0
         Gi_sum4_x = 0.0; Gi_sum4_y = 0.0; Gi_sum4_z = 0.0

         rij = r_ij_dev(nb1,ni)
         xr_ij(1) = x_ij_dev(nb1,ni)
         xr_ij(2) = y_ij_dev(nb1,ni)
         xr_ij(3) = z_ij_dev(nb1,ni)
         rij1 = 1.0/rij

         fsij = fsij_acc(nb1,n_set,ni)

!$ACC LOOP SEQ
        do nb2 = 1,maxn_nbrs     ! Loop I: i - j bond !

         fsik = 2.0*fsij_acc(nb2,n_set,ni)

         xr_ik(1) = x_ij_dev(nb2,ni)
         xr_ik(2) = y_ij_dev(nb2,ni)
         xr_ik(3) = z_ij_dev(nb2,ni)

       cos_ijk = xr_ij(1)*xr_ik(1)+xr_ij(2)*xr_ik(2)+xr_ij(3)*xr_ik(3)
       cos2_ijk = cos_ijk**2
       dcos_ijk(1) = (xr_ik(1) - xr_ij(1)*cos_ijk)*rij1
       dcos_ijk(2) = (xr_ik(2) - xr_ij(2)*cos_ijk)*rij1
       dcos_ijk(3) = (xr_ik(3) - xr_ij(3)*cos_ijk)*rij1
                   ! (xik/rik - xij/rij*cos_ijk) / rij  !

      Gi_sum0_x = Gi_sum0_x + fsik*dfs_dxij(nb1,ni,n_set)
      Gi_sum0_y = Gi_sum0_y + fsik*dfs_dyij(nb1,ni,n_set)
      Gi_sum0_z = Gi_sum0_z + fsik*dfs_dzij(nb1,ni,n_set)

      dgij(1) = cos_ijk*dfs_dxij(nb1,ni,n_set) + fsij*dcos_ijk(1)
      dgij(2) = cos_ijk*dfs_dyij(nb1,ni,n_set) + fsij*dcos_ijk(2)
      dgij(3) = cos_ijk*dfs_dzij(nb1,ni,n_set) + fsij*dcos_ijk(3)
      Gi_sum1_x = Gi_sum1_x + fsik*dgij(1)
      Gi_sum1_y = Gi_sum1_y + fsik*dgij(2)
      Gi_sum1_z = Gi_sum1_z + fsik*dgij(3)

      Pl2 = (1.50*cos2_ijk - 0.50); dPl2 = 3.0*cos_ijk
      dgij(1) = Pl2*dfs_dxij(nb1,ni,n_set) + fsij*dPl2*dcos_ijk(1)
      dgij(2) = Pl2*dfs_dyij(nb1,ni,n_set) + fsij*dPl2*dcos_ijk(2)
      dgij(3) = Pl2*dfs_dzij(nb1,ni,n_set) + fsij*dPl2*dcos_ijk(3)
      Gi_sum2_x = Gi_sum2_x + fsik*dgij(1)
      Gi_sum2_y = Gi_sum2_y + fsik*dgij(2)
      Gi_sum2_z = Gi_sum2_z + fsik*dgij(3)

      Pl4 = (4.3750*cos2_ijk - 3.750)*cos2_ijk + 0.3750
      dPl4 = (17.50*cos2_ijk - 7.50)*cos_ijk
      dgij(1) = Pl4*dfs_dxij(nb1,ni,n_set) + fsij*dPl4*dcos_ijk(1)
      dgij(2) = Pl4*dfs_dyij(nb1,ni,n_set) + fsij*dPl4*dcos_ijk(2)
      dgij(3) = Pl4*dfs_dzij(nb1,ni,n_set) + fsij*dPl4*dcos_ijk(3)
      Gi_sum3_x = Gi_sum3_x + fsik*dgij(1)
      Gi_sum3_y = Gi_sum3_y + fsik*dgij(2)
      Gi_sum3_z = Gi_sum3_z + fsik*dgij(3)

      Pl6 =((14.43750*cos2_ijk - 19.68750)*cos2_ijk + 6.5625)*cos2_ijk
     -    - 0.31250
      dPl6=(86.6250*cos2_ijk**2 - 78.750*cos2_ijk + 13.1250)*cos_ijk
      dgij(1) = Pl6*dfs_dxij(nb1,ni,n_set) + fsij*dPl6*dcos_ijk(1)
      dgij(2) = Pl6*dfs_dyij(nb1,ni,n_set) + fsij*dPl6*dcos_ijk(2)
      dgij(3) = Pl6*dfs_dzij(nb1,ni,n_set) + fsij*dPl6*dcos_ijk(3)
      Gi_sum4_x = Gi_sum4_x + fsik*dgij(1)
      Gi_sum4_y = Gi_sum4_y + fsik*dgij(2)
      Gi_sum4_z = Gi_sum4_z + fsik*dgij(3)

        enddo ! do nb2 = 1,maxn_nbrs !

        icol = n_set*5 - 4

        Gi_fact0 = 1.0/sqrt(Gi_acc(icol,ni)**2 + 1.0)
        Gi_fact1 = 1.0/sqrt(Gi_acc(icol+1,ni)**2 + 1.0)
        Gi_fact2 = 1.0/sqrt(Gi_acc(icol+2,ni)**2 + 1.0)
        Gi_fact3 = 1.0/sqrt(Gi_acc(icol+3,ni)**2 + 1.0)
        Gi_fact4 = 1.0/sqrt(Gi_acc(icol+4,ni)**2 + 1.0)

        Gi_3D_acc(icol,nb1,ni,1) = Gi_sum0_x*Gi_fact0
        Gi_3D_acc(icol,nb1,ni,2) = Gi_sum0_y*Gi_fact0
        Gi_3D_acc(icol,nb1,ni,3) = Gi_sum0_z*Gi_fact0

        Gi_3D_acc(icol+1,nb1,ni,1) = Gi_sum1_x*Gi_fact1
        Gi_3D_acc(icol+1,nb1,ni,2) = Gi_sum1_y*Gi_fact1
        Gi_3D_acc(icol+1,nb1,ni,3) = Gi_sum1_z*Gi_fact1

        Gi_3D_acc(icol+2,nb1,ni,1) = Gi_sum2_x*Gi_fact2
        Gi_3D_acc(icol+2,nb1,ni,2) = Gi_sum2_y*Gi_fact2
        Gi_3D_acc(icol+2,nb1,ni,3) = Gi_sum2_z*Gi_fact2

        Gi_3D_acc(icol+3,nb1,ni,1) = Gi_sum3_x*Gi_fact3
        Gi_3D_acc(icol+3,nb1,ni,2) = Gi_sum3_y*Gi_fact3
        Gi_3D_acc(icol+3,nb1,ni,3) = Gi_sum3_z*Gi_fact3

        Gi_3D_acc(icol+4,nb1,ni,1) = Gi_sum4_x*Gi_fact4
        Gi_3D_acc(icol+4,nb1,ni,2) = Gi_sum4_y*Gi_fact4
        Gi_3D_acc(icol+4,nb1,ni,3) = Gi_sum4_z*Gi_fact4

         enddo ! do n_set = 1,n_set_ann !
       enddo ! do nb1 = 1,maxn_nbrs !
      enddo ! do ni=1,natoms !

! --- Gis are done here for atom ni ---
! --- Start NN on atom ni ---

!$ACC PARALLEL LOOP PRIVATE(U1,U2, U3_vect) VECTOR_LENGTH(16)
      do ni=1,natoms
       
! -- Input Layer ---
!$ACC LOOP PRIVATE(U_vect)
       do iraw=1,Nodes_of_layer(2)  ! 1.. 20 !
        U_vect = B1_ann(iraw)
!!$ACC LOOP REDUCTION(+:U_vect)
!$ACC LOOP SEQ
        do icol=1,Nodes_of_layer(1)   ! 1.. 60 !
         U_vect = U_vect + log(Gi_acc(icol,ni) +
     +   sqrt(Gi_acc(icol,ni)**2+1.0))*W1_ann(icol,iraw)
        enddo
        TrFunc = 1.0 / (1.0 + exp(-U_vect))
        U2(iraw) = TrFunc + ActFunc_shift
        dfuN_acc(iraw,1,ni) = TrFunc*(1.0-TrFunc)
       enddo

! -- Hidden Layers ---
!$ACC LOOP SEQ
       do layer=2, net_layers-2
        NLcurr = Nodes_of_layer(layer)
        NLnext = Nodes_of_layer(layer+1)
        U1(1:NLcurr)=U2(1:NLcurr)
!$ACC LOOP PRIVATE(U_vect)
        do iraw=1,NLnext  ! 1.. 20 !
         U_vect = B2_ann(iraw,layer)
!!$ACC LOOP REDUCTION(+:U_vect)
!$ACC LOOP SEQ
         do icol=1,NLcurr  ! 1.. 20 !
          U_vect = U_vect + U1(icol)*W2_ann(icol,iraw,layer)
         enddo
         TrFunc = 1.0 / (1.0 + exp(-U_vect))
         U2(iraw) = TrFunc + ActFunc_shift
         dfuN_acc(iraw,layer,ni) = TrFunc*(1.0-TrFunc)
        enddo
       enddo ! do layer=2, net_layers-1 !

! -- Output Layer ---     
!$ACC LOOP
       do iraw=1,Nodes_of_layer(net_layers)  ! 1.. 8 !
        U_vect = B3_ann(iraw)
!!$ACC LOOP REDUCTION(+:U_vect)
!$ACC LOOP SEQ
        do icol=1,Nodes_of_layer(net_layers-1)   ! 1.. 20 !
         U_vect = U_vect + U2(icol)*W3_ann(icol,iraw)
        enddo
        U3_of(iraw,ni) = U_vect + base_pot_param(iraw)
       enddo

      enddo ! do ni=1,natoms !

!$ACC PARALLEL LOOP GANG COLLAPSE(2) VECTOR_LENGTH(32)
!$ACC& PRIVATE(U1f,U2f,NLcurr,NLnext,U_x,U_y,U_z)
      do ni=1,natoms
       do nb1 = 1,maxn_nbrs

! --- DO ANN for each (i-j) pair using Gis as input vectors ---

! -- Input Layer ---     
!ACC LOOP VECTOR 
       do iraw=1,Nodes_of_layer(2)  ! 1.. 20 !
        U_x = 0.0; U_y = 0.0; U_z = 0.0
!$ACC LOOP SEQ 
        do icol=1,Nodes_of_layer(1)   ! 1.. 60 !
         U_x = U_x + Gi_3D_acc(icol,nb1,ni,1)*W1_ann(icol,iraw)
         U_y = U_y + Gi_3D_acc(icol,nb1,ni,2)*W1_ann(icol,iraw)
         U_z = U_z + Gi_3D_acc(icol,nb1,ni,3)*W1_ann(icol,iraw)
        enddo
        U2f(1,iraw) = U_x*dfuN_acc(iraw,1,ni) 
        U2f(2,iraw) = U_y*dfuN_acc(iraw,1,ni) 
        U2f(3,iraw) = U_z*dfuN_acc(iraw,1,ni) 
       enddo                  ! exp(-u)/(1+exp(-u))^2 !

! -- Hidden Layers ---     
!$ACC LOOP SEQ
       do layer=2, net_layers-2
        NLcurr = Nodes_of_layer(layer)
        NLnext = Nodes_of_layer(layer+1)

        do i=1,NLcurr
         U1f(1,i)=U2f(1,i); U1f(2,i)=U2f(2,i); U1f(3,i)=U2f(3,i)
        enddo

        do iraw=1,NLnext  ! 1.. 20 !
         U_x = 0.0; U_y = 0.0; U_z = 0.0
!$ACC LOOP SEQ
         do icol=1,NLcurr  ! 1.. 20 !
          U_x = U_x + U1f(1,icol)*W2_ann(icol,iraw,layer)
          U_y = U_y + U1f(2,icol)*W2_ann(icol,iraw,layer)
          U_z = U_z + U1f(3,icol)*W2_ann(icol,iraw,layer)
         enddo
          U2f(1,iraw) = U_x*dfuN_acc(iraw,layer,ni)
          U2f(2,iraw) = U_y*dfuN_acc(iraw,layer,ni)
          U2f(3,iraw) = U_z*dfuN_acc(iraw,layer,ni)
        enddo               ! exp(-u)/(1+exp(-u))^2
       enddo ! do layer=2, net_layers-1 !

! -- Output Layer ---

       do iraw=1,Nodes_of_layer(net_layers)  ! 1.. 1 !
        U_x = 0.0; U_y = 0.0; U_z = 0.0
!$ACC LOOP SEQ
        do icol=1,Nodes_of_layer(net_layers-1)   ! 1.. 20 !
         U_x = U_x + U2f(1,icol)*W3_ann(icol,iraw)
         U_y = U_y + U2f(2,icol)*W3_ann(icol,iraw)
         U_z = U_z + U2f(3,icol)*W3_ann(icol,iraw)
        enddo
        dANN_dxij(1,iraw,nb1,ni) = U_x
        dANN_dxij(2,iraw,nb1,ni) = U_y
        dANN_dxij(3,iraw,nb1,ni) = U_z
       enddo

! --- End of ANN for each (i-j) pair using gij as input vectors ---

       enddo ! do nb1 = 1,maxn_nbrs !
      enddo ! do ni = 1, natoms !

      return
      end subroutine      ! dANN_ACC  !
c
c ---------------------------------------------------------------------
c Calculates Analytical derivatives and force calculation.
c ---------------------------------------------------------------------
c
      subroutine Frc_ANN_ACC(ecohe)

      use sys_ACC
      use atoms

      implicit double precision (a-h,o-z)
  
      double precision, intent(out) :: ecohe

C     write(6,*)'Frc_ANN_ACC: natoms=',natoms,' max_nbrs=',max_nbrs

      integer ialloc(16)
      ialloc(:) = 0

      max_ls = 1
      do i=2,net_layers-1
       if (Nodes_of_layer(i).gt.max_ls) max_ls=Nodes_of_layer(i)
      enddo

      allocate(r_ij_dev(max_nbrs,natoms), stat=ialloc(1))
      allocate(x_ij_dev(max_nbrs,natoms), stat=ialloc(2))
      allocate(y_ij_dev(max_nbrs,natoms), stat=ialloc(3))
      allocate(z_ij_dev(max_nbrs,natoms), stat=ialloc(4))
      allocate(fsij_acc(max_nbrs,n_set_ann,natoms), stat=ialloc(5))
      allocate(dfs_dxij(max_nbrs,natoms,n_set_ann), stat=ialloc(6))
      allocate(dfs_dyij(max_nbrs,natoms,n_set_ann), stat=ialloc(7))
      allocate(dfs_dzij(max_nbrs,natoms,n_set_ann), stat=ialloc(8))
      allocate(Gi_acc(mG_dim,natoms), stat=ialloc(9))
      allocate(Gi_3D_acc(mG_dim,max_nbrs,natoms,3), stat=ialloc(10))
      allocate(dfuN_acc(max_ls,net_layers-2,natoms), stat=ialloc(11))
      allocate(U1f(3,max_ls), stat=ialloc(12))
      allocate(U2f(3,max_ls), stat=ialloc(13))

      ierr = 0
      do i=1,16
       ierr = ierr + ialloc(i)
      enddo
      if(ierr.ne.0) then
       write(6,*)'ERROR deallocating x in Frc_ANN_ACC'
       ihalt = 1
       return
      endif

      h11 = h(1,1); h12 = h(1,2); h13 = h(1,3)
      h22 = h(2,2); h23 = h(2,3); h33 = h(3,3)

!$ACC DATA COPYIN(nbr_list,sx,sy,sz,r0Rc,r0pRc)
!$ACC& PRESENT(Nodes_of_layer,r0_value, base_pot_param)
!$ACC& PRESENT(W1_ann, W2_ann, W3_ann, B1_ann, B2_ann, B3_ann)
!$ACC& CREATE(r_ij_dev,x_ij_dev,y_ij_dev,z_ij_dev) 
!$ACC& CREATE(dfs_dxij, dfs_dyij, dfs_dzij,fsij_acc)
!$ACC& CREATE(Gi_acc, Gi_3D_acc, dfuN_acc)
!$ACC& CREATE(dANN_dxij)
!$ACC& COPYOUT(Ep_of, frr)

      call dANN_ACC(max_nbrs,Rc_ann)

      if(ihalt.ne.0) return

! --- Calc Actual Force Vectors ---

      ecohe = 0.0

!$ACC PARALLEL LOOP REDUCTION(+:ecohe)
      do ni = 1,natoms   ! loop over atoms !

       frr1 = 0.0; frr2 = 0.0; frr3 = 0.0

!$ACC LOOP REDUCTION(+:frr1,frr2,frr3)
       do nb1 = 1,max_nbrs     ! Loop I: i - j bond !
        nj = nbr_list(nb1,ni)
        ij_delta = min(abs(nj-ni),1) ! 0 if ni=nj; 1 if ni=/=nj !
        nbrs_of_j = ij_delta*max_nbrs  ! rij < Rc !

        nbi=nb1
!$ACC LOOP SEQ
        do nb2 = 1,nbrs_of_j  ! search for i as a neighbor of j !
         nj2 = nbr_list(nb2,nj)
         if(nj2.eq.ni) nbi=nb2       ! i is the nbi-th nbr of j !
        enddo
 
               ! Y3ij - Y3ji !
      fij1 = dANN_dxij(1,1,nb1,ni) - dANN_dxij(1,1,nbi,nj)
      fij2 = dANN_dxij(2,1,nb1,ni) - dANN_dxij(2,1,nbi,nj)
      fij3 = dANN_dxij(3,1,nb1,ni) - dANN_dxij(3,1,nbi,nj)
      frr1 = frr1 + fij1   ! x,y,z forces !
      frr2 = frr2 + fij2   ! x,y,z forces !
      frr3 = frr3 + fij3   ! x,y,z forces !

       enddo ! do nb1 = 1,max_nbrs !

       frr(1,ni) = frr1
       frr(2,ni) = frr2
       frr(3,ni) = frr3

       Ep = U3_of(1,ni)
       Ep_of(ni) = 2.0*Ep  ! Twice the atomic energy  !
       ecohe = ecohe + Ep   ! devided by 2 later in MSR for !
                            ! compatibility with other potentials !
      enddo ! do ni = 1,natoms !

!$ACC END DATA

      if (allocated(r_ij_dev)) deallocate(r_ij_dev, stat=ialloc(1))
      if (allocated(x_ij_dev)) deallocate(x_ij_dev, stat=ialloc(2))
      if (allocated(y_ij_dev)) deallocate(y_ij_dev, stat=ialloc(3))
      if (allocated(z_ij_dev)) deallocate(z_ij_dev, stat=ialloc(4))
      if (allocated(fsij_acc)) deallocate(fsij_acc, stat=ialloc(5))
      if (allocated(dfs_dxij))deallocate(dfs_dxij,stat=ialloc(6))
      if (allocated(dfs_dyij))deallocate(dfs_dyij,stat=ialloc(7))
      if (allocated(dfs_dzij))deallocate(dfs_dzij,stat=ialloc(8))
      if (allocated(Gi_acc)) deallocate(Gi_acc, stat=ialloc(9))
      if (allocated(Gi_3D_acc)) deallocate(Gi_3D_acc, stat=ialloc(10))
      if (allocated(dfuN_acc)) deallocate(dfuN_acc, stat=ialloc(11))
      if (allocated(U1f)) deallocate(U1f, stat=ialloc(12))
      if (allocated(U2f)) deallocate(U2f, stat=ialloc(13))

      ierr = 0
      do i=1,16
       ierr = ierr + ialloc(i)
      enddo
      if(ierr.ne.0) then
       write(6,*)'ERROR deallocating x in Frc_ANN_ACC'
       ihalt = 1
      endif

      return
      end subroutine      ! Frc_ANN_ACC  !
!
! ------------------------------------------------------------------
! ------------------------------------------------------------------
!
      END MODULE  ! ANN_ACC !
