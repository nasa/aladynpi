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
! Use a trained Physically Informed NN to get pot. param for a 
! specific potential form
!------------------------------------------------------------------
!
      MODULE PINN_ACC

      use sys_OMP
      use constants
      use sim_box
      use pot_module
      use ANN
      use ANN_ACC
      use PINN
      save

      double precision, dimension(:,:), allocatable :: FcSij_,dFcSij_, 
     1 Bb1Sij_,Pb1Sij_, Bb3Sij_,Pb3Sij_, Bb1Sji_,Pb1Sji_, PsidFij_dev

      double precision, dimension(:,:), allocatable :: dE7Sij_list_,
     1 dBZij_list_, PsidFij_list_, FSBb3C2_ij_list_, FSPb3C2_ij_list_,
     1 BZij_list_, BSij_list_, F2_rij_BOP_, F3_rij_BOP_
      double precision, dimension(:,:), allocatable :: FCxij_BOP_,
     1 FCyij_BOP_, FCzij_BOP_
      double precision, dimension(:,:,:), allocatable :: dE_dp

      CONTAINS
! ---------------------------------------------------------------------
c
c-----------------------------------------------------------------------
c Makes an invert list of neighbors nabors_num for each atom
c from nbr_list() to create nbr_list_inv()
c ACC version 
c-----------------------------------------------------------------------
c
      subroutine nbrs_BOP_ACC_inv

      use atoms

      implicit double precision (a-h,o-z)

!$OMP PARALLEL DO DEFAULT(NONE)
!$OMP& PRIVATE(l,ij_delta,nbrs_of_l,nbr_nr,n_nr,l2,nbrs2_of_l)
!$OMP& SHARED(natoms, max1_nbrs, max2_nbrs, nbr_list, nbr_list_inv)
!$OMP& SCHEDULE(STATIC)

       do nr=1,natoms

        do nbl = 1,max1_nbrs        ! loop over 1st nbrs of nr !
         l = nbr_list(nbl,nr)
         ij_delta = min(abs(l-nr),1) ! 0 if nr=l; 1 if nr=/=l !
         nbrs_of_l = ij_delta*max1_nbrs  ! rij < Rc !

c loop over 1st nbrs of l to find the nr position as a neighbor of l
         
         nbr_nr = nbl  ! assuming nr = l, so that nbr_inv = nbr !
         do n_nr = 1,nbrs_of_l
          l2 = nbr_list(n_nr,l)     ! l2 is a 1st neighbor of l !
          if(l2.eq.nr) nbr_nr=n_nr
         enddo
         nbr_list_inv(nbl,nr) = nbr_nr  ! or = nbl if nr=l !
        enddo  ! do nbl = 1,max1_nbrs !

        do nbl = max1_nbrs+1,max2_nbrs  ! loop over 2nd nbrs of nr !
         l = nbr_list(nbl,nr)
         ij_delta = min(abs(l-nr),1)    ! 0 if nr=l; 1 if nr=/=l !
         nbr_nr = nbl     ! assuming nr = l, so that nbr_inv = nbr !
         nbrs2_of_l = ij_delta*max2_nbrs  ! Rc < rij < 1.5Rc !
         do n_nr = max1_nbrs+1,nbrs2_of_l
          l2 = nbr_list(n_nr,l)        ! l2 is a 2nd neighbor of l !
          if(l2.eq.nr) nbr_nr=n_nr
         enddo
         nbr_list_inv(nbl,nr) = nbr_nr ! or = nbl if nr=l !

        enddo ! do nbl = 1,max2_nbrs !

                  ! nr = nbr_list(nbr_nr,l)
                  ! nr is nbr_list_inv(nbl,nr)-th nbr of l !
                  ! or
                  ! l = nbr_list(nbl,nr)
                  ! nr = nbr_list(nbr_list_inv(nbl,nr),l)
                  
       enddo ! do nr = 1,natoms !

!$OMP END PARALLEL DO

      return
      end subroutine ! nbrs_BOP_ACC_inv !
!
! ---------------------------------------------------------------------
! Calculates Analytical derivatives of the BOP function.
! ---------------------------------------------------------------------
!
      subroutine dBOP_ACC(ecohe)

      use atoms

      implicit double precision (a-h,o-z)

      double precision, intent(out) :: ecohe
 
       Rc = Rc_ann; Rc2 = Rc_ann**2; d4 = d4_ann 
       ecohe = 0.d0

c ------ START ENERGY BOP PART ------
c
!$ACC KERNELS LOOP WORKER COLLAPSE (2) VECTOR_LENGTH(32)
      do ni=1,natoms   ! II loop over atoms !
       do nb1 = 1, max2_nbrs     ! rij < Rc !
c
c II loop over neighbors of ni to get Screening factor, Sij(nb1)
c
c   o -------- o ........ o
c   j=nb1      i         k=j2=nb2
c
        nj = nbr_list(nb1,ni)
        rij = r_ij_dev(nb1,ni)  ! rij < Rc !
        xij = x_ij_dev(nb1,ni)*rij  ! x = x/r * r !
        yij = y_ij_dev(nb1,ni)*rij
        zij = z_ij_dev(nb1,ni)*rij

        Sij1 = 1.d0; dE7Sij = 0.d0

!$ACC LOOP VECTOR REDUCTION(+:dE7Sij) REDUCTION(*:Sij1)
        do nb2 = 1, max2_nbrs ! Screen factor: i-j bond; rik < 1.5Rc !

         k = nbr_list(nb2,ni)

         ! Exclude j=k !
         jk_delta = min(abs(nj-k),1) ! 0 if nj=k; 1 if nj=/=k !

         rik = r_ij_dev(nb2,ni)   ! rik < 1.5Rc !
         xik = x_ij_dev(nb2,ni)*rik
         yik = y_ij_dev(nb2,ni)*rik
         zik = z_ij_dev(nb2,ni)*rik

         xjk = xik - xij
         yjk = yik - yij
         zjk = zik - zij
         rjk = sqrt(xjk**2 + yjk**2 + zjk**2)

         rrr = rik + rjk - rij
         bl = U3_of(nlambda_bop,ni)
         rrc = jk_delta*max(Rc-rrr,0.d0)
         rrc4 = rrc**4
         Fc = rrc4 / (d4+rrc4)
         Exp_blr = exp(-bl*rrr)
         Sijk_term = 1.d0 - Fc*Exp_blr
         Sij1 = Sij1 * Sijk_term

         dE7Sijk = Fc*Exp_blr*rrr / Sijk_term
         dE7Sij = dE7Sij + dE7Sijk

        enddo ! do nb2 = 1, max2_nbrs  ! Screening: j-k bond !

        RcRij = max(Rc-rij, 0.d0) ! rij < Rc case !
        RcRij4 = RcRij**4
        fc_rij = RcRij4 / (d4+RcRij4)

        BSij_list_(nb1,ni) = Sij1
        dE7Sij_list_(nb1,ni) = fc_rij*Sij1*dE7Sij

       enddo ! do nb1 = 1,max2_nbrs   ! Screening: i - j bond !
      enddo ! do ni=1,natoms   ! I loop over atoms !

c --- Start another loop over neighbors of ni: Chem i - j bond ---

!$ACC PARALLEL LOOP COLLAPSE (2) VECTOR_LENGTH(32)
      do ni=1,natoms   ! II loop over atoms !
       do nb1 = 1,max2_nbrs  ! rij < Rc !

        nj = nbr_list(nb1,ni)

        xij = x_ij_dev(nb1,ni)  ! x/r; = 0 if j=i !
        yij = y_ij_dev(nb1,ni)  ! y/r !
        zij = z_ij_dev(nb1,ni)  ! z/r !

        Zij1 = 0.d0; dZij1 = 0.d0

!$ACC LOOP REDUCTION(+:Zij1,dZij1)
        do nb2 = 1,max2_nbrs   ! Three body: i-k bond; rik < Rc !
         k = nbr_list(nb2,ni)

         ! Exclude j=k !
         jk_delta = min(abs(nj-k),1) ! 0 if nj=k; 1 if nj=/=k !

         rik = r_ij_dev(nb2,ni)  ! rik < Rc !
         xik = x_ij_dev(nb2,ni)  ! = 0 if k=i !
         yik = y_ij_dev(nb2,ni)
         zik = z_ij_dev(nb2,ni)

         cos_ijk = xij*xik + yij*yik + zij*zik
         cosTh_h2 = (cos_ijk - U3_of(nh_bop,ni))**2

         RcRik = max(Rc-rik, 0.d0) ! rik < Rc case !
         RcRik4 = RcRik**4
         Fc_ik = jk_delta*RcRik4 / (d4+RcRik4)  ! 0 if j=k !

         S_ik = BSij_list_(nb2,ni)

         FcSikCh2 = Fc_ik*S_ik*cosTh_h2
         Zij1 = Zij1 + U3_of(naa_bop,ni)*FcSikCh2
         dZij1 = dZij1 + FcSikCh2
        enddo ! do nb2 = 1,max2_nbrs   ! Three body: i-k bond !

        BZij_list_(nb1,ni) = Zij1
        dBZij_list_(nb1,ni) = dZij1

       enddo ! do nb1 = 1,max2_nbrs   ! Screening: i - j bond !
      enddo ! do ni=1,natoms   ! I loop over atoms !

c --- Start Final loop over neighbors of ni: 2-body terms ---

!$ACC PARALLEL LOOP VECTOR_LENGTH(32) REDUCTION(+:ecohe)
      do ni=1,natoms   ! II loop over atoms !

       Ei = 0.d0; Psi = 0.d0

!$ACC LOOP REDUCTION(+:Ei,Psi)
       do nb1 = 1,max2_nbrs  ! rij < Rc !
        rij = r_ij_dev(nb1,ni)
        RcRij = max(Rc-rij, 0.d0) ! rij < Rc case !
        RcRij3 = RcRij**3
        RcRij4 = RcRij**4
        Fc_ij = RcRij4 / (d4+RcRij4)
        dFc_ij = -4.d0*RcRij3*d4 / (d4+RcRij4)**2    ! dFc/dr !

        S_ij = BSij_list_(nb1,ni)
        Z_ij = BZij_list_(nb1,ni)

        b_ij = 1.d0/sqrt(1.d0 + Z_ij)
        Psi_ij = S_ij*b_ij  ! Sij*bij !

        Psi = Psi + Fc_ij*Psi_ij       ! e_i = \sum_j(Fc*Sij*bij) !

        VA_ij =  ! Aij*exp(-alpha*rij) !
     =  exp(U3_of(nA_bop,ni)-U3_of(nalpha_bop,ni)*rij)

        VB_ij_exp =  ! Bij*exp(-beta.rij) !
     =  exp(U3_of(nB_bop,ni)-U3_of(nbeta_bop,ni)*rij)
        VB_ij = VB_ij_exp * Psi_ij ! Sij*bij*Bij*exp(-beta.rij) !

        U_ij = VA_ij - VB_ij
        Ei = Ei + U_ij*Fc_ij

        dVA_ij = -U3_of(nalpha_bop,ni)*VA_ij
        dVB_ij = -U3_of(nbeta_bop,ni)*VB_ij
        dU_ij = dVA_ij - dVB_ij
        ! 2 body dU_ij at fixed Sij, bij(zij) !

        F2_rij_BOP_(nb1,ni) = (dU_ij*Fc_ij + U_ij*dFc_ij)

        PsidFij_list_(nb1,ni) = Psi_ij*dFc_ij ! Sij*bij*dFc/dr !
       enddo ! do nb1 = 1,max2_nbrs   ! Screening: i - j bond !

       Epot = 0.50*Ei - U3_of(nsigma_bop,ni)*sqrt(Psi)

       ecohe = ecohe + Epot

       Ep_of(ni) = 2.d0*Epot  ! Twice Epot !
       Psi_of(ni) = max(Psi,1.d-8)

      enddo ! do ni=1,natoms !
c
c ------ START FORCE BOP PART ------
c
c ------------- Loop I: --------------
c
!$ACC PARALLEL LOOP COLLAPSE (2) VECTOR_LENGTH(32)
      do ni=1,natoms
       do nb1 = 1,max2_nbrs     ! Loop 0: i - j bond; rij < 1.5Rc !

        nj = nbr_list(nb1,ni)  ! VECTORIZED: speedup: 1.550      !
        rij = r_ij_dev(nb1,ni)

        RcRij = max(Rc-rij, 0.d0) ! rij < Rc case !
        RcRij3 = RcRij**3
        RcRij4 = RcRij**4
        Fcij = RcRij4 / (d4+RcRij4)
        dFcij = -4.d0*RcRij3*d4 / (d4+RcRij4)**2    ! dFc/dr !

        i_at_j = nbr_list_inv(nb1,ni)   ! get i as a nbr of j !

        S_ij = BSij_list_(nb1,ni)        ! Sij !
        S_ji = BSij_list_(i_at_j,nj)     ! Sji !
        Z_ij = BZij_list_(nb1,ni)        ! Zij !
        Z_ji = BZij_list_(i_at_j,nj)     ! Zji !
        b_ij = 1.d0/sqrt(1.d0 + Z_ij)   ! bij !
        b_ji = 1.d0/sqrt(1.d0 + Z_ji)   ! bji !

        FcSij = Fcij*S_ij                        ! Fc(rij)*Sij  !
        FcSij_(nb1,ni) = FcSij
        dFcSij_(nb1,ni) = dFcij*S_ij             ! dFc(rij)*Sij !
        Fb1Sij = FcSij*b_ij                      ! Fc*Sij*bij   !
        Fb1Sji = Fcij*S_ji*b_ji                  ! Fc*Sji*bji   !
        Fb3Sij = FcSij*b_ij**3                   ! Fc*Sij*bij^3 !
        dE7Sij = dE7Sij_list_(nb1,ni)

        Sigma_Psi_i = U3_of(nsigma_bop,ni)/sqrt(Psi_of(ni))
                  ! = sigma_i*(Psi(i)^-1/2) !
        Sigma_Psi_j = U3_of(nsigma_bop,nj)/sqrt(Psi_of(nj))
                  ! = sigma_j*(Psi(j)^-1/2) !

        VA_ij_exp =  ! Aij*exp(-alpha*rij) !
     =  exp(U3_of(nA_bop,ni)-U3_of(nalpha_bop,ni)*rij)
        VB_ij_exp =  ! Bij*exp(-beta.rij) !
     =  exp(U3_of(nB_bop,ni)-U3_of(nbeta_bop,ni)*rij)

        Bb1ij = VB_ij_exp*b_ij    ! Bij*exp(-beta.rij)*bij !
        Pb1ij = Sigma_Psi_i*b_ij  ! sigma_i*(Psi(i)^-1/2)*bij     !
        Bb1Sij = VB_ij_exp*Fb1Sij ! Bij*exp(-beta.rij)*Fc*Sij*bij !
        Bb1Sij_(nb1,ni) = Bb1Sij
        Pb1Sij_(nb1,ni) = Sigma_Psi_i*Fb1Sij 
                     ! sigma_i*(Psi(i)^-1/2)*Fc*Sij*bij !

        VB_ji_exp =  ! Bji*exp(-beta.rij) !
     =  exp(U3_of(nB_bop,nj)-U3_of(nbeta_bop,nj)*rij)

        Bb1Sji_(nb1,ni) = VB_ji_exp*Fb1Sji
                     ! Bji*exp(-beta.rij)*Fc*Sji*bji !
        Pb1Sji_(nb1,ni) = Sigma_Psi_j*Fb1Sji 
                     ! sigma_j*(Psi(j)^-1/2)*Fc*Sji*bji !

        Bb3Sij = VB_ij_exp*Fb3Sij   
        Bb3Sij_(nb1,ni) = Bb3Sij  
                     ! Bij*exp(-beta.rij)*Fc*Sij*bij^3
        Pb3Sij = Sigma_Psi_i*Fb3Sij 
        Pb3Sij_(nb1,ni) =  Pb3Sij
                     ! sigma_i*(Psi(i)^-1/2)*Fc*Sij*bij^3

        VA_ij_Fc = VA_ij_exp*Fcij   ! Fc*Aij*exp(-alpha*rij) !

        dE_dp(:,nb1,ni) = 0.d0

        dE_dp(1,nb1,ni) =  VA_ij_Fc 
        dE_dp(2,nb1,ni) = -VA_ij_Fc*rij 
        dE_dp(3,nb1,ni) = -Bb1Sij
        dE_dp(4,nb1,ni) =  Bb1Sij*rij
        dE_dp(naa_bop,nb1,ni) = 
     =  0.50*(Bb3Sij+Pb3Sij)*dBZij_list_(nb1,ni)
        dE_dp(nlambda_bop,nb1,ni) = -(Bb1ij+Pb1ij)*dE7Sij

       enddo ! do nb1 = 1,max2_nbrs !
      enddo ! do ni=1,natoms !

!$ACC PARALLEL LOOP COLLAPSE (2) VECTOR_LENGTH(32) ! Fastest !
      do ni=1,natoms
       do nb1 = 1,max2_nbrs     ! Loop I: i - j bond; rij < 1.5Rc !
        nj = nbr_list(nb1,ni)

        rij = r_ij_dev(nb1,ni)
        xr_ij = x_ij_dev(nb1,ni) ! x/r !
        yr_ij = y_ij_dev(nb1,ni) ! y/r !
        zr_ij = z_ij_dev(nb1,ni) ! z/r !
        rij1 = 1.d0/rij

        xij = xr_ij*rij   ! x !
        yij = yr_ij*rij   ! x !
        zij = zr_ij*rij   ! x !
c
c --- Calc. dSij_dxi and dSji_dxi ---
c --- Loop i - k bond to calc the 1/Sijk and dSij_dxi terms ---
c
        FcSij = FcSij_(nb1,ni)

        Bb1Sij = Bb1Sij_(nb1,ni)   ! Bij*exp(-beta.rij)*Fc*Sij*bij
        Pb1Sij = Pb1Sij_(nb1,ni)   ! sigma_i*(Psi(i)^-1/2)*Fc*Sij*bij

        Bb3Sij = Bb3Sij_(nb1,ni)   ! Bij*exp(-beta.rij)*Fc*Sij*bij^3
        Pb3Sij = Pb3Sij_(nb1,ni)   ! sigma_i*(Psi(i)^-1/2)*Fc*Sij*bij^3

        dE7Sij = dE7Sij_list_(nb1,ni)

        dSr_ij = 0.d0; dEsr_ij = 0.d0; dWsr_ij = 0.d0
        Bb3C2_ij = 0.d0; dWbfr_ij = 0.d0

        dE1bC1 = 0.d0; dE2bC1 = 0.d0; dW2bC1 = 0.d0
        dE1bC2 = 0.d0; dE2bC2 = 0.d0; dW2bC2 = 0.d0
        dE1bC3 = 0.d0; dE2bC3 = 0.d0; dW2bC3 = 0.d0
        dE6 = 0.d0; dW6 = 0.d0
        dE7 = 0.d0; dW7 = 0.d0

! This is faster for GTX1080, but slower for v100 !
! $ACC LOOP SEQ  ! Too many reductions to parallelize !

! This is faster for V100, but slower for GTX1080 !
!$ACC LOOP
!$ACC& REDUCTION(+:dSr_ij,dEsr_ij,dWsr_ij,Bb3C2_ij,dWbfr_ij)
!$ACC& REDUCTION(+:dE1bC1,dE1bC2,dE1bC3, dE2bC1,dE2bC2,dE2bC3)
!$ACC& REDUCTION(+:dW2bC1,dW2bC2,dW2bC3, dE6,dE7, dW6,dW7)
        do nbk = 1,max2_nbrs     ! rij < Rc; rik < Rc !
         k = nbr_list(nbk,ni)

         ! Exclude j=k !
         njk_delta = min(abs(nj-k),1) ! 0 if nj=k; 1 if nj=/=k !

         FcSik = FcSij_(nbk,ni)
         dE7Sik = dE7Sij_list_(nbk,ni)

         rik = r_ij_dev(nbk,ni)
         xik = x_ij_dev(nbk,ni)*rik
         yik = y_ij_dev(nbk,ni)*rik
         zik = z_ij_dev(nbk,ni)*rik
         rik1= 1.d0/rik

         xjk = xik - xij
         yjk = yik - yij
         zjk = zik - zij
         r2jk = xjk**2 + yjk**2 + zjk**2

         rjk = sqrt(r2jk)

         xr_ik = xik*rik1
         yr_ik = yik*rik1
         zr_ik = zik*rik1

         Bb1Sik = Bb1Sij_(nbk,ni) ! Bik*exp(-beta.rik)*Fc*Sik*bik !
         Pb1Sik = Pb1Sij_(nbk,ni) ! sigma_i*(Psi(i)^-1/2)*Fc*Sik*bik !

         Bb1Ski = Bb1Sji_(nbk,ni) ! Bki*exp(-beta.rik)*Fc*Ski*bki ! 
         Pb1Ski = Pb1Sji_(nbk,ni) ! sigma_k*(Psi(k)^-1/2)*Fc*Ski*bki ! 

         rrr = rij + rjk - rik
         bl_ikj = U3_of(nlambda_bop,ni)
         bl_kij = U3_of(nlambda_bop,k)
         rrc = njk_delta * max(Rc-rrr,0.d0)  ! 0 if j=k OR rrr > Rc !
         rrc3 = rrc**3
         rrc4 = rrc*rrc3
         Fc_rrr = rrc4 / (d4+rrc4)          ! Fc(r_ikj) = Fc(r_kij) !
         Exp_blr_ikj = exp(-bl_ikj*rrr)
         Exp_blr_kij = exp(-bl_kij*rrr)
         Sikj_term = 1.d0 - Fc_rrr*Exp_blr_ikj  ! S_ikj !
         Skij_term = 1.d0 - Fc_rrr*Exp_blr_kij  ! S_kij !

         dFc = -4.d0*rrc3*d4 / (d4+rrc4)**2     ! dFc/dr !
         dS_ikj = (bl_ikj*Fc_rrr - dFc)*Exp_blr_ikj
         dS_kij = (bl_kij*Fc_rrr - dFc)*Exp_blr_kij
         dS_ikj = dS_ikj/Sikj_term
         dS_kij = dS_kij/Skij_term

         Bb1dSr_ikj = Bb1Sik*dS_ikj
         ! Bik*exp(-beta.rik)*Fc*Sik*bik*dS_ikj/Sikj_term !
         Bb1dSr_kij = Bb1Ski*dS_kij
         ! Bki*exp(-beta.rik)*Fc*Ski*bki*dS_kij/Skij_term !
         Pb1dSr_ikj = Pb1Sik*dS_ikj
         ! sigma_i*(Psi(i)^-1/2)*Fc*Sik*bik*dS_ikj/Sikj_term !
         Pb1dSr_kij = Pb1Ski*dS_kij
         ! sigma_k*(Psi(k)^-1/2)*Fc*Ski*bki*dS_kij/Skij_term !

         rrr = rik + rjk - rij
         bl_ijk = U3_of(nlambda_bop,ni)
         rrc = njk_delta * max(Rc-rrr,0.d0)  ! 0 if j=k OR rrr > Rc !
         rrc3 = rrc**3
         rrc4 = rrc*rrc3
         Fc_rrr = rrc4 / (d4+rrc4)          ! Fc(r_ijk) = Fc(r_jik) !
         Exp_blr_ijk = exp(-bl_ijk*rrr)
         Sijk_term = 1.d0 - Fc_rrr*Exp_blr_ijk  ! S_ijk !

         dFc = -4.d0*rrc3*d4 / (d4+rrc4)**2     ! dFc/dr !
         dS_ijk = (bl_ijk*Fc_rrr - dFc)*Exp_blr_ijk
         dS_ijk = dS_ijk/Sijk_term

         Bb1dSr_ijk = Bb1Sij*dS_ijk
         ! Bij*exp(-beta.rij)*Fc*Sij*bij*dS_ijk/Sijk_term !
         Pb1dSr_ijk = Pb1Sij*dS_ijk
         ! sigma_i*(Psi(i)^-1/2)*Fc*Sij*bij*dS_ijk/Sijk_term !

         dSr_ij = dSr_ij + dS_ijk

         dEsr_ij = dEsr_ij + Bb1dSr_ikj + Bb1dSr_kij - Bb1dSr_ijk
         dWsr_ij = dWsr_ij + Pb1dSr_ikj + Pb1dSr_kij - Pb1dSr_ijk

         aa_ikj = njk_delta*U3_of(naa_bop,ni)
         h_ikj = U3_of(nh_bop,ni)
         aa_ijk = njk_delta*U3_of(naa_bop,ni)
         h_ijk = U3_of(nh_bop,ni)

       cos_ikj = xr_ij*xr_ik + yr_ij*yr_ik + zr_ij*zr_ik
       cos_h_ikj = cos_ikj - h_ikj
       cos_h_ijk = cos_ikj - h_ijk
       cos2_h_ikj = cos_h_ikj**2
       cos2_h_ijk = cos_h_ijk**2

         xr_ik_ij = (xr_ik - xr_ij*cos_ikj)*rij1
         yr_ik_ij = (yr_ik - yr_ij*cos_ikj)*rij1
         zr_ik_ij = (zr_ik - zr_ij*cos_ikj)*rij1
                     ! (xik/rik - xij/rij*cos_ikj)/rij !
         C1ikj = aa_ikj*cos_h_ikj
         C1ijk = aa_ijk*cos_h_ijk
         C2ikj = aa_ikj*cos2_h_ikj
         C2ijk = aa_ijk*cos2_h_ijk

         Bb3Sik =  Bb3Sij_(nbk,ni)
         BexpFb3SijC1 =  Bb3Sij*C1ijk
       ! Bij*exp(-beta.rij)*Fc*Sij*(bij**3)*aa_ijk*cos_h_ijk !
         BexpFb3SikC1 =  Bb3Sik*C1ikj
       ! Bik*exp(-beta.rik)*Fc*Sik*(bik**3)*aa_ikj*cos_h_ikj !
         BexpFb3SikC2 =  Bb3Sik*C2ikj
       ! Bik*exp(-beta.rik)*Fc*Sik*(bik**3)*a_ikj*cos2_h_ikj !
         BexpFb3SijC2 =  Bb3Sij*C2ijk
       ! Bij*exp(-beta.rij)*Fc*Sij*(bij**3)*a_ijk*cos2_h_ijk !

         Pb3Sik = Pb3Sij_(nbk,ni)
         PiFb3SijC1 = Pb3Sij*C1ijk
       ! sigma_i*(Psi(i)^-1/2)*Fc*Sij*(bij**3)*aa_ijk*cos_h_ijk    !
         PiFb3SikC1 = Pb3Sik*C1ikj
       ! sigma_i*(Psi(i)^-1/2)*Fc*Sik*(bik**3)*aa_ikj*cos_h_ikj    !
         PiFb3SikC2 = Pb3Sik*C2ikj
       ! sigma_i*(Psi(i)^-1/2)*Fc*Sik*bik**3*a_ikj*cos2_h_ikj !
         PiFb3SijC2 = Pb3Sij*C2ijk
       ! sigma_i*(Psi(i)^-1/2)*Fc*Sij*bij**3*a_ijk*cos2_h_ijk !

         Bb3C2_ij = Bb3C2_ij + BexpFb3SikC2
         dWbfr_ij = dWbfr_ij + PiFb3SikC2

         dE1bC1 = dE1bC1 + FcSik*C1ijk*xr_ik_ij
         dE1bC2 = dE1bC2 + FcSik*C1ijk*yr_ik_ij
         dE1bC3 = dE1bC3 + FcSik*C1ijk*zr_ik_ij
       ! += Fc*Sik*aa_ijk*cos_h_ijk*(xik/rik - xij/rij*cos_ikj)/rij

         dE2bC1 = dE2bC1 + BexpFb3SikC1*xr_ik_ij
         dE2bC2 = dE2bC2 + BexpFb3SikC1*yr_ik_ij
         dE2bC3 = dE2bC3 + BexpFb3SikC1*zr_ik_ij
       ! += Bik*exp(-beta.rik)*Fc*Sik*(bik**3)*aa_ikj*cos_h_ikj*
       !   *(xik/rik - xij/rij*cos_ikj)/rij

         dW2bC1 = dW2bC1 + PiFb3SikC1*xr_ik_ij
         dW2bC2 = dW2bC2 + PiFb3SikC1*yr_ik_ij
         dW2bC3 = dW2bC3 + PiFb3SikC1*zr_ik_ij
       ! += sigma_i*(Psi(i)^-1/2)*Fc*Sik*bik*aa_ijk*cos_h_ijk*
       !   *(xik/rik - xij/rij*cos_ikj)/rij

         dE6 = dE6 + BexpFb3SijC1*FcSik 
       ! Bij*exp(-beta.rij)*FcSij*(bij**3)*aa_ijk*cos_h_ijk*FcSik !

         dW6 = dW6 + PiFb3SijC1*FcSik 
       ! sigma_i*(Psi(i)^-1/2)*FcSij*(bij**3)*aa_ijk*cos_h_ijk*FcSik !

         dE7 = dE7 + BexpFb3SijC2*dE7Sik
         dW7 = dW7 + PiFb3SijC2*dE7Sik

        enddo ! do nbk = 1,max2_nbrs  ! Loop B: i - j - j2 bond !

        dE_dp(nh_bop,nb1,ni) = -(dE6 + dW6)
        dE_dp(nlambda_bop,nb1,ni) = 
     =  dE_dp(nlambda_bop,nb1,ni) + 0.50*(dE7 + dW7)

        Bb3Sij =  Bb3Sij_(nb1,ni)  ! Bij*exp(-beta.rij)*Fc*Sij*bij^3 !
        Pb3Sij =  Pb3Sij_(nb1,ni)  ! sigma_i*(Psi(i)^-1/2)*Fc*Sij*bij^3
        FSBb3C2ij = FcSij*Bb3C2_ij ! Fc*Sij*Bb3C2_ij !
        FSPb3C2ij = FcSij*dWbfr_ij ! Fc*Sij*Pb3C2_ij !

        FSBb3C2_ij_list_(nb1,ni) = FSBb3C2ij
! Fij*Sij*Sum_k{Bik*exp(-beta.rik)*Fc*Sik*bik**3*a_ikj*cos2_h_ikj} !
        FSPb3C2_ij_list_(nb1,ni) = FSPb3C2ij
! Fij*Sij*Sum_k{sigma_i*(Psi(i)^-1/2)*Fc*Sik*bik**3*a_ikj*cos2_h_ikj} !

        dEbS3r_ij = FSBb3C2ij*dSr_ij  ! 3rd term of dEbSr_ij !
        dWbS3r_ij = FSPb3C2ij*dSr_ij  ! 3rd term of dWbSr_ij !
        dEWbfr_ij = dFcSij_(nb1,ni)*(Bb3C2_ij + dWbfr_ij)
        F3_rij_BOP_(nb1,ni) = dEsr_ij + dWsr_ij - 
     -  0.50*(dEWbfr_ij - dEbS3r_ij - dWbS3r_ij) 

        FCxij_BOP_(nb1,ni) = -((Bb3Sij + Pb3Sij)*dE1bC1 + 
     +                   FcSij*(dE2bC1 + dW2bC1))
        FCyij_BOP_(nb1,ni) = -((Bb3Sij + Pb3Sij)*dE1bC2 + 
     +                   FcSij*(dE2bC2 + dW2bC2))
        FCzij_BOP_(nb1,ni) = -((Bb3Sij + Pb3Sij)*dE1bC3 + 
     +                   FcSij*(dE2bC3 + dW2bC3))

       enddo ! do nb1 = 1,max2_nbrs !
      enddo ! do ni=1,natoms !
c
c ------------- End Loop I: ---------------------
c
c ------------- Collect dE_dparam ---------------------
c
!$ACC PARALLEL LOOP COLLAPSE (2)
      do ni=1,natoms
       do i=1,nBOP_params
        dE_dp_sum = 0.d0
!$ACC LOOP SEQ
        do nb1 = 1,max1_nbrs
         dE_dp_sum = dE_dp_sum + dE_dp(i,nb1,ni)
        enddo
        dE_dparam_list(i,ni) = 0.50*dE_dp_sum
       enddo
      enddo
 
!$ACC PARALLEL LOOP
      do ni=1,natoms           ! dE/dSigma = -Psi(i)^1/2 !
       dE_dparam_list(nsigma_bop,ni) = -sqrt(Psi_of(ni))
      enddo
c
c ------------- Loop II: --------------
c
!$ACC PARALLEL LOOP COLLAPSE (2) VECTOR_LENGTH(32)
      do ni=1,natoms
       do nb1 = 1,max2_nbrs     ! 0 < rij < 1.5Rc !
        nj = nbr_list(nb1,ni)

        rij = r_ij_dev(nb1,ni)
        xij = x_ij_dev(nb1,ni)*rij  ! x !
        yij = y_ij_dev(nb1,ni)*rij
        zij = z_ij_dev(nb1,ni)*rij

        dEbS1r_ij = 0.d0; dEbS2r_ij = 0.d0
        dWbS1r_ij = 0.d0; dWbS2r_ij = 0.d0

!$ACC LOOP REDUCTION(+:dEbS1r_ij,dWbS1r_ij,dEbS2r_ij,dWbS2r_ij)
        do nbk = 1,max2_nbrs     ! rij < 1.5Rc; rik < Rc !
         k = nbr_list(nbk,ni)
         ! Exclude j=k !
         jk_delta = min(abs(nj-k),1) ! 0 if nj=k; 1 if nj=/=k !

         rik = r_ij_dev(nbk,ni)
         xik = x_ij_dev(nbk,ni)*rik
         yik = y_ij_dev(nbk,ni)*rik
         zik = z_ij_dev(nbk,ni)*rik

         xjk = xik - xij
         yjk = yik - yij
         zjk = zik - zij
         rjk = sqrt(xjk**2 + yjk**2 + zjk**2)

          ! Calc. dSr_kij !

          rrr = rjk + rij - rik
          bl_kij = U3_of(nlambda_bop,k)
          bl_ikj = U3_of(nlambda_bop,ni)
          rrc = jk_delta * max(Rc-rrr,0.d0)
          rrc3 = rrc**3
          rrc4 = rrc*rrc3
          Fc_rrr = rrc4 / (d4+rrc4)              ! Fc(r_kij) !
          Exp_blr_kij = exp(-bl_kij*rrr)
          Exp_blr_ikj = exp(-bl_ikj*rrr)
          Skij_term = 1.d0 - Fc_rrr*Exp_blr_kij  ! S_kij !
          Sikj_term = 1.d0 - Fc_rrr*Exp_blr_ikj  ! S_ikj !

          dFc = -4.d0*rrc3*d4 / (d4+rrc4)**2     ! dFc/dr !
          dS_kij = (bl_kij*Fc_rrr - dFc)*Exp_blr_kij/Skij_term
          dS_ikj = (bl_ikj*Fc_rrr - dFc)*Exp_blr_ikj/Sikj_term

          i_at_k = nbr_list_inv(nbk,ni)  ! get i as a nbr of k !
          dEbS1r_ij = dEbS1r_ij + dS_kij*FSBb3C2_ij_list_(i_at_k,k)
                              ! + dS_kij*Fc_ik*Ski*Bb3C2_ki !
          dWbS1r_ij = dWbS1r_ij + dS_kij*FSPb3C2_ij_list_(i_at_k,k)
                              ! + dS_kij*Fc_ik*Ski*Pb3C2_ki !
         
          FSBb3C2_ik = FSBb3C2_ij_list_(nbk,ni)
          dEbS2r_ij = dEbS2r_ij + dS_ikj*FSBb3C2_ik
                              ! + dS_ikj*Fc_ik*Sik*Bb3C2_ik !
          FSPb3C2_ik = FSPb3C2_ij_list_(nbk,ni)
          dWbS2r_ij = dWbS2r_ij + dS_ikj*FSPb3C2_ik
                              ! + dS_ikj*Fc_ik*Sik*Pb3C2_ik !
        enddo ! do nbk = 1,max2_nbrs     ! rij < 1.5Rc; rik < Rc !

        ! Add the 1st term of dEbSr_ij and dWbSr_ij !
        F3_rij_BOP_(nb1,ni) = F3_rij_BOP_(nb1,ni) - 
     -  0.50*(dEbS1r_ij + dEbS2r_ij + dWbS1r_ij + dWbS2r_ij)

       enddo ! do nb1 = 1,max2_nbrs !
      enddo ! do ni=1,natoms !

      return
      end subroutine      ! dBOP_ACC !
c
c ---------------------------------------------------------------------
c Calculates Analytical derivatives and force calculation.
c ---------------------------------------------------------------------
c
      subroutine Frc_PINN_ACC(ecohe)

      use sys_ACC
      use atoms

      implicit double precision (a-h,o-z)
  
      double precision, intent(out) :: ecohe

      double precision, dimension(3) :: xr_ik_ij, fr

      integer ialloc(36)
      ialloc(:) = 0

C     write(6,10)natoms,max1_nbrs,max2_nbrs
C 10  format('ACC: Frc_PINN_ACC: natoms=',i8,' max1_nbrs=',i4,
C    & ' max2_nbrs=',i4)

      max_ls = 1
      do i=2,net_layers-1
       if (Nodes_of_layer(i).gt.max_ls) max_ls=Nodes_of_layer(i)
      enddo
 
      ecohe = 0.d0
  
      allocate(r_ij_dev(max2_nbrs,natoms), stat=ialloc(1))
      allocate(x_ij_dev(max2_nbrs,natoms), stat=ialloc(2))
      allocate(y_ij_dev(max2_nbrs,natoms), stat=ialloc(3))
      allocate(z_ij_dev(max2_nbrs,natoms), stat=ialloc(4))
      allocate(fsij_acc(max2_nbrs,n_set_ann,natoms), stat=ialloc(5))
      allocate(dfs_dxij(max2_nbrs,natoms,n_set_ann), stat=ialloc(6))
      allocate(dfs_dyij(max2_nbrs,natoms,n_set_ann), stat=ialloc(7))
      allocate(dfs_dzij(max2_nbrs,natoms,n_set_ann), stat=ialloc(8))
      allocate(Gi_acc(mG_dim,natoms), stat=ialloc(9))
      allocate(Gi_3D_acc(mG_dim,max2_nbrs,natoms,3), stat=ialloc(10))
      allocate(dfuN_acc(max_ls,net_layers-2,natoms), stat=ialloc(11))
      allocate(U1f(3,max_ls), stat=ialloc(12))
      allocate(U2f(3,max_ls), stat=ialloc(13))

      allocate(FcSij_(max2_nbrs,natoms), stat=ialloc(14))
      allocate(dFcSij_(max2_nbrs,natoms), stat=ialloc(15))
      allocate(Bb1Sij_(max2_nbrs,natoms), stat=ialloc(16))
      allocate(Pb1Sij_(max2_nbrs,natoms), stat=ialloc(17))
      allocate(Bb1Sji_(max2_nbrs,natoms), stat=ialloc(18))
      allocate(Pb1Sji_(max2_nbrs,natoms), stat=ialloc(19))
      allocate(Bb3Sij_(max2_nbrs,natoms), stat=ialloc(20))
      allocate(Pb3Sij_(max2_nbrs,natoms), stat=ialloc(21))
      allocate(PsidFij_dev(max2_nbrs,natoms), stat=ialloc(22))
      allocate(dE_dp(nBOP_params,max2_nbrs,natoms), stat=ialloc(23))

      allocate(dE7Sij_list_(max2_nbrs,natoms),stat=ialloc(24))
      allocate(dBZij_list_(max2_nbrs,natoms),stat=ialloc(25))
      allocate(PsidFij_list_(max2_nbrs,natoms),stat=ialloc(26))
      allocate(FSBb3C2_ij_list_(max2_nbrs,natoms),stat=ialloc(27))
      allocate(FSPb3C2_ij_list_(max2_nbrs,natoms),stat=ialloc(28))
      allocate(BZij_list_(max2_nbrs,natoms),stat=ialloc(29))
      allocate(BSij_list_(max2_nbrs,natoms),stat=ialloc(30))
      allocate(F2_rij_BOP_(max2_nbrs,natoms),stat=ialloc(31))
      allocate(F3_rij_BOP_(max2_nbrs,natoms),stat=ialloc(32))
      allocate(FCxij_BOP_(max2_nbrs,natoms),stat=ialloc(33))
      allocate(FCyij_BOP_(max2_nbrs,natoms),stat=ialloc(34))
      allocate(FCzij_BOP_(max2_nbrs,natoms),stat=ialloc(35))

      ierr = 0
      do i=1,36
       ierr = ierr + ialloc(i)
      enddo
      if(ierr.ne.0) then
       write(6,*)'ERROR deallocating x in Frc_PINN_ACC'
       ihalt = 1
       return
      endif

      h11 = h(1,1); h12 = h(1,2); h13 = h(1,3)
      h22 = h(2,2); h23 = h(2,3); h33 = h(3,3)

      Rc = Rc_ann; Rc2 = Rc_ann**2; d4 = d4_ann
      Rc15 = 1.50*Rc
      Rc15_sqr = Rc15**2

!$ACC DATA COPYIN(nbr_list,nbr_list_inv,sx,sy,sz,r0Rc,r0pRc)
!$ACC& PRESENT(Nodes_of_layer,r0_value, base_pot_param)
!$ACC& PRESENT(W1_ann, W2_ann, W3_ann, B1_ann, B2_ann, B3_ann)
!$ACC& COPYIN(nbrs_2nd_of, FcSij_, dFcSij_)
!$ACC& CREATE(r_ij_dev,x_ij_dev,y_ij_dev,z_ij_dev) 
!$ACC& CREATE(dfs_dxij, dfs_dyij, dfs_dzij,fsij_acc)
!$ACC& CREATE(Gi_acc, Gi_3D_acc, dfuN_acc)
!$ACC& CREATE(dANN_dxij, dE_dp)
!$ACC& CREATE(Bb1Sij_, Pb1Sij_, Bb1Sji_, Pb1Sji_, Bb3Sij_, Pb3Sij_)
!$ACC& CREATE(F2_rij_BOP_,F3_rij_BOP_,FCxij_BOP_,FCyij_BOP_,FCzij_BOP_)
!$ACC& CREATE(dE_dparam_list, dfs_dxij, dfs_dyij, dfs_dzij)
!$ACC& CREATE(dE7Sij_list_, dBZij_list_, PsidFij_list_)
!$ACC& CREATE(FSBb3C2_ij_list_, FSPb3C2_ij_list_, U3_of)
!$ACC& COPYOUT(BZij_list_, BSij_list_, Psi_of, Ep_of, frr)
!$ACC& COPYOUT(Ep_of, frr)

      call dANN_ACC(max2_nbrs,Rc15) ! aladyn_pi_PINN_ACC.f, get dANN_dxij 
      call dBOP_ACC(ecohe)          ! Differentiates BOP part !

! --- Calc Actual Force Vectors ---

! $ACC PARALLEL LOOP
!$ACC KERNELS LOOP WORKER
      do ni=1,natoms !  ! Loop III: Full forces !

       fr1 = 0.d0; fr2 = 0.d0; fr3 = 0.d0

!$ACC LOOP VECTOR REDUCTION(+: fr1,fr2,fr3)
       do nb1 = 1,max2_nbrs     ! Loop I: i - j bond !
        nj = nbr_list(nb1,ni)

        xr_ij = x_ij_dev(nb1,ni) ! x/r !
        yr_ij = y_ij_dev(nb1,ni) ! y/r !
        zr_ij = z_ij_dev(nb1,ni) ! z/r !

        i_at_j = nbr_list_inv(nb1,ni)

        Sr_ij = F3_rij_BOP_(nb1,ni)
        Sr_ji = F3_rij_BOP_(i_at_j,nj)

        Sigma_Psi_i = U3_of(nsigma_bop,ni)/sqrt(Psi_of(ni))
                  ! = sigma_i*(Psi(i)^-1/2) !
        Sigma_Psi_j = U3_of(nsigma_bop,nj)/sqrt(Psi_of(nj))
                   ! = sigma_j*(Psi(j)^-1/2) !

        ! Loop of parameter count ip= 1,.. 8 !
        Fa_ij1=0.d0; Fa_ij2=0.d0; Fa_ij3=0.d0
!$ACC LOOP SEQ
        do ip = 1, nBOP_params
         Fa_ij1 = Fa_ij1 +
     =   dE_dparam_list(ip,ni)*dANN_dxij(1,ip,nb1,ni) - 
     -   dE_dparam_list(ip,nj)*dANN_dxij(1,ip,i_at_j,nj)
         Fa_ij2 = Fa_ij2 +
     =   dE_dparam_list(ip,ni)*dANN_dxij(2,ip,nb1,ni) - 
     -   dE_dparam_list(ip,nj)*dANN_dxij(2,ip,i_at_j,nj)
         Fa_ij3 = Fa_ij3 +
     =   dE_dparam_list(ip,ni)*dANN_dxij(3,ip,nb1,ni) - 
     -   dE_dparam_list(ip,nj)*dANN_dxij(3,ip,i_at_j,nj)
        enddo

        Fw_ij = Sigma_Psi_i*PsidFij_list_(nb1,ni)
            ! = sigma_i*Sij*bij*dFc*(Psi(i)^-1/2) !
        Fw_ji = Sigma_Psi_j*PsidFij_list_(i_at_j,nj)
            ! = sigma_j*Sji*bji*dFc*(Psi(j)^-1/2) !

        Fr_ij = 0.5d0*(F2_rij_BOP_(nb1,ni) + F2_rij_BOP_(i_at_j,nj) - 
     -                (Fw_ij + Fw_ji))
        Fr_Sr_ij_ji = (Fr_ij - 0.5d0*(Sr_ij+Sr_ji))
        ff1 = Fr_Sr_ij_ji*xr_ij -
     -  0.5d0*(FCxij_BOP_(nb1,ni)-FCxij_BOP_(i_at_j,nj)) + Fa_ij1
        ff2 = Fr_Sr_ij_ji*yr_ij -
     -  0.5d0*(FCyij_BOP_(nb1,ni)-FCyij_BOP_(i_at_j,nj)) + Fa_ij2
        ff3 = Fr_Sr_ij_ji*zr_ij -
     -  0.5d0*(FCzij_BOP_(nb1,ni)-FCzij_BOP_(i_at_j,nj)) + Fa_ij3

        fr1 = fr1 + ff1   ! x force !
        fr2 = fr2 + ff2   ! y force !
        fr3 = fr3 + ff3   ! z force !

       enddo ! do nb1 = 1,max2_nbrs !

       frr(1,ni) = fr1; frr(2,ni) = fr2; frr(3,ni) = fr3

      enddo ! do ni=1,natoms !  ! Loop III: no stress !


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

      if (allocated(FcSij_)) deallocate(FcSij_, stat=ialloc(14))
      if (allocated(dFcSij_)) deallocate(dFcSij_, stat=ialloc(15))
      if (allocated(Bb1Sij_)) deallocate(Bb1Sij_, stat=ialloc(16))
      if (allocated(Pb1Sij_)) deallocate(Pb1Sij_, stat=ialloc(17))
      if (allocated(Bb1Sji_)) deallocate(Bb1Sji_, stat=ialloc(18))
      if (allocated(Pb1Sji_)) deallocate(Pb1Sji_, stat=ialloc(19))
      if (allocated(Bb3Sij_)) deallocate(Bb3Sij_, stat=ialloc(20))
      if (allocated(Pb3Sij_)) deallocate(Pb3Sij_, stat=ialloc(21))
      if (allocated(PsidFij_dev)) deallocate(PsidFij_dev, 
     1 stat=ialloc(22))
      if (allocated(dE_dp)) deallocate(dE_dp, stat=ialloc(23))
      if (allocated(dE7Sij_list_)) deallocate(dE7Sij_list_, 
     1 stat=ialloc(24))
      if (allocated(dBZij_list_)) deallocate(dBZij_list_, 
     1 stat=ialloc(25))
      if (allocated(PsidFij_list_)) deallocate(PsidFij_list_, 
     1 stat=ialloc(26))
      if (allocated(FSBb3C2_ij_list_)) deallocate(FSBb3C2_ij_list_, 
     1 stat=ialloc(27))
      if (allocated(FSPb3C2_ij_list_)) deallocate(FSPb3C2_ij_list_, 
     1 stat=ialloc(28))
      if (allocated(BSij_list_)) deallocate(BSij_list_, stat=ialloc(29))
      if (allocated(BZij_list_)) deallocate(BZij_list_, stat=ialloc(30))
      if (allocated(F2_rij_BOP_)) deallocate(F2_rij_BOP_, 
     1 stat=ialloc(31))
      if (allocated(F3_rij_BOP_)) deallocate(F3_rij_BOP_, 
     1 stat=ialloc(32))
      if (allocated(FCxij_BOP_)) deallocate(FCxij_BOP_, stat=ialloc(33))
      if (allocated(FCyij_BOP_)) deallocate(FCyij_BOP_, stat=ialloc(34))
      if (allocated(FCzij_BOP_)) deallocate(FCzij_BOP_, stat=ialloc(35))

      ierr = 0
      do i=1,36
       ierr = ierr + ialloc(i)
      enddo
      if(ierr.ne.0) then
       write(6,*)'ERROR deallocating x in Frc_PINN_ACC'
       ihalt = 1
      endif

      return
      end subroutine      ! Frc_PINN_ACC  !
!
! ------------------------------------------------------------------
! ------------------------------------------------------------------
!
      END MODULE  ! PINN_ACC !
