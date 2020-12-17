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
! Use a trained Physically Informed NN to get pot. param for a specific 
! potential form
!------------------------------------------------------------------
!
      MODULE PINN

      use sys_OMP
      use constants
      use sim_box
      use pot_module
      use ANN
      save

      integer, dimension(:), allocatable :: nbrs_all_of,nbrs_2nd_of,
     & id_nbr
      !dir$ attributes align:64 :: nbrs_all_of
      !dir$ attributes align:64 :: nbrs_2nd_of
      !dir$ attributes align:64 :: id_nbr

      integer, dimension(:,:), allocatable :: nbr_list_inv
      !dir$ attributes align:64 :: nbr_list_inv

      double precision, dimension(:), allocatable :: Sijn,Zijn,FcSij_ni,
     1 VB_of, Psi_of, Psi_j, bop_lambda_of, bop_sigma_of
      !dir$ attributes align:64 :: VB_of
      !dir$ attributes align:64 :: Psi_of
      !dir$ attributes align:64 :: Psi_J

      double precision, dimension(:,:), allocatable :: BSij_list,
     1 BZij_list, dBZij_list, dE7Sij_list, U3_of_ni
      !dir$ attributes align:64 :: BSij_list
      !dir$ attributes align:64 :: BZij_list
      !dir$ attributes align:64 :: dBZij_list
      !dir$ attributes align:64 :: dE7Sij_list
      !dir$ attributes align:64 :: U3_of_ni

      double precision, dimension(:), allocatable :: rij_ni_nbr
      double precision, dimension(:,:), allocatable :: 
     1 FSBb3C2_ij_list, FSPb3C2_ij_list, PsidFij_list, dE_dparam_list
      !dir$ attributes align:64 :: FSBb3C2_ij_list
      !dir$ attributes align:64 :: FSPb3C2_ij_list
      !dir$ attributes align:64 :: PsidFij_list ! Sij*bij*dFc !
      !dir$ attributes align:64 :: dE_dparam_list

      double precision, dimension(:,:), allocatable :: F2_rij_BOP,
     1 F3_rij_BOP
      !dir$ attributes align:64 :: F2_rij_BOP
      !dir$ attributes align:64 :: F3_rij_BOP
      double precision, dimension(:,:,:), allocatable :: FCxij_BOP
      !dir$ attributes align:64 :: FCxij_BOP

      integer, parameter :: nA_bop = 1
      integer, parameter :: nalpha_bop = 2
      integer, parameter :: nB_bop = 3
      integer, parameter :: nbeta_bop = 4
      integer, parameter :: nh_bop = 5
      integer, parameter :: nsigma_bop = 6
      integer, parameter :: naa_bop = 7
      integer, parameter :: nlambda_bop = 8
      
      CONTAINS
!
! ---------------------------------------------------------------------
! Calculates Analytical derivatives of the BOP function.
! ---------------------------------------------------------------------
!
      subroutine dBOP_OMP(ecohe)

      use atoms

      implicit double precision (a-h,o-z)

      double precision, intent(out) :: ecohe
  
      double precision, dimension(3,nbrs_per_atom) :: xr_ij_ni,
     1 x_ij_ni, xij_
      double precision, dimension(nbrs_per_atom) :: rij_
      double precision, dimension(nbrs1_per_atom) :: 
     1 FSBb3C2_ij, FSPb3C2_ij, FSBb3C2_ji, FSPb3C2_ji,
     1 Bb1Sij_, Pb1Sij_, Bb1Sji_, Pb1Sji_, Bb3Sij_, Pb3Sij_,dE7Sij_
      double precision, dimension(nbrs_per_atom) :: dFcSij_ni,
     1 dZijn
      double precision, dimension(3) :: ff, xr_ij, xr_ik, xr_ik_ij,
     1 xr_ij_ik, dcos_ijk, x_ij, x_ik, x_jk

      double precision, dimension(3) :: dE1bC, dE2bC, dW2bC
      double precision, dimension(Nodes_of_layer(net_layers)) :: 
     1 dE_dparam

      Rc = Rc_ann
      Rc15 = 1.5d0*Rc
      Sqr_Rc15 = Rc15**2
      d4 = d4_ann

      ecohe = 0.d0

      h11 = h(1,1); h12 = h(1,2); h13 = h(1,3)
      h22 = h(2,2); h23 = h(2,3); h33 = h(3,3)
 
!$OMP  PARALLEL DO DEFAULT(NONE) REDUCTION(+:ecohe)
!$OMP& PRIVATE(ni, nbrs_1st, nbrs_2nd, nj, nb1,nb2)
!$OMP& PRIVATE(x_ij,rij,rij1, x_ik,rik, x_jk,rjk, sx0,sy0,sz0)
!$OMP& PRIVATE(cos_ijk, j,k,l, indk, rrc_ijk, rrr,rrc,rrc3,rrc4, Fc_ij)
!$OMP& PRIVATE(Exp_blr_ijk, Sijk_term, Sij1, dE7Sij, dE7Sijk, Fc_ijk)
!$OMP& PRIVATE(rrc_ikj, Fc_ikj, Exp_blr_ikj, Sikj_term, dE7Sikj, r2ij)
!$OMP& PRIVATE(Zij1, dZij1, x_ij_ni, bl, U3)
!$OMP& PRIVATE(b_ij, Psi_ij, Ei, Psi, VB_of_i, xr_ij,xr_ik)
!$OMP& PRIVATE(Sijn,Fcijn, Zijn,dZijn, FcSij_ni, dE7Sij_)
!$OMP& PRIVATE(FcSij, VA_ij,VB_ij, VB_ij_exp, U_ij, Epot, xr_ij_ni)
!$OMP& PRIVATE(cosTh_h2, FcSik, Sij_term, dFc_ij, dVA_ij,dVB_ij, dU_ij)
!$OMP& SHARED(natoms, nbrs_all_of, nbrs_2nd_of, nbr_list)
!$OMP& SHARED(sx,sy,sz, h11,h22,h33,h12,h13,h23, Rc, d4)
!$OMP& SHARED(BSij_list,BZij_list)
!$OMP& SHARED(F2_rij_BOP, PsidFij_list, dBZij_list, dE7Sij_list) 
!$OMP& SHARED(Ep_of, U3_of, Psi_of, dfuN)
!$OMP& SCHEDULE(DYNAMIC, CHUNK)

      do ni=1,natoms

       Ei = 0.d0
       Psi = 0.d0
       VB_of_i = 0.d0

       nbrs_1st = nbrs_all_of(ni)
       nbrs_2nd = nbrs_2nd_of(ni)

       U3(:) = U3_of(:,ni)

       do nb1 = 1,nbrs_2nd
        nj = nbr_list(nb1,ni)

        sx0 = sx(nj) - sx(ni)
        sy0 = sy(nj) - sy(ni)
        sz0 = sz(nj) - sz(ni)
        sx0 = sx0 - dnint(sx0)     ! make periodic along X !
        sy0 = sy0 - dnint(sy0)     ! make periodic along Y !
        sz0 = sz0 - dnint(sz0)     ! make periodic along Z !

        x_ij(1) = h11*sx0 + h12*sy0 + h13*sz0
        x_ij(2) = h22*sy0 + h23*sz0
        x_ij(3) = h33*sz0

        r2ij = x_ij(1)**2 + x_ij(2)**2 + x_ij(3)**2
        rij = sqrt(r2ij)
        rij1 = 1.0/rij

        xr_ij(:) = x_ij(:)*rij1
        x_ij_ni(:,nb1) = x_ij(:)
        xr_ij_ni(:,nb1) = xr_ij(:)
       enddo ! do nb1 = 1,nbrs_2nd !

       do nb1 = 1,nbrs_1st     ! Loop I: i - j bond !

        x_ij(:) = x_ij_ni(:,nb1)
        rij = sqrt(x_ij(1)**2 + x_ij(2)**2 + x_ij(3)**2)

        Sij1 = 1.d0
        dE7Sij = 0.d0

        do nb2 = 1,nbrs_2nd   ! Screening: i-k bond !
          if(nb2.ne.nb1) then
         x_ik(:) = x_ij_ni(:,nb2)
         rik = sqrt(x_ik(1)**2 + x_ik(2)**2 + x_ik(3)**2)
         x_jk(:) = x_ik(:) - x_ij(:)
         rjk = sqrt(x_jk(1)**2 + x_jk(2)**2 + x_jk(3)**2)

         rrr = rik + rjk - rij
         if(rrr.lt.Rc) then
          bl = U3(nlambda_bop)
          rrc_ijk = Rc-rrr
          rrc4 = rrc_ijk**4
          Fc_ijk = rrc4 / (d4+rrc4)
          Exp_blr_ijk = exp(-bl*rrr)
          Sijk_term = 1.d0 - Fc_ijk*Exp_blr_ijk
          Sij1 = Sij1 * Sijk_term

          dE7Sijk = Fc_ijk*Exp_blr_ijk*rrr / Sijk_term
          dE7Sij = dE7Sij + dE7Sijk
         endif
          endif ! if(nb2.ne.nb1)... !
        enddo ! do nb2 = 1,nbrs_2nd   ! Screening: i-k bond !
 
        rrc = max(Rc-rij,0.d0)
        rrc4 = rrc**4
        Fc_ij = rrc4 / (d4+rrc4)

        Sijn(nb1) = Sij1
        FcSij = Sij1*Fc_ij
        FcSij_ni(nb1) = FcSij

        dE7Sij_(nb1) = FcSij*dE7Sij

       ! Sijn(nb1) is done, but not Sijn(nb > nb1) !

       enddo ! do nb1 = 1,nbrs_1st   ! Screening: i - j bond !

c --- Start another loop over neighbors of ni: Chem i - j bond ---

       Zijn(1:nbrs_1st) = 0.d0
       dZijn(1:nbrs_1st) = 0.d0

       do nb1 = 1,nbrs_1st
        x_ij(:) = x_ij_ni(:,nb1)  ! xij !
        rij = sqrt(x_ij(1)**2 + x_ij(2)**2 + x_ij(3)**2)
        rij1 = 1.d0/rij

        FcSij = FcSij_ni(nb1)

        Zij1 = 0.d0; dZij1 = 0.d0

        do nb2 = nb1+1,nbrs_1st   ! Three body: i-k bond !

         x_ik(:) = xr_ij_ni(:,nb2)  ! xik/rik !

       cos_ijk = 
     1 (x_ij(1)*x_ik(1) + x_ij(2)*x_ik(2) + x_ij(3)*x_ik(3)) * rij1
       cosTh_h2 = (cos_ijk - U3(nh_bop))**2

       ! Zij1 = Zij1 + Fcik*Sik*a_ijk*cos_ijk !
       ! dZij1 = dZij1 + Fcik*Sik*cos_ijk*da_ijk !
       FcSik = FcSij_ni(nb2)
       Zij1 = Zij1 + U3(naa_bop)*FcSik*cosTh_h2
       Zijn(nb2) = Zijn(nb2) + U3(naa_bop)*FcSij*cosTh_h2
       dZij1 = dZij1 + FcSik*cosTh_h2
       dZijn(nb2) = dZijn(nb2) + FcSij*cosTh_h2
        enddo ! do nb2 = nb1+1,nbrs_1st   ! Three body: i-k bond !

        Zijn(nb1) = Zijn(nb1) + Zij1

        ! 2-body potential energy terms !

        rrc = max(Rc-rij,0.d0)
        rrc4 = rrc**4
        Fc_ij = rrc4 / (d4+rrc4)

        b_ij = 1.d0/sqrt(1.d0 + Zijn(nb1))
        Psi_ij = Sijn(nb1)*b_ij  ! Sij*bij !

        Psi = Psi + Fc_ij*Psi_ij       ! e_i = \sum_j(Fc*Sij*bij) !

        VA_ij =  ! Aij*exp(-alpha*rij) !
     =  exp(U3(nA_bop)-U3(nalpha_bop)*rij)

        VB_ij_exp =  ! Bij*exp(-beta.rij) !
     =  exp(U3(nB_bop)-U3(nbeta_bop)*rij)
        VB_ij = VB_ij_exp * Psi_ij ! Fc*Sij*bij*Bij*exp(-beta.rij) !

        U_ij = VA_ij - VB_ij
        Ei = Ei + U_ij*Fc_ij

        VB_of_i = VB_of_i + VB_ij*Fc_ij   ! Fc*VB     !

        Sij_term = Sijn(nb1)

        BSij_list(nb1,ni) = Sij_term   ! store Sij only !
        BZij_list(nb1,ni) = Zijn(nb1)     ! store Zij !

        rrc3 = rrc**3
        dFc_ij = -4.d0*rrc3*d4 / (d4+rrc4)**2    ! dFc/dr !

        dVA_ij = -U3(nalpha_bop)*VA_ij
        dVB_ij = -U3(nbeta_bop)*VB_ij
        dU_ij = dVA_ij - dVB_ij
        ! 2 body dU_ij at fixed Sij, bij(zij) !

        F2_rij_BOP(nb1,ni) = (dU_ij*Fc_ij + U_ij*dFc_ij)
        PsidFij_list(nb1,ni) = Psi_ij*dFc_ij ! Sij*bij*dFc/dr !

        dZijn(nb1) = dZijn(nb1) + dZij1
        dBZij_list(nb1,ni) = dZijn(nb1)
        dE7Sij_list(nb1,ni) = dE7Sij_(nb1)

       enddo ! do nb1 = 1,nbrs_1st   ! Screening: i - j bond !

       Epot = 0.5d0*Ei - U3(nsigma_bop)*sqrt(Psi)

       ecohe = ecohe + Epot

       Ep_of(ni) = 2.d0*Epot  ! Twice Epot !
       Psi_of(ni) = max(Psi,1.d-8)

      enddo ! do ni=1,natoms !

!$OMP END PARALLEL DO

C     do ni=1,5
C      write(6,*)'Ep_of(',ni,')=',0.5d0*Ep_of(ni), Psi_of(ni)
C     enddo
C     write(6,*)'1:ecohe=',ecohe/natoms

c
c ------------- Force start here --------------
c
c ------------- Start dU/dr: Loop I: --------------
c
!$OMP PARALLEL DO DEFAULT(NONE)
!$OMP& PRIVATE(ni,nb1, nbrs_1st, nbrs_2nd, nj,nbrs1_of_j, i_at_j,ntp_k)
!$OMP& PRIVATE(ntp_i,ntp_i1,ntp_i2,ntp_j,ntp_j1,ntp_j2,ntp_ijj,ntp_jii)
!$OMP& PRIVATE(nbk,k, ntp_k1,ntp_k2, ntp_ikj,ntp_ijk, ntp_kij,ntp_jik)
!$OMP& PRIVATE(id_nbr, Sigma_Psi_i,Sigma_Psi_j, rrr)
!$OMP& PRIVATE(x_ij,rij,rij1, xr_ij, xij_,rij_, sx0,sy0,sz0)
!$OMP& PRIVATE(x_ik,rik,rik1, xr_ik, x_jk,rjk,r2jk)
!$OMP& PRIVATE(rrc_ij,rrc3,rrc4, Fcij,dFcij, Fb1Sij,Fb1Sji,Fb3Sij)
!$OMP& PRIVATE(VB_ij_exp,VB_ji_exp, Bb1ij,Pb1ij, Bb1Sij,Bb3Sij,Pb3Sij)
!$OMP& PRIVATE(rrc, S_ij,S_ji, Z_ij,Z_ji, b_ij,b_ji, FcSij, Pb1Sij,dFc)
!$OMP& PRIVATE(FcSij_ni,dFcSij_ni,FcbijSij,FcbijSji,Fc_rrr)
!$OMP& PRIVATE(Bb1Sij_,Bb1Sji_,Bb1Sji, Pb1Sij_,Pb1Sji_,Bb3Sij_,Pb3Sij_)
!$OMP& PRIVATE(dE1bC,dE2bC,dW2bC, xr_ij_ik,xr_ik_ij)
!$OMP& PRIVATE(dSr_ij,dEsr_ij,dWsr_ij, Bb3C2_ij,dWbfr_ij, FcSik)
!$OMP& PRIVATE(dE7Sij,dE7Sij_, VA_ij_exp, VA_ij_Fc, Pb1Sji)
!$OMP& PRIVATE(bl_ijk,bl_ikj,bl_jik,bl_kij, aa_ikj,aa_ijk, h_ikj,h_ijk)
!$OMP& PRIVATE(Exp_blr_ijk, Exp_blr_ikj, Exp_blr_jik, Exp_blr_kij)
!$OMP& PRIVATE(Sijk_term, Sjik_term, Sikj_term, Skij_term)
!$OMP& PRIVATE(dS_ijk,dS_jik, dS_ikj,dS_kij, Pb1dSr_ikj,Pb1dSr_kij)
!$OMP& PRIVATE(Bb1Sik,Bb1Ski, Pb1Sik,Pb1Ski, Bb1dSr_ikj,Bb1dSr_kij)
!$OMP& PRIVATE(Bb1dSr_ijk, Bb1dSr_jik, Pb1dSr_ijk, Pb1dSr_jik)
!$OMP& PRIVATE(cos_ikj,cos_h_ikj,cos2_h_ikj, cos_h_ijk,cos2_h_ijk)
!$OMP& PRIVATE(C1ijk,C1ikj,C2ijk,C2ikj, Bb3Sik,Pb3Sik, dE7Sik)
!$OMP& PRIVATE(PiFb3SijC1, PiFb3SikC1, PiFb3SijC2,PiFb3SikC2)
!$OMP& PRIVATE(BexpFb3SijC1, BexpFb3SikC1, BexpFb3SikC2, BexpFb3SijC2)
!$OMP& PRIVATE(FSBb3C2ij, FSPb3C2ij, dEbS3r_ij, dWbS3r_ij, dEWbfr_ij)
!$OMP& PRIVATE(dE_dparam, bl,dE6,dE7, dW6,dW7,dE6_,dE7_,dW6_,dW7_)
!$OMP& PRIVATE(U3, U3_of_ni)
!$OMP& SHARED(natoms, Rc, Rc15, Sqr_Rc15, d4)
!$OMP& SHARED(nbrs_all_of, nbrs_2nd_of, nbr_list_inv)
!$OMP& SHARED(sx,sy,sz, h11,h22,h33,h12,h13,h23, nbr_list)
!$OMP& SHARED(BSij_list, BZij_list, Psi_of, U3_of, FCxij_BOP)
!$OMP& SHARED(FSBb3C2_ij_list, FSPb3C2_ij_list, F3_rij_BOP)
!$OMP& SHARED(dBZij_list, dE7Sij_list, dE_dparam_list) 
!$OMP& SCHEDULE(DYNAMIC,CHUNK)

      do ni=1,natoms

       nbrs_1st = nbrs_all_of(ni)
       nbrs_2nd = nbrs_2nd_of(ni)  ! rij < 1.5Rc !

       ! get BOP parameters for the current atom i !

       U3(:) = U3_of(:,ni)  ! Get BOP params for each atom !

       if(nbrs_2nd.gt.0) then
        Sigma_Psi_i = U3(nsigma_bop)/sqrt(Psi_of(ni))
                  ! = sigma_i*(Psi(i)^-1/2) !
       else
        Sigma_Psi_i = 0.d0
       endif

       dE_dparam(:) = 0.d0

       do nb1 = 1,nbrs_2nd     ! Loop 0: i - j bond; rij < 1.5Rc !
        nj = nbr_list(nb1,ni)  ! VECTORIZED: speedup: 1.550      !
        id_nbr(nb1) = nj

        sx0 = sx(nj) - sx(ni)
        sy0 = sy(nj) - sy(ni)
        sz0 = sz(nj) - sz(ni)
        sx0 = sx0 - dnint(sx0)     ! make periodic along X !
        sy0 = sy0 - dnint(sy0)     ! make periodic along Y !
        sz0 = sz0 - dnint(sz0)     ! make periodic along Z !

        x_ij(1) = h11*sx0 + h12*sy0 + h13*sz0
        x_ij(2) = h22*sy0 + h23*sz0
        x_ij(3) = h33*sz0

        xij_(:,nb1) = x_ij(:)
        rij = sqrt(x_ij(1)**2 + x_ij(2)**2 + x_ij(3)**2)
        rij_(nb1) = rij

         ! get BOP parameters for neighbor j !
         ! Only these are needed !
      U3_of_ni(nB_bop,nb1) = U3_of(nB_bop,nj) 
      U3_of_ni(nbeta_bop,nb1) = U3_of(nbeta_bop,nj)
      U3_of_ni(nh_bop,nb1) = U3_of(nh_bop,nj)
      U3_of_ni(nsigma_bop,nb1) = U3_of(nsigma_bop,nj)
      U3_of_ni(nlambda_bop,nb1)= U3_of(nlambda_bop,nj)

         nbrs1_of_j = nbrs_2nd_of(nj)

          if(nb1.le.nbrs_1st) then
         rrc_ij = max(Rc-rij,0.d0)
         rrc3 = rrc_ij**3
         rrc4 = rrc_ij**4
         Fcij = rrc4 / (d4+rrc4)            ! Fc(rij) !
         dFcij = -4.d0*rrc3*d4 / (d4+rrc4)**2     ! dFc/dr !

         i_at_j = nbr_list_inv(nb1,ni)        ! get i as a nbr of j !

         S_ij = BSij_list(nb1,ni)       ! Sij !
         S_ji = BSij_list(i_at_j,nj)    ! Sji !
         Z_ij = BZij_list(nb1,ni)       ! Zij !
         Z_ji = BZij_list(i_at_j,nj)    ! Zji !
         b_ij = 1.d0/sqrt(1.d0 + Z_ij)  ! bij !
         b_ji = 1.d0/sqrt(1.d0 + Z_ji)  ! bji !

         FcSij = Fcij*S_ij                        ! Fc(rij)*Sij  !
         FcSij_ni(nb1) = FcSij
         dFcSij_ni(nb1) = dFcij*S_ij              ! dFc(rij)*Sij !
         Fb1Sij = FcSij*b_ij                      ! Fc*Sij*bij   !
         Fb1Sji = Fcij*S_ji*b_ji                  ! Fc*Sji*bji   !
         Fb3Sij = FcSij*b_ij**3                   ! Fc*Sij*bij^3 !
         dE7Sij = dE7Sij_list(nb1,ni)
         dE7Sij_(nb1) = dE7Sij

         if(nbrs1_of_j.gt.0) then
          Sigma_Psi_j = U3_of_ni(nsigma_bop,nb1)/sqrt(Psi_of(nj))
                    ! = sigma_j*(Psi(j)^-1/2) !
         else
          Sigma_Psi_j = 0.d0
         endif

         VA_ij_exp = exp(U3(nA_bop)-U3(nalpha_bop)*rij) ! Aij*exp(-alpha*rij)
         VB_ij_exp = exp(U3(nB_bop)-U3(nbeta_bop)*rij)  ! Bij*exp(-beta*rij)

         Bb1ij = VB_ij_exp*b_ij    ! Bij*exp(-beta.rij)*bij !
         Pb1ij = Sigma_Psi_i*b_ij  ! sigma_i*(Psi(i)^-1/2)*bij     !
         Bb1Sij = VB_ij_exp*Fb1Sij ! Bij*exp(-beta.rij)*Fc*Sij*bij !
         Bb1Sij_(nb1) = Bb1Sij
         Pb1Sij_(nb1) = Sigma_Psi_i*Fb1Sij 
                      ! sigma_i*(Psi(i)^-1/2)*Fc*Sij*bij !

         ! Bji*exp(-beta.rij) !
         VB_ji_exp = 
     =   exp(U3_of_ni(nB_bop,nb1) - U3_of_ni(nbeta_bop,nb1)*rij)
         Bb1Sji_(nb1) = VB_ji_exp*Fb1Sji ! Bji*exp(-beta.rij)*FcSji*bji 
         Pb1Sji_(nb1) = Sigma_Psi_j*Fb1Sji ! sigma_j*(Psi(j)^-1/2)*Fc*Sji*bji

         Bb3Sij = VB_ij_exp*Fb3Sij   
         Bb3Sij_(nb1) = Bb3Sij  ! Bij*exp(-beta.rij)*Fc*Sij*bij^3
         Pb3Sij = Sigma_Psi_i*Fb3Sij 
         Pb3Sij_(nb1) =  Pb3Sij ! sigma_i*(Psi(i)^-1/2)*Fc*Sij*bij^3

         VA_ij_Fc = VA_ij_exp*Fcij   ! Fc*Aij*exp(-alpha*rij) !
       dE_dparam(1) = dE_dparam(1) + VA_ij_Fc 
       dE_dparam(2) = dE_dparam(2) - VA_ij_Fc*rij 
       dE_dparam(3) = dE_dparam(3) - Bb1Sij
       dE_dparam(4) = dE_dparam(4) + Bb1Sij*rij
       dE_dparam(naa_bop) = dE_dparam(naa_bop) + 
     + 0.5d0*(Bb3Sij+Pb3Sij)*dBZij_list(nb1,ni)
       dE_dparam(nlambda_bop) = 
     = dE_dparam(nlambda_bop) - (Bb1ij+Pb1ij)*dE7Sij

          endif ! if(nb1.le.nbrs_1st)... !
       enddo ! do nb1 = 1,nbrs_2nd !

       dE_dparam(nsigma_bop) = -2.d0*sqrt(Psi_of(ni))  
                                    ! -2(Psi(i)^1/2) !
c
c --- Calc. dSij_dxi and dSji_dxi ---
c
       do nb1 = 1,nbrs_1st     ! Loop I: i - j bond; rij < 1.5Rc !
        nj = id_nbr(nb1)

        rij = rij_(nb1)
        rij1 = 1.d0/rij

        x_ij(:) = xij_(:,nb1)
        xr_ij(:) = x_ij(:)*rij1
c
c --- Loop i - k bond to calc the 1/Sijk and dSij_dxi terms ---
c
        FcSij = FcSij_ni(nb1)

        Bb1Sij = Bb1Sij_(nb1)   ! Bij*exp(-beta.rij)*Fc*Sij*bij
        Pb1Sij = Pb1Sij_(nb1)   ! sigma_i*(Psi(i)^-1/2)*Fc*Sij*bij

        Bb1Sji = Bb1Sji_(nb1)   ! Bji*exp(-beta.rij)*Fc*Sji*bji
        Pb1Sji = Pb1Sji_(nb1)   ! sigma_j*(Psi(j)^-1/2)*Fc*Sji*bji

        Bb3Sij = Bb3Sij_(nb1)   ! Bij*exp(-beta.rij)*Fc*Sij*bij^3
        Pb3Sij = Pb3Sij_(nb1)   ! sigma_i*(Psi(i)^-1/2)*Fc*Sij*bij^3

        dE7Sij = dE7Sij_(nb1)

        dSr_ij = 0.d0; dEsr_ij = 0.d0; dWsr_ij = 0.d0
        Bb3C2_ij = 0.d0; dWbfr_ij = 0.d0
        dE1bC(:) = 0.d0; dE2bC(:) = 0.d0; dW2bC(:) = 0.d0
        dE6 = 0.d0; dW6 = 0.d0
        dE7 = 0.d0; dW7 = 0.d0

        do nbk = 1,nbrs_1st  ! rij < 1.5Rc; rik < Rc !
         if(nbk.ne.nb1) then

         k = id_nbr(nbk)

         FcSik = FcSij_ni(nbk)
         dE7Sik = dE7Sij_(nbk)

         x_ik(:) = xij_(:,nbk)

         x_jk(:) = xij_(:,nbk) - x_ij(:)
         r2jk = x_jk(1)**2 + x_jk(2)**2 + x_jk(3)**2

         rjk = sqrt(r2jk)
         rik = rij_(nbk)
         rik1= 1.d0/rik

         xr_ik(:) = x_ik(:)*rik1

         if(r2jk.lt.Sqr_Rc15) then   ! rjk < 1.5Rc !

          Bb1Sik = Bb1Sij_(nbk) ! Bik*exp(-beta.rik)*Fc*Sik*bik !
          Pb1Sik = Pb1Sij_(nbk) ! sigma_i*(Psi(i)^-1/2)*Fc*Sik*bik !

          Bb1Ski = Bb1Sji_(nbk) ! Bki*exp(-beta.rik)*Fc*Ski*bki ! 
          Pb1Ski = Pb1Sji_(nbk) ! sigma_k*(Psi(k)^-1/2)*Fc*Ski*bki ! 

          rrr = rij + rjk - rik
           if(rrr.lt.Rc) then
          bl_ikj = U3(nlambda_bop)
          bl_kij = U3_of_ni(nlambda_bop,nbk)
          rrc = max(Rc-rrr,0.d0)
          rrc3 = rrc**3
          rrc4 = rrc*rrc3
          Fc_rrr = rrc4 / (d4+rrc4)   ! Fc(r_ikj) = Fc(r_kij) !
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

           else
          Bb1dSr_ikj = 0.d0
          Bb1dSr_kij = 0.d0
          Pb1dSr_ikj = 0.d0
          Pb1dSr_kij = 0.d0
           endif ! if(rrr.lt.Rc)... !

          rrr = rik + rjk - rij
           if(rrr.lt.Rc) then
          bl_ijk = U3(nlambda_bop)
          rrc = max(Rc-rrr,0.d0)
          rrc3 = rrc**3
          rrc4 = rrc*rrc3
          Fc_rrr = rrc4 / (d4+rrc4)
          Exp_blr_ijk = exp(-bl_ijk*rrr)
          Sijk_term = 1.d0 - Fc_rrr*Exp_blr_ijk  ! S_ijk !

          dFc = -4.d0*rrc3*d4 / (d4+rrc4)**2   ! dFc/dr !
          dS_ijk = (bl_ijk*Fc_rrr - dFc)*Exp_blr_ijk
          dS_ijk = dS_ijk/Sijk_term

          Bb1dSr_ijk = Bb1Sij*dS_ijk
          ! Bij*exp(-beta.rij)*Fc*Sij*bij*dS_ijk/Sijk_term !
          Pb1dSr_ijk = Pb1Sij*dS_ijk
          ! sigma_i*(Psi(i)^-1/2)*Fc*Sij*bij*dS_ijk/Sijk_term !

          dSr_ij = dSr_ij + dS_ijk

           else
          Bb1dSr_ijk = 0.d0
          Pb1dSr_ijk = 0.d0
           endif ! if(rrr.lt.Rc)... !

          dEsr_ij = dEsr_ij + Bb1dSr_ikj + Bb1dSr_kij - Bb1dSr_ijk
          dWsr_ij = dWsr_ij + Pb1dSr_ikj + Pb1dSr_kij - Pb1dSr_ijk

         endif ! if(r2jk.lt.Sqr_Rc15)... rjk < 1.5Rc !   

          aa_ikj = U3(naa_bop)
          h_ikj = U3(nh_bop)
          aa_ijk = U3(naa_bop)
          h_ijk = U3(nh_bop)

       cos_ikj = xr_ij(1)*xr_ik(1)+xr_ij(2)*xr_ik(2)+xr_ij(3)*xr_ik(3)
       cos_h_ikj = cos_ikj - h_ikj
       cos_h_ijk = cos_ikj - h_ijk
       cos2_h_ikj = cos_h_ikj**2
       cos2_h_ijk = cos_h_ijk**2

          xr_ik_ij(:) = (xr_ik(:) - xr_ij(:)*cos_ikj)*rij1
                      ! (xik/rik - xij/rij*cos_ikj)/rij !
          xr_ij_ik(:) = (xr_ij(:) - xr_ik(:)*cos_ikj)*rik1
                      ! (xij/rij - xik/rik*cos_ikj)/rik !
          C1ikj = aa_ikj*cos_h_ikj
          C1ijk = aa_ijk*cos_h_ijk
          C2ikj = aa_ikj*cos2_h_ikj
          C2ijk = aa_ijk*cos2_h_ijk

          Bb3Sik =  Bb3Sij_(nbk)
          BexpFb3SijC1 =  Bb3Sij*C1ijk
        ! Bij*exp(-beta.rij)*Fc*Sij*(bij**3)*aa_ijk*cos_h_ijk !
          BexpFb3SikC1 =  Bb3Sik*C1ikj
        ! Bik*exp(-beta.rik)*Fc*Sik*(bik**3)*aa_ikj*cos_h_ikj !
          BexpFb3SikC2 =  Bb3Sik*C2ikj
        ! Bik*exp(-beta.rik)*Fc*Sik*(bik**3)*a_ikj*cos2_h_ikj !
          BexpFb3SijC2 =  Bb3Sij*C2ijk
        ! Bij*exp(-beta.rij)*Fc*Sij*(bij**3)*a_ijk*cos2_h_ijk !

          Pb3Sik = Pb3Sij_(nbk)
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

          dE1bC(:) = dE1bC(:) + FcSik*C1ijk*xr_ik_ij(:)
        ! += Fc*Sik*aa_ijk*cos_h_ijk*(xik/rik - xij/rij*cos_ikj)/rij
          dE2bC(:) = dE2bC(:) + BexpFb3SikC1*xr_ik_ij(:)
        ! += Bik*exp(-beta.rik)*Fc*Sik*(bik**3)*aa_ikj*cos_h_ikj*
        !   *(xik/rik - xij/rij*cos_ikj)/rij
          dW2bC(:) = dW2bC(:) + PiFb3SikC1*xr_ik_ij(:)
        ! += sigma_i*(Psi(i)^-1/2)*Fc*Sik*bik*aa_ijk*cos_h_ijk*
        !   *(xik/rik - xij/rij*cos_ikj)/rij
          dE6 = dE6 + BexpFb3SijC1*FcSik 
        ! Bij*exp(-beta.rij)*FcSij*(bij**3)*aa_ijk*cos_h_ijk*FcSik !
          dW6 = dW6 + PiFb3SijC1*FcSik 
        ! sigma_i*(Psi(i)^-1/2)*FcSij*(bij**3)*aa_ijk*cos_h_ijk*FcSik !
          dE7 = dE7 + BexpFb3SijC2*dE7Sik
          dW7 = dW7 + PiFb3SijC2*dE7Sik

         endif ! if(nbk.ne.nb1)... !
        enddo ! do nbk = 1,nbrs_1st  ! rij < 1.5Rc; rik < Rc !

        do nbk = nbrs_1st+1,nbrs_2nd  ! rij < 1.5Rc; Rc < rik < 1.5Rc !
         x_ik(:) = xij_(:,nbk)
         x_jk(:) = xij_(:,nbk) - x_ij(:)
         r2jk = x_jk(1)**2 + x_jk(2)**2 + x_jk(3)**2

         if(r2jk.lt.Sqr_Rc15) then   ! rjk < 1.5Rc !

          rjk = sqrt(r2jk)
          rik = rij_(nbk)

          rrr = rik + rjk - rij
           if(rrr.lt.Rc) then
          bl_ijk = U3(nlambda_bop)
          rrc = max(Rc-rrr,0.d0)
          rrc3 = rrc**3
          rrc4 = rrc*rrc3
          Fc_rrr = rrc4 / (d4+rrc4)
          Exp_blr_ijk = exp(-bl_ijk*rrr)
          Sijk_term = 1.d0 - Fc_rrr*Exp_blr_ijk  ! S_ijk !

          dFc = -4.d0*rrc3*d4 / (d4+rrc4)**2   ! dFc/dr !
          dS_ijk = (bl_ijk*Fc_rrr - dFc)*Exp_blr_ijk
          dS_ijk = dS_ijk/Sijk_term

          Bb1dSr_ijk = Bb1Sij*dS_ijk
          ! Bij*exp(-beta.rij)*Fc*Sij*bij*dS_ijk/Sijk_term !
          Pb1dSr_ijk = Pb1Sij*dS_ijk
          ! sigma_i*(Psi(i)^-1/2)*Fc*Sij*bij*dS_ijk/Sijk_term !

          dSr_ij = dSr_ij + dS_ijk
          dEsr_ij = dEsr_ij - Bb1dSr_ijk
          dWsr_ij = dWsr_ij - Pb1dSr_ijk
           endif ! if(rrr.lt.Rc)... !

         endif ! if(r2jk.lt.Sqr_Rc15)... rjk < 1.5Rc !   
        enddo ! do nbk = nbrs_1st+1,nbrs_2nd ! rij<1.5Rc; Rc<rik<1.5Rc !

        FSBb3C2ij = FcSij*Bb3C2_ij ! Fc*Sij*Bb3C2_ij !
        FSPb3C2ij = FcSij*dWbfr_ij ! Fc*Sij*Pb3C2_ij !

        FSBb3C2_ij_list(nb1,ni) = FSBb3C2ij
! Fij*Sij*Sum_k{Bik*exp(-beta.rik)*Fc*Sik*bik**3*a_ikj*cos2_h_ikj} !
        FSPb3C2_ij_list(nb1,ni) = FSPb3C2ij
! Fij*Sij*Sum_k{sigma_i*(Psi(i)^-1/2)*Fc*Sik*bik**3*a_ikj*cos2_h_ikj} !

        dEbS3r_ij = FSBb3C2ij*dSr_ij  ! 3rd term of dEbSr_ij !
        dWbS3r_ij = FSPb3C2ij*dSr_ij  ! 3rd term of dWbSr_ij !
        dEWbfr_ij = dFcSij_ni(nb1)*(Bb3C2_ij + dWbfr_ij)
        F3_rij_BOP(nb1,ni) = dEsr_ij + dWsr_ij - 
     -  0.5d0*(dEWbfr_ij - dEbS3r_ij - dWbS3r_ij) 

        FCxij_BOP(:,nb1,ni) = -((Bb3Sij+Pb3Sij)*dE1bC(:) + 
     +                   FcSij*(dE2bC(:) + dW2bC(:)))

        dE_dparam(nh_bop) = dE_dparam(nh_bop) - (dE6+dW6)
        dE_dparam(nlambda_bop) = 
     =  dE_dparam(nlambda_bop) + 0.5d0*(dE7+dW7)

       enddo ! do nb1 = 1,nbrs_1st  ! rij < Rc !

       do nb1 = nbrs_1st+1,nbrs_2nd  ! Rc < rij < 1.5Rc !
        nj = id_nbr(nb1)

        rij = rij_(nb1)
        rij1 = 1.d0/rij

        x_ij(:) = xij_(:,nb1)
        xr_ij(:) = x_ij(:)*rij1
c
c --- Loop i - k bond to calc the 1/Sijk and dSij_dxi terms ---
c
        dSr_ij = 0.d0; dEsr_ij = 0.d0; dWsr_ij = 0.d0
        Bb3C2_ij = 0.d0; dWbfr_ij = 0.d0
        dE1bC(:) = 0.d0; dE2bC(:) = 0.d0; dW2bC(:) = 0.d0

        do nbk = 1,nbrs_1st  ! Rc < rij < 1.5Rc; rik < Rc !
         k = id_nbr(nbk)

         FcSik = FcSij_ni(nbk)
         dE7Sik = dE7Sij_(nbk)

         x_ik(:) = xij_(:,nbk)

         x_jk(:) = xij_(:,nbk) - x_ij(:)
         r2jk = x_jk(1)**2 + x_jk(2)**2 + x_jk(3)**2

         rjk = sqrt(r2jk)
         rik = rij_(nbk)
         rik1= 1.d0/rik

         xr_ik(:) = x_ik(:)*rik1

         if(r2jk.lt.Sqr_Rc15) then   ! rjk < 1.5Rc !

          Bb1Sik = Bb1Sij_(nbk) ! Bik*exp(-beta.rik)*Fc*Sik*bik !
          Pb1Sik = Pb1Sij_(nbk) ! sigma_i*(Psi(i)^-1/2)*Fc*Sik*bik !

          Bb1Ski = Bb1Sji_(nbk) ! Bki*exp(-beta.rik)*Fc*Ski*bki ! 
          Pb1Ski = Pb1Sji_(nbk) ! sigma_k*(Psi(k)^-1/2)*Fc*Ski*bki ! 

          rrr = rij + rjk - rik
           if(rrr.lt.Rc) then
          bl_ikj = U3(nlambda_bop)
          bl_kij = U3_of_ni(nlambda_bop,nbk)
          rrc = max(Rc-rrr,0.d0)
          rrc3 = rrc**3
          rrc4 = rrc*rrc3
          Fc_rrr = rrc4 / (d4+rrc4)   ! Fc(r_ikj) = Fc(r_kij) !
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

           else
          Bb1dSr_ikj = 0.d0
          Bb1dSr_kij = 0.d0
          Pb1dSr_ikj = 0.d0
          Pb1dSr_kij = 0.d0
           endif ! if(rrr.lt.Rc)... !

          rrr = rik + rjk - rij
           if(rrr.lt.Rc) then
          bl_ijk = U3(nlambda_bop)
          rrc = max(Rc-rrr,0.d0)
          rrc3 = rrc**3
          rrc4 = rrc*rrc3
          Fc_rrr = rrc4 / (d4+rrc4)
          Exp_blr_ijk = exp(-bl_ijk*rrr)
          Sijk_term = 1.d0 - Fc_rrr*Exp_blr_ijk  ! S_ijk !

          dFc = -4.d0*rrc3*d4 / (d4+rrc4)**2   ! dFc/dr !
          dS_ijk = (bl_ijk*Fc_rrr - dFc)*Exp_blr_ijk
          dS_ijk = dS_ijk/Sijk_term

          dSr_ij = dSr_ij + dS_ijk

           endif ! if(rrr.lt.Rc)... !

          dEsr_ij = dEsr_ij + Bb1dSr_ikj + Bb1dSr_kij
          dWsr_ij = dWsr_ij + Pb1dSr_ikj + Pb1dSr_kij

         endif ! if(r2jk.lt.Sqr_Rc15)... rjk < 1.5Rc !   

          aa_ikj = U3(naa_bop)
          h_ikj = U3(nh_bop)
          aa_ijk = U3(naa_bop)
          h_ijk = U3(nh_bop)

       cos_ikj = xr_ij(1)*xr_ik(1)+xr_ij(2)*xr_ik(2)+xr_ij(3)*xr_ik(3)
       cos_h_ikj = cos_ikj - h_ikj
       cos_h_ijk = cos_ikj - h_ijk

          xr_ik_ij(:) = (xr_ik(:) - xr_ij(:)*cos_ikj)*rij1
                      ! (xik/rik - xij/rij*cos_ikj)/rij !
          C1ikj = aa_ikj*cos_h_ikj
          C1ijk = aa_ijk*cos_h_ijk

          Bb3Sik =  Bb3Sij_(nbk)
          BexpFb3SikC1 =  Bb3Sik*C1ikj
        ! Bik*exp(-beta.rik)*Fc*Sik*(bik**3)*aa_ikj*cos_h_ikj !

          dE1bC(:) = dE1bC(:) + FcSik*C1ijk*xr_ik_ij(:)
        ! += Fc*Sik*aa_ijk*cos_h_ijk*(xik/rik - xij/rij*cos_ikj)/rij

        enddo ! do nbk = 1,nbrs_1st  ! rij < 1.5Rc; rik < Rc !

        F3_rij_BOP(nb1,ni) = dEsr_ij + dWsr_ij

       enddo ! do nb1 = nbrs_1st+1,nbrs_2nd  ! Rc < rij < 1.5Rc !

       dE_dparam_list(:,ni) = 0.5d0*dE_dparam(:)

C      write(100+mynod,*)'2d: dE_dparam_list(ni=',ni,')=',
C    1 dE_dparam_list(1,ni)

      enddo ! do ni=1,natoms !

!$OMP END PARALLEL DO
c
c ------------- Loop II: --------------
c
!$OMP PARALLEL DO DEFAULT(NONE)
!$OMP& PRIVATE(ni,nb1, nbrs_1st, nbrs_2nd, nj, i_at_j, nbk)
!$OMP& PRIVATE(xij,yij,zij,rij, x_ij, xij_,rij_)
!$OMP& PRIVATE(x_ik,rik, x_jk,rjk, sx0,sy0,sz0)
!$OMP& PRIVATE(bop_lambda_pair, bop_lambda_of, dE_dparam)
!$OMP& PRIVATE(FSBb3C2_ij, FSPb3C2_ij, FSBb3C2_ji, FSPb3C2_ji)
!$OMP& PRIVATE(dEbS1r_ij, dEbS2r_ij, dWbS1r_ij, dWbS2r_ij)
!$OMP& PRIVATE(rrr,rrc,rrc3,rrc4, bl_kij,bl_ikj, Fc_rrr, dFc)
!$OMP& PRIVATE(Exp_blr_kij, Exp_blr_ikj, Skij_term, Sikj_term)
!$OMP& PRIVATE(dS_kij, dS_ikj, FSBb3C2_ik, FSPb3C2_ik)
!$OMP& SHARED(natoms, nbr_list)
!$OMP& SHARED(sx,sy,sz, h11,h22,h33,h12,h13,h23, Rc, d4)
!$OMP& SHARED(nbrs_all_of, nbrs_2nd_of, nbr_list_inv)
!$OMP& SHARED(FSBb3C2_ij_list, FSPb3C2_ij_list, F3_rij_BOP, U3_of)
!$OMP& SCHEDULE(DYNAMIC,CHUNK)

      do ni=1,natoms

       nbrs_1st = nbrs_all_of(ni)
       nbrs_2nd = nbrs_2nd_of(ni)  ! rij < 1.5Rc !

       ! get BOP parameters for the current atom i !

       bop_lambda_pair = U3_of(nlambda_bop,ni)

       do nb1 = 1,nbrs_2nd     ! Loop 0: i - j bond; rij < 1.5Rc !
        nj = nbr_list(nb1,ni)       ! VECTORIZED: speedup: 1.620 !

        sx0 = sx(nj) - sx(ni)
        sy0 = sy(nj) - sy(ni)
        sz0 = sz(nj) - sz(ni)
        sx0 = sx0 - dnint(sx0)     ! make periodic along X !
        sy0 = sy0 - dnint(sy0)     ! make periodic along Y !
        sz0 = sz0 - dnint(sz0)     ! make periodic along Z !

        xij = h11*sx0 + h12*sy0 + h13*sz0
        yij = h22*sy0 + h23*sz0
        zij = h33*sz0

        rij = sqrt(xij**2 + yij**2 + zij**2)
        rij_(nb1) = rij
        xij_(1,nb1) = xij
        xij_(2,nb1) = yij
        xij_(3,nb1) = zij
       enddo ! do nb1 = 1,nbrs_2nd !

       do nb1 = 1,nbrs_1st
        nj = nbr_list(nb1,ni)       ! VECTORIZED: speedup: 1.620 !
        i_at_j = nbr_list_inv(nb1,ni)  ! get i as a nbr of j !
        FSBb3C2_ij(nb1) = FSBb3C2_ij_list(nb1,ni)
        FSPb3C2_ij(nb1) = FSPb3C2_ij_list(nb1,ni)
        FSBb3C2_ji(nb1) = FSBb3C2_ij_list(i_at_j,nj)
        FSPb3C2_ji(nb1) = FSPb3C2_ij_list(i_at_j,nj)
        bop_lambda_of(nb1) = U3_of(nlambda_bop,nj)
       enddo ! do nb1 = 1,nbrs_1st !

       do nb1 = 1,nbrs_2nd     ! Loop 0: i - j bond; rij < 1.5Rc !
        rij = rij_(nb1)
        x_ij(:) = xij_(:,nb1)

        dEbS1r_ij = 0.d0; dEbS2r_ij = 0.d0
        dWbS1r_ij = 0.d0; dWbS2r_ij = 0.d0

        do nbk = 1,nbrs_1st     ! rij < 1.5Rc; rik < Rc !
          if(nbk.ne.nb1) then
         rik = rij_(nbk)
         x_ik(:) = xij_(:,nbk)
         x_jk(:) = x_ik(:) - x_ij(:)
         rjk = sqrt(x_jk(1)**2 + x_jk(2)**2 + x_jk(3)**2)

          ! Calc. dSr_kij !

          rrr = rjk + rij - rik
           if(rrr.lt.Rc) then
          bl_kij = bop_lambda_of(nbk)
          bl_ikj = bop_lambda_pair
          rrc = max(Rc-rrr,0.d0)
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

          dEbS1r_ij = dEbS1r_ij + dS_kij*FSBb3C2_ji(nbk) 
                              ! + dS_kij*Fc_ik*Ski*Bb3C2_ki !
          dWbS1r_ij = dWbS1r_ij + dS_kij*FSPb3C2_ji(nbk) 
                              ! + dS_kij*Fc_ik*Ski*Pb3C2_ki !
         
          FSBb3C2_ik = FSBb3C2_ij(nbk)
          dEbS2r_ij = dEbS2r_ij + dS_ikj*FSBb3C2_ik
                              ! + dS_ikj*Fc_ik*Sik*Bb3C2_ik !
          FSPb3C2_ik = FSPb3C2_ij(nbk)
          dWbS2r_ij = dWbS2r_ij + dS_ikj*FSPb3C2_ik
                              ! + dS_ikj*Fc_ik*Sik*Pb3C2_ik !
           endif ! if(rrr.lt.Rc)... !
         
          endif ! if(nbk.ne.nb1)... ! 
        enddo ! do nbk = 1,nbrs_1st     ! rij < 1.5Rc; rik < Rc !

        ! Add the 1st term of dEbSr_ij and dWbSr_ij !
        F3_rij_BOP(nb1,ni) = F3_rij_BOP(nb1,ni) - 
     -  0.5d0*(dEbS1r_ij + dEbS2r_ij + dWbS1r_ij + dWbS2r_ij)

       enddo ! do nb1 = 1,nbrs_2nd !
      enddo ! do ni=1,natoms !

!$OMP END PARALLEL DO

      return
      end subroutine      ! dBOP_OMP !
!
! ---------------------------------------------------------------------
! Calculates Analytical derivatives and force calculation.
! ---------------------------------------------------------------------
!
      subroutine Frc_PINN_OMP(ecohe)

      use atoms

      implicit double precision (a-h,o-z)
  
      double precision, intent(out) :: ecohe

      double precision, dimension(3) :: x_ij,xr_ij,ff,fr,Fa_ij
      real(kind=4) Rc15

      ecohe = 0.d0

      if(ihalt.ne.0) return

      Rc = Rc_ann
      Rc15 = 1.5d0*Rc

      h11 = h(1,1); h12 = h(1,2); h13 = h(1,3)
      h22 = h(2,2); h23 = h(2,3); h33 = h(3,3)

      call dANN_OMP(max2_nbrs, Rc15) ! aladyn_pi_PINN_OMP, get dANN_dxij
      call dBOP_OMP(ecohe)           ! Get BOP function derivatives !

! --- Calc Actual Force Vectors and Energy ---

!$OMP PARALLEL DO DEFAULT(NONE)
!$OMP& PRIVATE(ni,nb1, nbrs_1st, nbrs_2nd, nj, nbrs_of_j, i_at_j)
!$OMP& PRIVATE(x_ij,rij,rij1, xr_ij, sx0,sy0,sz0)
!$OMP& PRIVATE(Sr_ij,Sr_ji, ff,fr, Fa_ij, Fw_ij,Fw_ji, Fr_ij)
!$OMP& PRIVATE(bop_sigma_pair,bop_sigma_of, Sigma_Psi_i,Sigma_Psi_j)
!$OMP& SHARED(natoms, iatom_types, Nodes_of_layer, net_layers)
!$OMP& SHARED(sx,sy,sz, h11,h22,h33,h12,h13,h23)
!$OMP& SHARED(ntype, nbr_list, nbrs_all_of,nbrs_2nd_of, nbr_list_inv)
!$OMP& SHARED(F2_rij_BOP, F3_rij_BOP, FCxij_BOP, PsidFij_list)
!$OMP& SHARED(dE_dparam_list, dANN_dxij)
!$OMP& SHARED(Psi_of, U3_of, frr)
!$OMP& SCHEDULE(DYNAMIC,CHUNK)

      do ni=1,natoms !  ! Loop III: Full forces !

       ! get BOP parameters for the current atom i !
       bop_sigma_pair = U3_of(nsigma_bop,ni)

       nbrs_1st = nbrs_all_of(ni)
       nbrs_2nd = nbrs_2nd_of(ni)

       fr(1:3) = 0.d0

       if(nbrs_2nd.gt.0) then
        Sigma_Psi_i = bop_sigma_pair/sqrt(Psi_of(ni))
                  ! = sigma_i*(Psi(i)^-1/2) !
       else
        Sigma_Psi_i = 0.d0
       endif

       do nb1 = 1,nbrs_2nd     ! Loop Ia: i - j bond !
        nj = nbr_list(nb1,ni)
        nbrs_of_j = nbrs_2nd_of(nj)  ! rij < 1.5Rc !

        sx0 = sx(nj) - sx(ni)
        sy0 = sy(nj) - sy(ni)
        sz0 = sz(nj) - sz(ni)
        sx0 = sx0 - dnint(sx0)     ! make periodic along X !
        sy0 = sy0 - dnint(sy0)     ! make periodic along Y !
        sz0 = sz0 - dnint(sz0)     ! make periodic along Z !

        x_ij(1) = h11*sx0 + h12*sy0 + h13*sz0
        x_ij(2) = h22*sy0 + h23*sz0
        x_ij(3) = h33*sz0
        rij = sqrt(x_ij(1)**2 + x_ij(2)**2 + x_ij(3)**2)
        rij1 = 1.d0/rij

        xr_ij(:) = x_ij(:)*rij1  ! xij/rij !

        i_at_j = nbr_list_inv(nb1,ni)

        Sr_ij = F3_rij_BOP(nb1,ni)
        Sr_ji = F3_rij_BOP(i_at_j,nj)

        bop_sigma_of(nb1) = U3_of(nsigma_bop,nj)

        if(nbrs_of_j.gt.0) then
         Sigma_Psi_j = bop_sigma_of(nb1)/sqrt(Psi_of(nj))
                   ! = sigma_j*(Psi(j)^-1/2) !
        else
         Sigma_Psi_j = 0.d0
        endif

        ! Initialize with the first BOP parameter !
        Fa_ij(:) = 
     =  dE_dparam_list(1,ni)*dANN_dxij(:,1,nb1,ni) - 
     -  dE_dparam_list(1,nj)*dANN_dxij(:,1,i_at_j,nj)

        ! Loop of parameter count ip= 2,.. 8 !
        do ip = 2, Nodes_of_layer(net_layers)
         Fa_ij(:) = Fa_ij(:) +
     =   dE_dparam_list(ip,ni)*dANN_dxij(:,ip,nb1,ni) - 
     -   dE_dparam_list(ip,nj)*dANN_dxij(:,ip,i_at_j,nj)
        enddo

        if(nb1.le.nbrs_1st) then
         Fw_ij = Sigma_Psi_i*PsidFij_list(nb1,ni)
            ! = sigma_i*Sij*bij*dFc*(Psi(i)^-1/2) !
         Fw_ji = Sigma_Psi_j*PsidFij_list(i_at_j,nj)
            ! = sigma_j*Sji*bji*dFc*(Psi(j)^-1/2) !
         Fr_ij = 0.5d0*(F2_rij_BOP(nb1,ni) + F2_rij_BOP(i_at_j,nj) - 
     -                 (Fw_ij + Fw_ji))
         ff(:) = (Fr_ij - 0.5d0*(Sr_ij+Sr_ji))*xr_ij(:) -
     -   0.5d0*(FCxij_BOP(:,nb1,ni)-FCxij_BOP(:,i_at_j,nj)) + Fa_ij(:)

        else
         ff(:) = -0.5d0*(Sr_ij+Sr_ji)*xr_ij(:) + Fa_ij(:)
        endif

        fr(1:3) = fr(1:3) + ff(1:3)   ! x,y,z forces !

       enddo ! do nb1 = 1,nbrs_2nd !

       frr(1:3,ni) = fr(1:3)

      enddo ! do ni=1,natoms !  ! Loop III: no stress !

!$OMP END PARALLEL DO

C     write(6,*)'Frc_PINN_OMP: ecohe=',ecohe

C     write(300+mynod,*)' istress=',istress
C     do n=1,natoms
C       if(ident(n).eq.1) then
C      write(300+mynod,20) ident(n),0.5d0*Ep_of(n),frr(1:4,n)
C       endif
C     enddo
C 20  format('frr: ',i6,' ',5f16.12)

      return
      end subroutine      ! Frc_PINN_OMP  !
c
c ---------------------------------------------------------------------
c
      subroutine alloc_atoms_PINN(ierror)

      integer, intent(out) :: ierror
      integer :: ialloc(16)

      ialloc(:) = 0

      allocate(id_nbr(nbrs_per_atom), stat=ialloc(1))
      allocate(bop_lambda_of(nbrs_per_atom), stat=ialloc(2))
      allocate(bop_sigma_of(nbrs_per_atom), stat=ialloc(3))

      allocate(VB_of(natoms_alloc), stat=ialloc(4))
      allocate(Psi_of(natoms_alloc), stat=ialloc(5))
      allocate(Psi_j(nbrs_per_atom), stat=ialloc(6))

      allocate(nbrs_all_of(natoms_alloc), stat=ialloc(7))
      allocate(nbrs_2nd_of(natoms_alloc), stat=ialloc(8))

      allocate(nbr_list_inv(nbrs_per_atom,natoms_alloc), 
     1 stat=ialloc(9))

      allocate(Sijn(nbrs1_per_atom), stat=ialloc(1))
      allocate(Zijn(nbrs1_per_atom), stat=ialloc(3))
      allocate(U3_of_ni(nBOP_params, nbrs_per_atom), stat=ialloc(6))

      allocate(BSij_list(nbrs1_per_atom,natoms_alloc), stat=ialloc(2))
      allocate(BZij_list(nbrs1_per_atom,natoms_alloc), stat=ialloc(3))
      allocate(dBZij_list(nbrs1_per_atom,natoms_alloc), stat=ialloc(9))
      allocate(dE7Sij_list(nbrs1_per_atom,natoms_alloc),stat=ialloc(10))
      allocate(dE_dparam_list(nBOP_params,natoms_alloc),stat=ialloc(16))

       if(I_have_GPU.eq.0) then ! OMP only arrays !
      allocate(F2_rij_BOP(nbrs1_per_atom,natoms_alloc),stat=ialloc(10))
      allocate(F3_rij_BOP(nbrs_per_atom,natoms_alloc), stat=ialloc(11))
      allocate(FCxij_BOP(3,nbrs1_per_atom,natoms_alloc),stat=ialloc(12))
      allocate(FSBb3C2_ij_list(nbrs1_per_atom,natoms_alloc),
     1 stat=ialloc(13))
      allocate(FSPb3C2_ij_list(nbrs1_per_atom,natoms_alloc),
     1 stat=ialloc(14))
      allocate(PsidFij_list(nbrs1_per_atom,natoms_alloc), 
     1 stat=ialloc(15))
      allocate(FcSij_ni(nbrs1_per_atom), stat=ialloc(16))
       endif

      ierror = 0
      do i=1,16
       ierror = ierror + ialloc(i)
      enddo

      return
      end subroutine        ! alloc_atoms_PINN !
c
c ---------------------------------------------------------------------
c
      subroutine deall_atoms_PINN(ierror)

      integer, intent(out) :: ierror
      integer ialloc(16)

      ialloc(:) = 0

      if (allocated(id_nbr)) deallocate(id_nbr,stat=ialloc(1))
      if (allocated(bop_lambda_of)) deallocate(bop_lambda_of,
     & stat=ialloc(2))
      if(allocated(bop_sigma_of))deallocate(bop_sigma_of,stat=ialloc(3))

      if (allocated(VB_of)) deallocate(VB_of,stat=ialloc(4))
      if (allocated(Psi_of)) deallocate(Psi_of,stat=ialloc(5))
      if (allocated(Psi_j)) deallocate(Psi_j,stat=ialloc(6))

      if (allocated(nbrs_all_of)) deallocate(nbrs_all_of,stat=ialloc(7))
      if (allocated(nbrs_2nd_of)) deallocate(nbrs_2nd_of,stat=ialloc(8))

      if (allocated(Sijn)) deallocate(Sijn, stat=ialloc(2))
      if (allocated(Zijn)) deallocate(Zijn, stat=ialloc(2))
      if(allocated(nbr_list_inv))deallocate(nbr_list_inv,stat=ialloc(9))
      if (allocated(F2_rij_BOP)) deallocate(F2_rij_BOP, stat=ialloc(10))
      if (allocated(FCxij_BOP)) deallocate(FCxij_BOP, stat=ialloc(11))
      if (allocated(FcSij_ni)) deallocate(FcSij_ni, stat=ialloc(11))
      if (allocated(U3_of_ni)) deallocate(U3_of_ni, stat=ialloc(6))

      if (allocated(BSij_list)) deallocate(BSij_list, stat=ialloc(2))
      if (allocated(BZij_list)) deallocate(BZij_list, stat=ialloc(3))
      if (allocated(dBZij_list)) deallocate(dBZij_list, stat=ialloc(9))
      if (allocated(dE7Sij_list))deallocate(dE7Sij_list,stat=ialloc(10))

      if (allocated(FSBb3C2_ij_list)) deallocate(FSBb3C2_ij_list, 
     & stat=ialloc(12))
      if (allocated(FSPb3C2_ij_list)) deallocate(FSPb3C2_ij_list, 
     & stat=ialloc(13))

      if (allocated(F3_rij_BOP)) deallocate(F3_rij_BOP, stat=ialloc(14))
      if (allocated(PsidFij_list)) deallocate(PsidFij_list,
     & stat=ialloc(15))
      if (allocated(dE_dparam_list)) deallocate(dE_dparam_list,
     & stat=ialloc(16))

      ierror = 0
      do i=1,16
       ierror = ierror + ialloc(i)
      enddo

      return
      end subroutine        ! deall_atoms_PINN !
!
! ------------------------------------------------------------------
!
      END MODULE  ! PINN !
