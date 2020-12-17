c------------------------------------------------------------------
c 11-20-2020
c
c Miolecular Dynamics Module Unit for aladyn_pi.f code
c
c Vesselin Yamakov
c National Institute of Aerospace
c 100 Exploration Way,
c Hampton, VA 23666 
c phone: (757)-864-2850
c fax:   (757)-864-8911
c e-mail: yamakov@nianet.org
c------------------------------------------------------------------
c
      MODULE MD   ! NNNN !

      use constants
      use sim_box
      use pot_module
      use atoms
      save

      ! Predictor-corrector coefficients !

      double precision :: f02 = 3.0d0/20.0d0
      double precision :: f02viscous = 3.0d0/16.0d0
      double precision :: f12 = 251.0d0/360.0d0
      double precision :: f32 = 11.0d0/18.0d0
      double precision :: f42 = 1.0d0/6.0d0
      double precision :: f52 = 1.0d0/60.0d0
      double precision :: f02_wall = 3.0d0/20.0d0
      double precision :: f02_atom = 3.0d0/20.0d0


      CONTAINS
c
c--------------------------------------------------------------------
c
      subroutine init_MD

      implicit double precision (a-h,o-z)

      ! Fifth order predictor-corrector coefficients !
       f02 = 3.0d0/20.0d0        ! in sim_box module !
       f02viscous = 3.0d0/16.0d0
       f12 = 251.0d0/360.0d0
       ! f22 = 1.0d0 !
       f32 = 11.0d0/18.0d0
       f42 = 1.0d0/6.0d0
       f52 = 1.0d0/60.0d0
       f02_wall = f02
       f02_atom = f02

      return
      end subroutine         ! init_MD !
c
c-----------------------------------------------------------------------
c Called either from read_structure with or without velocities,
c or from SIM_run when MD starts
c-----------------------------------------------------------------------
c
      subroutine init_vel(T_set0)
      implicit double precision (a-h,o-z)
c
c   ***  assign initial velocities to atoms  ***************************
c
      hi11 = hi(1,1); hi12 = hi(1,2); hi13 = hi(1,3)
      hi22 = hi(2,2); hi23 = hi(2,3); hi33 = hi(3,3)

      tfc=2.0/(3.0*Boltz_Kb)
c
c   ***  use ranmar number gen. of Marsaglia and Zaman FSU-SCRI-87-50
c   ***  irr must be integer and originally set to an odd large number
c   ***  xx,yy,zz range between 0.0 and 1.0
c
      irr=6751  ! VVVV for test with escm_init only !
      call rmarin(irr+17*mynod) ! node depndnt init. random generator !

      do i=1,iatom_types
       sumPx(i)=0.0d0
       sumPy(i)=0.0d0
       sumPz(i)=0.0d0
       E_kin(i)=0.0d0
      enddo

      sumx=0.0
      sumy=0.0
      sumz=0.0
      call ranmar(3*natoms)
      ii=1 
      do i=1,natoms
       ntp = ntype(i)
       Am = amass(ntp)
       xx=buf(ii)
       yy=buf(ii+1)  ! get random velocity numbers !
       zz=buf(ii+2)
       ii=ii+3 
       xyz= sqrt(xx*xx+yy*yy+zz*zz)
       x1(i)=xx/xyz
       y1(i)=yy/xyz
       z1(i)=zz/xyz
       sumPx(ntp) = sumPx(ntp) + Am*x1(i)
       sumPy(ntp) = sumPy(ntp) + Am*y1(i) ! in Real space !
       sumPz(ntp) = sumPz(ntp) + Am*z1(i)
      enddo

      tot_Px = 0.d0; tot_Py = 0.d0; tot_Pz = 0.d0; tot_M = 0.d0
      do ntp=1,iatom_types
       Am = amass(ntp)
       na = max(1,natoms_of_type(ntp))
       tot_Px = tot_Px + sumPx(ntp)
       tot_Py = tot_Py + sumPy(ntp)
       tot_Pz = tot_Pz + sumPz(ntp)
       tot_M  = tot_M  + Am*na
       sumPx(ntp)=0.d0; sumPy(ntp)=0.d0; sumPz(ntp)=0.d0
      enddo
c
c    ***  adjust velocities such that total linear momentum is zero.
c
      do i=1,natoms   ! V_cm = Sum_i{mi*Vi}/Sum_i{mi} !
       ntp = ntype(i)
       Am = amass(ntp)
       x1(i) = x1(i) - tot_Px/tot_M   ! vi = vi - V_cm !
       y1(i) = y1(i) - tot_Py/tot_M
       z1(i) = z1(i) - tot_Pz/tot_M
       vx = x1(i)
       vy = y1(i)    ! R - space [A/ps] !
       vz = z1(i)
       sumPx(ntp) = sumPx(ntp) + Am*x1(i) ! New momentums !
       sumPy(ntp) = sumPy(ntp) + Am*y1(i) ! in Real space !
       sumPz(ntp) = sumPz(ntp) + Am*z1(i)
       E_kin(ntp) = E_kin(ntp) + Am*(vx*vx + vy*vy + vz*vz)

C     if(ident(i).eq.1) then
C      write(50,110)i,ident(i),x1(i),y1(i),z1(i),sumPx(ntp),
C    1 sumPy(ntp),sumPz(ntp),Am,E_kin(ntp)
C     endif

      enddo

      pTemp(0) = 0.d0
      do i=1,iatom_types
       na = natoms_of_type(i)
       if(na.gt.0) then 
        pTemp(i) = 0.5d0*tfc*E_kin(i)/na  ! = 2/3*Ek/Kb !
        if(pTemp(i).gt.0.001d0) then
         Tscale(i) = sqrt(T_set0/pTemp(i))
        else
         Tscale(i) = 0.0d0  ! if T is too low, freeze this atom !
        endif
       else  ! when na=0 !
        pTemp(i) = 0.0d0  ! No atom is rescaled to V=0, E_kin=0 !
        Tscale(i) = 0.0d0
        na = 1
       endif ! if(na.gt.0)... !

C     write(50,*)'init_vel: natoms_of_type(',i,')=',
C    1 natoms_of_type(i),' T=',pTemp(i),' Ek=',E_kin(i)/na,' atu=',atu,
C    1' T_0=',T_set0,' Tscale=',Tscale(i),' tfc=',tfc,' mass=',amass(i),
C    1' na=',na

      enddo ! do i=1,iatom_types !

c
c  ***   scale velocities to desired temp. ***********************
c
       do i=1,natoms
        ntp = ntype(i)
        x1(i)=x1(i)*Tscale(ntp) ! R - space !
        y1(i)=y1(i)*Tscale(ntp) 
        z1(i)=z1(i)*Tscale(ntp)
       enddo
 
      ! Convert to S-space !

      do i=1,natoms

C      if(ident(i).eq.1) then
C       write(50,15)ident(i),x1(i),y1(i),z1(i)
C      endif
C 15   format('init_vel1: id:',i5,' R x1=',3E18.10)

       x3(i) = (hi11*x1(i) + hi12*y1(i) + hi13*z1(i))
       y3(i) = (hi22*y1(i) + hi23*z1(i))
       z3(i) = hi33*z1(i)

       x1(i) = x3(i)*dt_step
       y1(i) = y3(i)*dt_step
       z1(i) = z3(i)*dt_step

       x2(i) = 0.d0; y2(i) = 0.d0; z2(i) = 0.d0
       x3(i) = 0.d0; y3(i) = 0.d0; z3(i) = 0.d0
       x4(i) = 0.d0; y4(i) = 0.d0; z4(i) = 0.d0
       x5(i) = 0.d0; y5(i) = 0.d0; z5(i) = 0.d0

C      if(ident(i).eq.1) then
C       write(50,16)ident(i),x1(i),y1(i),z1(i)
C      endif
C 16   format('init_vel2: id:',i5,' S x1*dt=',3E18.10)

      enddo

      do i=1,3
      do j=1,3
       h1(i,j)=0.d0; h2(i,j)=0.d0; h3(i,j)=0.d0; h4(i,j)=0.d0;
       h5(i,j)=0.d0
      enddo
      enddo
     
      A_fr = 0.d0 ! Reset dissipated friction energy !
      Q_heat = 0.d0

      call get_T ! calc. pTemp(ntp), E_kin(ntp), T_sys and Ek_sys !

      return
      end subroutine       ! init_vel !
c
c-----------------------------------------------------------------------
c  Sets second derivatives (accelerations) according to forces,   
c  so that the predictor-corrector can start with supplied 
c  first and second derivatives.
c-----------------------------------------------------------------------
c
      subroutine initaccel

      implicit double precision (a-h,o-z)

      dtsqh = 0.5d0*dt_step**2

      hi11 = hi(1,1); hi12 = hi(1,2); hi13 = hi(1,3)
      hi22 = hi(2,2); hi23 = hi(2,3); hi33 = hi(3,3)

!$OMP PARALLEL DO PRIVATE(n,ntp,fxsn,fysn,fzsn,dtsqhM)
!$OMP& SCHEDULE(DYNAMIC,CHUNK)

      do n=1,natoms
       ntp = ntype(n)
       dtsqhM = dtsqh/amass(ntp)  ! for predict - correct !
!      fxsn = hi11*frr1(n) + hi12*frr2(n) + hi13*frr3(n)
!      fysn = hi22*frr2(n) + hi23*frr3(n)  ! S-space !
!      fzsn = hi33*frr3(n)

       fxsn = hi11*frr(1,n) + hi12*frr(2,n) + hi13*frr(3,n)
       fysn = hi22*frr(2,n) + hi23*frr(3,n)  ! S-space !
       fzsn = hi33*frr(3,n)

       x2(n) = dtsqhM*fxsn      ! Taylor expansion in s-space : !
       y2(n) = dtsqhM*fysn      ! vel = dt * dr/dt  !
       z2(n) = dtsqhM*fzsn
      enddo                       ! acc = 1/2*dt^2 * f/m !

!$OMP END PARALLEL DO

      return
      end subroutine  ! initaccel !
c
c ---------------------------------------------------------------------
c Uses x1(n), y1(n), and z1(n) in S-space only to calculate
c sumPxyz(n), pTemp(ntp), E_kin(ntp), T_sys and Ek_sys
c ---------------------------------------------------------------------
c
      subroutine get_T

      implicit double precision (a-h,o-z)

      dt = dt_step
      tfc=2.0/(3.0*Boltz_Kb)
    
      do i=1,iatom_types
       Am_of_type(i) = 0.0d0
       sumPx(i)=0.0d0
       sumPy(i)=0.0d0
       sumPz(i)=0.0d0
       E_kin(i)=0.0d0
      enddo

!$OMP PARALLEL DO DEFAULT(SHARED)
!$OMP& PRIVATE(ntp,Am,Vx,Vy,Vz)
!$OMP& REDUCTION(+:Am_of_type,sumPx,sumPy,sumPz,E_kin)
!$OMP& SCHEDULE(DYNAMIC,CHUNK)

      do i=1,natoms
       ntp = ntype(i)
       Am = amass(ntp)
       Am_of_type(ntp) = Am_of_type(ntp) + Am
       Vx = h(1,1)*x1(i) + h(1,2)*y1(i) + h(1,3)*z1(i)
       Vy = h(2,2)*y1(i) + h(2,3)*z1(i)  ! Real units !
       Vz = h(3,3)*z1(i)
       sumPx(ntp) = sumPx(ntp) + Am*Vx
       sumPy(ntp) = sumPy(ntp) + Am*Vy
       sumPz(ntp) = sumPz(ntp) + Am*Vz
       E_kin(ntp) = E_kin(ntp) + Am*(Vx**2 + Vy**2 + Vz**2) ! mv^2 !
      enddo

!$OMP END PARALLEL DO

       Ek_sys = 0.d0
       avr_mass = 0.0d0
       do k=1,iatom_types
        na = natoms_of_type(k)  ! do not use max(1,natoms_of_type(k)) !
        if(na.gt.0) then
         sumPx(k) = sumPx(k)/na
         sumPy(k) = sumPy(k)/na
         sumPz(k) = sumPz(k)/na
         Ek_sys = Ek_sys + E_kin(k)  ! Sum_(mv^2) !
         E_kin(k) = E_kin(k)/na      ! mv^2 !
         sum_mass(k) = Am_of_type(k)
         avr_mass = avr_mass + sum_mass(k)
        else
         sumPx(k) = 0.d0; sumPy(k) = 0.d0; sumPz(k) = 0.d0
         E_kin(k) = 0.d0; sum_mass(k) = 0.d0
        endif
       enddo

       Ek_sys = Ek_sys/(2.0d0*natoms*dt**2)  ! mv^2 / 2 !
       avr_mass = avr_mass/natoms

        T_sys = 0.0d0
        Ttemp = 1.0d0/(3.0d0*Boltz_Kb*dt**2)
       do ntp=1,iatom_types
        if(natoms_of_type(ntp).gt.1) then  ! more than 1 atom of type !
         pTemp(ntp) = Ttemp*(E_kin(ntp) -
     -   (sumPx(ntp)**2 + sumPy(ntp)**2 + sumPz(ntp)**2)/sum_mass(ntp))
        else if(natoms_of_type(ntp).eq.1) then  ! only 1 atom of type !
         pTemp(ntp) = Ttemp*E_kin(ntp)     ! (2/3kT)*E_k = mv^2 / 3kT !
        else                               ! No atoms of type ntp     !
         pTemp(ntp) = 0.d0
        endif
        T_sys = T_sys + pTemp(ntp)*sum_mass(ntp)/amass(ntp)/natoms

C      write(1000+mynod,10)ntp,natoms_of_type(ntp),pTemp(ntp),T_set,
C    1 amass(ntp)/atu,sum_mass(ntp),T_sys,Ek_sys
C 10   format('get_T: natoms_of_type(',i2,')=',i6,' T=',f10.2,
C    1 ' T_set=',f10.2,' mass=',f10.3,' sum_mass=',f10.3,
C    1 ' T_sys=',f10.2,' Ek_sys=',f12.8)

        enddo  ! do ntp=1,iatom_types !

      return
      end subroutine       ! get_T !
c
c-----------------------------------------------------------------------
c Collects Px,y,z and Ek from all nodes after correct_(slow,fast)_atoms
c-----------------------------------------------------------------------
c
      subroutine T_broadcast

      implicit double precision (a-h,o-z)

       dt = dt_step     !  uses the basic time step !
       dtsq = dt**2
       Ttemp = 1.0d0/(3.0d0*Boltz_Kb*dtsq)  ! 1 / (3kB.T) !

       Q_heat = Q_heat + A_fr/max(1,natoms)
       A_fr = 0.0d0

       do ntp=1,iatom_types
       na = natoms_of_type(ntp)  ! do not use max(1,natoms_of_type()) !
        if(na.gt.0) then
         sumPx(ntp) = sumPx(ntp)/na
         sumPy(ntp) = sumPy(ntp)/na
         sumPz(ntp) = sumPz(ntp)/na
         E_kin(ntp) = E_kin(ntp)/na   ! mv^2 !
         if(na.gt.1) then
          pTemp(ntp) = Ttemp*(E_kin(ntp) - 
     -    (sumPx(ntp)**2 + sumPy(ntp)**2 + sumPz(ntp)**2)/sum_mass(ntp))
         else
          pTemp(ntp) = Ttemp*E_kin(ntp)    ! (2/3kT)*E_k = mv^2 / 3kT !
         endif
        endif       ! 2Ek - P^2/M = 2(Ek - P^2/2M) !
       enddo

c      ! T = (2Ek - P^2/M)/3Kb = 2(Ek - P^2/2M)/3Kb !
c      ! Ek, P, and M are the TOTAL kin. energy, momentum, and mass !

C      do ntp=1,iatom_types
C       write(1000+mynod,10)ntp,pTemp(ntp),sum_mass(ntp),
C    1  E_kin(ntp)/(2.0d0*dt**2), natoms_of_type(ntp)
C      enddo
C 10  format('TB: pTemp(',i2,')=',f7.2,' sum_mass=',f12.8,
C    1 ' E_kin=',f12.8,' na_of_type=',i5)

      return
      end subroutine     ! T_broadcast !
c
c-----------------------------------------------------------------------
c predicts pos., vel., and higher deriv. 
c Fast atoms and their neighbors are using the smallest basic dt,
c Slow atoms that do not neighbor a fast atom use the long time step.
c-----------------------------------------------------------------------
c
      subroutine predict_atoms(ndof_fl)

      implicit double precision (a-h,o-z)

      integer, intent(in) :: ndof_fl

      double precision :: a1i,a2i,a3i,a4i,a5i,a24,a45,a2345
      double precision, dimension(3) :: dsx
      integer :: i,kdof

c   ***  start predictor step of integration scheme  ******************
c   ***  all x,y,z coordinates are s-scaled coordinates
c        and are scale invariant (do not change with length units).
c
      h11 = h(1,1); h12 = h(1,2); h13 = h(1,3)
      h22 = h(2,2); h23 = h(2,3); h33 = h(3,3)
      hi11 = hi(1,1); hi12 = hi(1,2); hi13 = hi(1,3)
      hi22 = hi(2,2); hi23 = hi(2,3); hi33 = hi(3,3)

!$OMP PARALLEL DO PRIVATE(i,a1i,a2i,a3i,a4i,a5i,a24,a45,a2345)
!$OMP& SCHEDULE(DYNAMIC,CHUNK)

      do i=1,natoms
       a1i = x1(i); a2i = x2(i); a3i = x3(i); a4i = x4(i); a5i = x5(i)
       sx(i) = sx(i) + a1i + a2i + a3i + a4i + a5i
       a24 = 2.d0*a4i
       a45 = a24 + 5.0d0*a5i
       a2345 = a2i + 3.0d0*a3i + a45 + a24
       x1(i) = a1i + a2i + a2345
       x2(i) = a2345 + a45
       x3(i) = a3i + 2.0d0*a45
       x4(i) = a45 - a4i

       a1i = y1(i); a2i = y2(i); a3i = y3(i); a4i = y4(i); a5i = y5(i)
       sy(i) = sy(i) + a1i + a2i + a3i + a4i + a5i
       a24 = 2.d0*a4i
       a45 = a24 + 5.0d0*a5i
       a2345 = a2i + 3.0d0*a3i + a45 + a24
       y1(i) = a1i + a2i + a2345
       y2(i) = a2345 + a45
       y3(i) = a3i + 2.0d0*a45
       y4(i) = a45 - a4i

       a1i = z1(i); a2i = z2(i); a3i = z3(i); a4i = z4(i); a5i = z5(i)
       sz(i) = sz(i) + a1i + a2i + a3i + a4i + a5i
       a24 = 2.d0*a4i
       a45 = a24 + 5.0d0*a5i
       a2345 = a2i + 3.0d0*a3i + a45 + a24
       z1(i) = a1i + a2i + a2345
       z2(i) = a2345 + a45
       z3(i) = a3i + 2.0d0*a45
       z4(i) = a45 - a4i
      enddo ! do i=1,natoms !

!$OMP END PARALLEL DO

      return
      end subroutine      ! predict_atoms !
c
c-----------------------------------------------------------------------
c *** this is the SLOW atoms corrector step ***
c-----------------------------------------------------------------------
c
      subroutine correct_atoms(ndof_fl)

      implicit double precision (a-h,o-z)

      integer, intent(in) :: ndof_fl

      integer :: n,ntp
      double precision :: dt,dtsqh,dtsqhM,sai
      double precision :: am,fr_ntp
      double precision :: x1i,xmp,afrx,Vxt,fx
      double precision :: y1i,ymp,afry,Vyt,fy
      double precision :: z1i,zmp,afrz,Vzt,fz
      double precision, dimension(3) :: dsx,sum_dcm0

       dt = dt_step
       dtsqh = 0.5d0*(dt**2)

       do n=1,iatom_types
        sumPx(n)=0.0d0
        sumPy(n)=0.0d0
        sumPz(n)=0.0d0
        E_kin(n)=0.0d0
       enddo

       h11 = h(1,1); h12 = h(1,2); h13 = h(1,3)
       h22 = h(2,2); h23 = h(2,3); h33 = h(3,3)

       hi11 = hi(1,1); hi12 = hi(1,2); hi13 = hi(1,3)
       hi22 = hi(2,2); hi23 = hi(2,3); hi33 = hi(3,3)

!$OMP PARALLEL DO PRIVATE(n,ntp,kdof,am,dtsqhM,sai)
!$OMP& PRIVATE(x1i,xmp,Vxt,fx,fxsn)
!$OMP& PRIVATE(y1i,ymp,Vyt,fy,fysn)
!$OMP& PRIVATE(z1i,zmp,Vzt,fz,fzsn)
!$OMP& REDUCTION(+: sumPx, sumPy, sumPz, E_kin)
!$OMP& SCHEDULE(DYNAMIC,CHUNK)

      do n=1,natoms

       ntp = ntype(n)
       am = amass(ntp)
       dtsqhM = dtsqh/am  ! = 0.5*dt*dt/m with dt=dt_step !
                          ! a = - 1/2*dt^2 * f/m; dt=dt_step !

!      fxsn = hi11*frr1(n) + hi12*frr2(n) + hi13*frr3(n)
!      fysn = hi22*frr2(n) + hi23*frr3(n)  ! S-space !
!      fzsn = hi33*frr3(n)

       fxsn = hi11*frr(1,n) + hi12*frr(2,n) + hi13*frr(3,n)
       fysn = hi22*frr(2,n) + hi23*frr(3,n)  ! S-space !
       fzsn = hi33*frr(3,n)

        x1i=x1(n)
        y1i=y1(n)
        z1i=z1(n)

       xmp=x2(n)-dtsqhM*fxsn
       ymp=y2(n)-dtsqhM*fysn
       zmp=z2(n)-dtsqhM*fzsn

c ***  Nose-Hoover Thermostat, 

      sai = sx(n) - xmp*f02_atom
      x1i = x1i - xmp*f12
      sx(n) = sai
      x1(n) = x1i
      x2(n) = x2(n) - xmp
      x3(n) = x3(n) - xmp*f32
      x4(n) = x4(n) - xmp*f42
      x5(n) = x5(n) - xmp*f52

      sai = sy(n) - ymp*f02_atom
      y1i = y1i - ymp*f12
      sy(n) = sai
      y1(n) = y1i
      y2(n) = y2(n) - ymp
      y3(n) = y3(n) - ymp*f32
      y4(n) = y4(n) - ymp*f42
      y5(n) = y5(n) - ymp*f52

      sai = sz(n) - zmp*f02_atom
      z1i = z1i - zmp*f12
      sz(n) = sai
      z1(n) = z1i
      z2(n) = z2(n) - zmp
      z3(n) = z3(n) - zmp*f32
      z4(n) = z4(n) - zmp*f42
      z5(n) = z5(n) - zmp*f52

      Vxt = h11*x1i + h12*y1i + h13*z1i
      Vyt = h22*y1i + h23*z1i
      Vzt = h33*z1i

      sumPx(ntp) = sumPx(ntp) + am*Vxt
      sumPy(ntp) = sumPy(ntp) + am*Vyt
      sumPz(ntp) = sumPz(ntp) + am*Vzt

      Ek_atm = am*(Vxt**2 + Vyt**2 + Vzt**2)
      E_kin(ntp) = E_kin(ntp) + Ek_atm   ! mv^2 !

      enddo ! do n=1,natoms !

!$OMP END PARALLEL DO

      return
      end subroutine  ! correct_atoms !
c
c ------------------------------------------------------------------
c
      END MODULE  ! MD !
