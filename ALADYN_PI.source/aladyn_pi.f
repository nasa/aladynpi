!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! 11-20-2020
!
! aladyn_pi code:
!
! MiniApp for a simple molecular dynamics simulation of an atomistic 
! structure using Physically Informed Neural Network (PINN) potential.
!
! References:
! 1. Vesselin I. Yamakov, Edward H. Glaessgen, NASA/TM-2019-220431
!    (https://ntrs.nasa.gov/citations/20200000069)
! 2. G.P.Purja Pun, V.Yamakov, J.Hickman, E.H.Glaessgen, Y.Mishin, 
!    Physical Review Materials, 4, 113807 (2020).
!
! Source Files: 
!    aladyn_pi_sys.f   - system modul
!    aladyn_pi_sys_OMP.f    - system modul for OpneMP compilation
!    aladyn_pi_sys_NO_OMP.f - system modul without OpneMP compilation
!    aladyn_pi_sys_ACC.f    - system modul for OpneACC compilation
!    aladyn_pi_mods.f  - contains general purpose moduli
!    aladyn_pi_IO.f    - I/O operations
!    aladyn_pi_ANN_OMP.f  - Artificial Neural Network OpneMP code
!    aladyn_pi_ANN_ACC.f  - Artificial Neural Network OpneACC code
!    aladyn_pi_PINN_OMP.f - Physically Informed NN OpneMP code
!    aladyn_pi_PINN_ACC.f - Physically Informed NN OpneACC code
!    aladyn_pi_MD.f    - Molecular Dynamics integrator and subroutines
!    aladyn_pi.f       - Main program
!
! Compilation: 
! make -f makefile.intel OMP=TRUE   ! compile with INTEL with OpenMP
! make -f makefile.pgi ACC=TRUE     ! compile with PGI with OpenACC
!
! Execution:   
! ./aladyn_pi -n 100 -m 10  ! Run 100 MD steps, reporting every 10 steps
!
! For further information contact:
!
! Vesselin Yamakov
! National Institute of Aerospace
! 100 Exploration Way,
! Hampton, VA 23666
! phone: (757)-864-2850
! fax:   (757)-864-8911
! e-mail: yamakov@nianet.org
!
!------------------------------------------------------------------------
! Notices:
! Copyright 2020 United States Government as represented by the 
! Administrator of the National Aeronautics and Space Administration. 
! All Rights Reserved.
! 
! Disclaimers:
! No Warranty: THE SUBJECT SOFTWARE IS PROVIDED "AS IS" WITHOUT ANY 
! WARRANTY OF ANY KIND, EITHER EXPRESSED, IMPLIED, OR STATUTORY, 
! INCLUDING, BUT NOT LIMITED TO, ANY WARRANTY THAT THE SUBJECT SOFTWARE 
! WILL CONFORM TO SPECIFICATIONS, ANY IMPLIED WARRANTIES OF 
! MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, OR FREEDOM FROM 
! INFRINGEMENT, ANY WARRANTY THAT THE SUBJECT SOFTWARE WILL BE ERROR 
! FREE, OR ANY WARRANTY THAT DOCUMENTATION, IF PROVIDED, WILL CONFORM 
! TO THE SUBJECT SOFTWARE. THIS AGREEMENT DOES NOT, IN ANY MANNER, 
! CONSTITUTE AN ENDORSEMENT BY GOVERNMENT AGENCY OR ANY PRIOR RECIPIENT 
! OF ANY RESULTS, RESULTING DESIGNS, HARDWARE, SOFTWARE PRODUCTS OR ANY 
! OTHER APPLICATIONS RESULTING FROM USE OF THE SUBJECT SOFTWARE.  
! FURTHER, GOVERNMENT AGENCY DISCLAIMS ALL WARRANTIES AND LIABILITIES 
! REGARDING THIRD-PARTY SOFTWARE, IF PRESENT IN THE ORIGINAL SOFTWARE, 
! AND DISTRIBUTES IT "AS IS." 
! 
! Waiver and Indemnity:  
! RECIPIENT AGREES TO WAIVE ANY AND ALL CLAIMS AGAINST THE UNITED STATES 
! GOVERNMENT, ITS CONTRACTORS AND SUBCONTRACTORS, AS WELL AS ANY PRIOR 
! RECIPIENT. IF RECIPIENT'S USE OF THE SUBJECT SOFTWARE RESULTS IN ANY 
! LIABILITIES, DEMANDS, DAMAGES, EXPENSES OR LOSSES ARISING FROM SUCH 
! USE, INCLUDING ANY DAMAGES FROM PRODUCTS BASED ON, OR RESULTING FROM, 
! RECIPIENT'S USE OF THE SUBJECT SOFTWARE, RECIPIENT SHALL INDEMNIFY AND 
! HOLD HARMLESS THE UNITED STATES GOVERNMENT, ITS CONTRACTORS AND 
! SUBCONTRACTORS, AS WELL AS ANY PRIOR RECIPIENT, TO THE EXTENT PERMITTED 
! BY LAW. RECIPIENT'S SOLE REMEDY FOR ANY SUCH MATTER SHALL BE THE 
! IMMEDIATE, UNILATERAL TERMINATION OF THIS AGREEMENT.
!------------------------------------------------------------------------
!

      program aladyn_pi

      use sys_OMP
      use sys_ACC
      use constants
      use sim_box
      use IO

      mynod = 0; mpi_nodes = 1

       call read_Args

       call node_management ! Manages node architecture and devices !

       call read_pot       ! POTENTIAL from GMU pot set of files !

       call read_com       ! READ SIMULATION OPTIONS !

       call read_structure ! READ ATOMIC STRUCTURE after read_pot !
                           ! and keeps s-coord. from here on      !
                           ! first call alloc_atoms            !
       call init_param     ! Calls init_vel

       call SIM_run

       write(6,*)'NORMAL termination of the program.'
       call PROGRAM_END(0)

      END     ! MAIN PROGRAM !
!
!---------------------------------------------------------------------
! Check node resources 
!---------------------------------------------------------------------
!
      subroutine check_resources(node_info)

      use sys_OMP
      use sys_ACC
      use sim_box

       type(node_conf), intent(out) :: node_info

       ! Those are replacements of ACC_* equivalents    !
       ! redefined in pgmc_sys_ACC.f and pgmc_sys_OMP.f !

       devicetype = GET_DEVICE_TYPE()
       nACC_devices = GET_NUM_DEVICES(devicetype)
       I_have_GPU = nACC_devices

C      devicetype=1; I_have_GPU=1  ! Test GPU code without GPU VVVV !

       My_GPU_id = -1
       if(nACC_devices.gt.0) My_GPU_id=0

       ! Those are replacements of OMP_* equivalents    !
       ! redefined in pgmc_sys_OMP.f and pgmc_sys_ACC.f !

       MP_procs = GET_NUM_PROCS()
       MP_max_threads = GET_MAX_THREADS()
!      MP_threads = GET_NUM_THREADS()
       MP_threads = MP_max_threads

       node_info%MP_procs = MP_procs
       node_info%MP_threads = MP_threads
       node_info%I_have_GPU = I_have_GPU
       node_info%My_GPU_id = My_GPU_id
       node_info%devicetype = devicetype
       node_info%nACC_devices = nACC_devices

      return
      end ! check_resources !
!
!---------------------------------------------------------------------
! Check for host node resources 
!---------------------------------------------------------------------
!
      subroutine get_node_config(node_info)

      use sys_OMP
      use sys_ACC  ! defines devicetype !
      use sim_box

       type(node_conf), intent(out) :: node_info

       call check_resources(node_info)
       node_info%node_name = 'my_node'

      return
      end ! get_node_config !
!
!---------------------------------------------------------------------
! Report node resources 
!---------------------------------------------------------------------
!
      subroutine report_node_config(node_info)

      use sys_OMP
      use sys_ACC  ! defines devicetype !
      use sim_box

       type(node_conf), intent(in) :: node_info

       MP_procs = node_info%MP_procs
       MP_threads = node_info%MP_threads
       I_have_GPU = node_info%I_have_GPU
       My_GPU_id = node_info%My_GPU_id
       devicetype = node_info%devicetype
       nACC_devices = node_info%nACC_devices

       write(6,*)' '
       write(6,*)'-----------------------------------------------'
       write(6,*)'| '
       write(6,*)'| Node Resources:'

       ! Those are replacements of ACC_* equivalents    !
       ! redefined in pgmc_sys_ACC.f and pgmc_sys_OMP.f !

       if(I_have_GPU.gt.0) then
        call SET_DEVICE_NUM(My_GPU_id, devicetype)
        call GPU_Init(acc_device_current)
        My_GPU_mem = GET_GPU_MEM(My_GPU_id)
        My_GPU_free_mem = GET_GPU_FREE_MEM(My_GPU_id)
        write(6,20)I_have_GPU, My_GPU_id, My_GPU_mem, My_GPU_free_mem
       else 
        write(6,21)
       endif

  20   format(' | GPUs detected',i3,/,' | My_GPU_id=',i3,
     1 ' with memory of:',i12,' bytes, free:',i12)
  21   format(' | No GPU detected.')

        select case(My_GPU_id)
         case(-1)
       write(6,10) MP_procs,MP_threads
         case(0)
       write(6,11) MP_procs,MP_threads,My_GPU_id+1,devicetype
         case(1)
       write(6,12) MP_procs,MP_threads,My_GPU_id+1,devicetype
         case(2)
       write(6,13) MP_procs,MP_threads,My_GPU_id+1,devicetype
         case default
       write(6,14) MP_procs,MP_threads,My_GPU_id+1,devicetype
        end select

  10   format(' | CPUs:',i4,' using ',i4,' threads',
     1 ' and no GPU devices')
  11   format(' | CPUs:',i4,' using ',i4,' threads',
     1 /,' | and the ',i1,'-st GPU of devicetype=',i3)
  12   format(' | CPUs:',i4,' using ',i4,' threads',
     1 /,' | and the ',i1,'-nd GPU of devicetype=',i3)
  13   format(' | CPUs:',i4,' using ',i4,' threads',
     1 /,' | and the ',i1,'-rd GPU of devicetype=',i3)
  14   format(' | CPUs:',i4,' using ',i4,' threads',
     1 /,' | and the ',i2,'-th GPU of devicetype=',i3)

       write(6,*)'|' 
       write(6,*)'-----------------------------------------------'
       write(6,*)' '

      return
      end ! report_node_config !
!
!---------------------------------------------------------------------
! Pre-calculates the many-body part of the interatomic potential
!---------------------------------------------------------------------
!
      subroutine node_management

      use sys_OMP
      use sys_ACC
      use sim_box
      use IO

      type(node_conf) :: node_info

      character(40) my_node_name, node_name, name_of_node,
     1 name_in_group
      integer :: ngroup_of_node(1)

      call get_node_config(node_info)
      call report_node_config(node_info)

      call alloc_nodes(1,ierror) 
      call error_check(ierror,'ERROR in alloc_nodes...')

      nxold=0; nyold=0; nzold=0

      nodes_on_Y = 0
      nodes_on_Z = 0

      return    ! node_management !
      end
!
!---------------------------------------------------------------------
!  ! DO NOT CALL ALONE, only through force_global !
!---------------------------------------------------------------------
!
      subroutine force(ienergy)

      use sys_OMP
      use sys_ACC
      use sim_box
      use pot_module
      use IO
      use ANN
      use ANN_ACC
      use PINN
      use PINN_ACC

      implicit double precision (a-h,o-z)

      integer, intent(in) :: ienergy

      select case(iPOT_file_ver)
       case(5) ! ANN potential !
        if(I_have_GPU.gt.0) then
         call Frc_ANN_ACC(ecoh)  ! Analytical derivatives !
        else
         call Frc_ANN_OMP(ecoh)
        endif

       case(6) ! PINN potential !
        if(I_have_GPU.gt.0) then
         call nbrs_BOP_ACC_inv
         call Frc_PINN_ACC(ecoh)  ! Analytical derivatives !
        else
         call get_nbrs_inv
         call Frc_PINN_OMP(ecoh)
        endif
       end select

      call error_check(ihalt,'ERROR in force().')

      return    ! force !
      end
!
!---------------------------------------------------------------------
! At start up:
!
! call read_pot 
!      call read_pot_dat ! sets iatom_types, ipair_types,
!         	       ! iPOT_file_ver
!
!      ! for ANN potential types !	
!           call input_pot_ANN ! read ANN pot. file !
!           call alloc_types_ANN
!
!      call calc_pot_param
!           call init_param_ANN
!
!---------------------------------------------------------------------
!
      subroutine read_pot

      use pot_module
      use IO
      use ANN

      implicit double precision (a-h,o-z)

      call init_elements

      ierror = 1
      call read_pot_dat ! get iatom_types,ifile_numbs,call alloc_types !

      call input_pot_ANN(ierror)
         ! alloc_ates ANN types, alloc_ates ANN_BOP types !

      call error_check(ierror,'ERROR in read_pot...')

      call init_param_ANN ! rescales or initializes pot. param. !
                          ! for ANN_BOP, alloc_ates BOP types too!
      call error_check(ihalt,'ERROR in read_pot.')

      return     
      end        ! read_pot !
!
!---------------------------------------------------------------------
!
      subroutine init_param

      use constants
      use sim_box
      use atoms
      use pot_module
      use MD
      use IO

      implicit double precision (a-h,o-z)

      character seed_string*15

      istep = 0
      acc_V_rate = 0.D0
      PotEnrg_glb = 0.d0
      ecoh = 0.d0
!
!   ***  define the inverse, hi(1..3,1..3), of h(1..3,1..3)
!
! now add the periodic images of the input geometry
!
      ibox_error = 0
      if(abs(2.d0*h(1,2)).gt.h(1,1)) ibox_error = 1
      if(abs(2.d0*h(1,3)).gt.h(1,1)) ibox_error = 2
      if(abs(2.d0*h(1,3)).gt.h(2,2)) ibox_error = 3
      if(abs(2.d0*h(2,3)).gt.h(1,1)) ibox_error = 4
      if(abs(2.d0*h(2,3)).gt.h(2,2)) ibox_error = 5
      if(ibox_error.ne.0) then 
       write(6,*)'ERROR! Unacceptable h-matrix!'
       select case(ibox_error)
        case (1) 
         write(6,*)'h(1,2) or xy must be less than (xhi-xlo)/2'
        case (2) 
         write(6,*)'h(1,3) or xz must be less than (xhi-xlo)/2'
        case (3) 
         write(6,*)'h(1,3) or xz must be less than (yhi-ylo)/2'
        case (4) 
         write(6,*)'h(2,3) or xz must be less than (xhi-xlo)/2'
        case (5) 
         write(6,*)'h(2,3) or xz must be less than (yhi-ylo)/2'
       end select
       call PROGRAM_END(1)
      endif

      call matinv(h,hi,dh) ! h*hi=I

       do ntp = 1,iatom_types  
        Am_of_type(ntp) = 0.0d0
        pTemp(ntp) = T_sys
        E_kin(ntp) = 1.5d0*Boltz_Kb*T_sys
       enddo

       do n=1,natoms
        ntp = ntype(n)
        Am_of_type(ntp) = Am_of_type(ntp) + amass(ntp)
       enddo   

       sum_mass(1:iatom_types) = Am_of_type(1:iatom_types)

       avr_mass = 0.0d0
       do iatom = 1,iatom_types   ! possible but inefficient vect. !
        avr_mass = avr_mass + sum_mass(iatom)
       enddo
       avr_mass = avr_mass/natoms

! *** set up the random number generator

      iseed=6751
      ! Convert the string 'str_filename' to a numeric value !
       n = len(str_filename)
       do i=1,n 
        iseed = iseed + ICHAR(str_filename(i:i))  
       enddo          ! randomize with str_filename !

       iseed0 = iseed+17
       call rmarin(iseed0) ! node depndnt init. random generator

      call init_MD

      return
      end                 ! init_param !
!
! -------------------------------------------------------------------
!
      subroutine SIM_run

      use constants
      use sim_box
      use pot_module
      use atoms
      use MD
      use IO
      use ANN

      implicit double precision (a-h,o-z)
      
      call init_vel(T_set)  ! VVVV TEST !
      call force_global(1) ! err. check node_config finder   !
      call initaccel               ! sets accelerations using forces !

      write(6,*)' '
      write(6,*)'PotEnrg_atm=',PotEnrg_atm
      write(6,*)'Sys. Pot.En=',PotEnrg_atm*natoms

      call report(0) ! Initial structure measurement !

      BkT = 1.0d0/(Boltz_Kb*T_sys)
!     
!  ******************************************************************
!  ***  begin do loop for MD steps                                  *
!  ******************************************************************
!   
      istep = 0; 

      do kstep=1,nstep    ! MD loop !

       istep = istep + 1

       E1 = PotEnrg_glb/natoms

! --- MD step start ---
        real_time = real_time + 1000.d0*dt_step  ! [fs] MD run !

        call predict_atoms(ndof_flag)

        call force_global(0) ! no node_config !

        call correct_atoms(ndof_flag)   ! calc. sumPxyz() !
        call T_broadcast ! Send A_fr, sumPxyz() and calc. pTemp(ntp) !

! --- MD step end ---

       if ((mod(kstep,measure_step).eq.0).AND.(kstep.lt.nstep)) then 
        call force_global(0)
        call report(kstep) 
       endif

      enddo ! do kstep=1,nstep ! end of MD loop !

      call get_chem

      ! Calc. Final Energy w stress !
      call force_global(0)  ! calc atm stress !
      call report(kstep-1)
      call write_structure_plt
!
!   ******************************************************************
!   ***  end of do loop for time steps  ******************************
!   ******************************************************************
!
      return  
      end        ! SIM_run !
!
! -------------------------------------------------------------------
! Finds the best number of cells in a given direction (nnx,nny,nnz),
! which commensurates with MC_rank and nodes_on_X,Y,Z
! -------------------------------------------------------------------
!
      integer function nnd_fit(nodes_on_D, iD, nnd_min, MC_rank_D)

      use sim_box

      implicit double precision (a-h,o-z)

      integer, intent(in) :: nodes_on_D, iD
      integer, intent(out) :: nnd_min, MC_rank_D

      nnd_min = 3
      nnd_max = int(h(iD,iD)/size)

!     write(50,*)h
!     write(50,*)'size=',size

      nnd_tmp = 0; MC_rank_D = MC_rank

      do nnd = nnd_max, nnd_min, -1

       nnd_nodes = mod(nnd,nodes_on_D) ! check nnd vs nodes_on_D !

       if(nnd_nodes.eq.0) then         ! nnd fits on nodes_on_D !
        nnd_rank = mod(nnd,MC_rank_D)  ! check nnd with MC_rank_D !
        if(nnd_rank.eq.0) nnd_tmp=nnd  ! good nnd found !
        do while(nnd_rank.gt.0)        ! if doesn't do...         !
         if((MC_rank_D.lt.nnd_min).AND.(MC_rank_D.lt.MC_rank_max)) then
          MC_rank_D = MC_rank_D + 1
          nnd_rank = mod(nnd,MC_rank_D) ! check with new MC_rank_D !
          if(nnd_rank.eq.0) nnd_tmp=nnd ! good nnd found !
         else
          nnd_rank = 0  ! stop do while() !
         endif
        enddo ! do while(irepeat.gt.0) !

        if(nnd_tmp.gt.0) exit          ! exit do nnd loop !
       endif  ! if(nnd_nodes.eq.0)... !

      enddo  ! do nnd = nnd_max, nnd_min, -1 !

      nnd_fit = nnd_tmp

!     write(50,*)'d:',iD,' nnd_fit=',nnd_fit,' nnd_min=',nnd_min
!     write(50,*)'      MC_rank_D=',MC_rank_D

      return
      end function   ! nnd_fit(nnn) !
!
! -------------------------------------------------------------------
!  Looks for optimal node architecture configuration at a given
!  number of nodes on X:nodes_X, on Y:nodes_Y, and on Z:nodes_Z
! -------------------------------------------------------------------
!
      subroutine get_config(i1,i2,i3, nodes_X, nodes_Y, nodes_Z,
     1           nnx_min,nny_min,nnz_min, nnx_cell,nny_cell,nnz_cell, 
     1           ierror)

      use sim_box

      implicit double precision (a-h,o-z)

      integer, intent(in) :: i1,i2,i3, nodes_X, nodes_Y, nodes_Z
      integer, intent(out) :: nnx_min,nny_min,nnz_min,
     1 nnx_cell,nny_cell,nnz_cell,ierror

!      write(50,10) i1,i2,i3, nodes_X, nodes_Y, nodes_Z
!  10  format('get_config(',6i3,')')

       ierror = 0
       nnx = nnd_fit(nodes_X, i1, nnx_min, MC_rank_X)
       if(nnx.eq.0) ierror=IOR(ierror,1)
       nny = nnd_fit(nodes_Y, i2, nny_min, MC_rank_Y)
       if(nny.eq.0) ierror=IOR(ierror,2)
       nnz = nnd_fit(nodes_Z, i3, nnz_min, MC_rank_Z)
       if(nnz.eq.0) ierror=IOR(ierror,4)

!      write(50,*)'get_conf: nnx,y,z=',nnx,nny,nnz,ierror
!      write(50,*)'get_conf: MC_rank_X,Y,Z=',
!    1 MC_rank_X,MC_rank_Y,MC_rank_Z

      nnx_cell = nnx / nodes_X
      nny_cell = nny / nodes_Y
      nnz_cell = nnz / nodes_Z

      ! Check if cells per node commensurate with MC_ranks !

       nn_mod = mod(nnx_cell,MC_rank_X)
       if(nn_mod.ne.0) ierror=IOR(ierror,1)
       nn_mod = mod(nny_cell,MC_rank_Y)
       if(nn_mod.ne.0) ierror=IOR(ierror,2)
       nn_mod = mod(nnz_cell,MC_rank_Z)
       if(nn_mod.ne.0) ierror=IOR(ierror,4)

!      write(50,*)'get_conf: nnx,y,z_cell=',nnx_cell,nny_cell,
!    1 nnz_cell,ierror
!      write(50,*)' '

      return
      end       ! get_config !
!
! -------------------------------------------------------------------
!   Looks for optimal node architecture configuration
! -------------------------------------------------------------------
!
      subroutine node_config(nflag)
! 
!   ***  updates hij matrix according to the fardest atoms in the
!        system in Y and Z directions, and
!   ***  looks for optimal node architecture configuration
!
      use sim_box
      use pot_module
      use IO
      use atoms

      implicit double precision (a-h,o-z)

      integer, intent(out) :: nflag
      integer natoms_alloc_new ! local !
      integer nny_min,nnz_min

      nflag = 0
      ierror = 0
      nodes_on_Yo=nodes_on_Y; nodes_on_Zo=nodes_on_Z
      cell_size_Y = 0.d0; cell_size_Z = 0.d0

       i1=1; i2=2; i3=3

       nodes_on_Y=1; nodes_on_Z=1
       call get_config(i1,i2,i3, 1,1,1, nnx_min,nny_min,nnz_min,
     1      nnx_cell,nny_cell,nnz_cell,ierror)
            nny_try = 1; nnz_try = 1

      if(ierror.gt.0) then

        write(6,*)' '
      write(6,*)'ERROR: Unable to construct a suitable link-cell grid!'
        write(6,*)' '
        write(6,99) h(i1,i1),h(i2,i2),h(i3,i3),size
        if(IAND(ierror,1).gt.0) then
         write(6,47)nnx
        else
         write(6,37)nnx
        endif
        if(IAND(ierror,2).gt.0) then
         write(6,48)nny,nny_min-1
        else
         write(6,38)nny,nny_min-1
        endif
        if(IAND(ierror,4).gt.0) then
         write(6,49)nnz,nnz_min-1
        else
         write(6,39)nnz,nnz_min-1
        endif
        write(6,*)'Decrease number of nodes or increase system size...'

       call PROGRAM_END(1)

      else  ! if(ierror.gt.0)... !

       cell_size_X=h(i1,i1)/nnx
       cell_size_Y=h(i2,i2)/nny
       cell_size_Z=h(i3,i3)/nnz

      endif ! if(ierror.gt.0)... !

 37   format(' Cells per node on X =',i4,'  must be  >   2: YES')
 47   format(' Cells per node on X =',i4,'  must be  >   2: NO')
 38   format(' Cells per node on Y =',i4,'  must be  > ',i3,': YES')
 48   format(' Cells per node on Y =',i4,'  must be  > ',i3,': NO')
 39   format(' Cells per node on Z =',i4,'  must be  > ',i3,': YES')
 49   format(' Cells per node on Z =',i4,'  must be  > ',i3,': NO')
 99   format('System Box size:',2(f12.6,' x '),f12.6,
     1 ';   min. cell size=',f12.8)

      ncell=nnx*nny*nnz

      call matinv(h,hi,sys_vol)
      atom_vol1 = sys_vol/natoms

      ! Get the maximum possible number of atoms per link cell ! 

      rZ_min = 10.0d0 ! Start max value in Ang. !

C     write(6,*)'nelem_in_com=',nelem_in_com

      do i=1,nelem_in_com  ! read el. types as listed in aladyn_pi.com !
       iZ = iZ_elem_in_com(i)
       rad_of_Z = elem_radius(iZ) ! [Ang] !

C      write(6,*)'rad_of_Z(iZ=',iZ,')=',rad_of_Z

       if(rad_of_Z.lt.rZ_min) rZ_min=rad_of_Z
      enddo

      atom_vol2 = 4.0/3.0*3.141592*rZ_min**3

C     write(6,*)'sys_vol=',sys_vol,' r_max=',r_max,
C    1 ' rZ_min=',rZ_min
C     write(6,*)'atom_vol1=',atom_vol1,' atom_vol2',atom_vol2

      if(atom_vol1.lt.atom_vol2) then
        atom_vol = atom_vol1
      else
        atom_vol = atom_vol2
      endif

C     write(6,*)'atom_vol=',atom_vol,' cell_volume=',
C    1 cell_volume,' rZ_min=',rZ_min

      cell_volume = (cell_size_X+2.d0*rZ_min)*
     1 (cell_size_Y+2.d0*rZ_min)*(cell_size_Z+2.d0*rZ_min)
      cell3_volume = (3.d0*cell_size_X+2.d0*rZ_min)*
     1 (3.d0*cell_size_Y+2.d0*rZ_min)*(3.d0*cell_size_Z+2.d0*rZ_min)
      natoms_per_cell = int(cell_volume/atom_vol) + 1
      natoms_per_cell = (natoms_per_cell/8 + 1)*8 
      natoms_per_cell3 = int(cell3_volume/atom_vol) + 1
      natoms_per_cell3 = (natoms_per_cell3/8+1)*8 

      nflag=abs(nnx-nxold)+abs(nny-nyold)+abs(nnz-nzold)+
     +      abs(nodes_on_Y-nodes_on_Yo)+abs(nodes_on_Z-nodes_on_Zo)

!  reset cell grid if necessary
      if (nflag.gt.0) then
       call link_cell_setup

        write(6,101)ncell

        write(6,102) 1,nnx,cell_size_X
        write(6,103) nodes_on_Y,nny_cell,cell_size_Y
        write(6,104) nodes_on_Z,nnz_cell,cell_size_Z

        write(6,*)' '
      endif    ! if (nflag.gt.0)... !

      natoms_alloc_new = natoms + 100

      if(natoms_alloc_new.gt.natoms_alloc) then ! update natoms_alloc !
       natoms_alloc = (natoms_alloc_new/64 + 1)*64 
      endif

      cut_off_vol = 4.0/3.0*3.141592*(r_cut_off+rZ_min/2.0)**3
      nbrs_per_atom = nint(cut_off_vol/atom_vol) ! Correct one !
      nbrs1_per_atom = (int(nbrs_per_atom*0.3)/8 + 1)*8
      nbrs_per_atom = (nbrs_per_atom/8 + 1)*8
      nbrs_alloc = nbrs_per_atom*natoms_alloc

      if(iPOT_file_ver.eq.6) then ! PINN !
        cut_off_vol1 = 4.0/3.0*3.141592*(r_cut_short+rZ_min)**3
        nbrs1_per_atom = (int(nbrs_per_atom*0.3)/8 + 1)*8
      endif

C     write(6,110)nbrs1_per_atom,nbrs_per_atom,cut_off_vol

      nxold=nnx
      nyold=nny
      nzold=nnz

101   format(/,'Link cell configuration:',/,
     1 ' axis nodes cells/n thickness; total cell:',i6)
102   format('On X: ',i3,' x ',i3,' x ',f6.3)
103   format('On Y: ',i3,' x ',i3,' x ',f6.3)
104   format('On Z: ',i3,' x ',i3,' x ',f6.3)
105   format(/,'Allocated atoms:',i8)
110   format(/,'nbrs1_per_atom=',i4,' nbrs_per_atom=',i4,
     1 ' cut_off_vol=',f12.4,/)

      if (NOT(allocated(ident))) then
       call alloc_atoms ! alloc_ates natoms_alloc atoms in aladyn_pi_IO!
      endif

      return
      end       ! node_config !
!
! -------------------------------------------------------------------
!
      subroutine force_global(ilong)
!
!   ***  subroutine for doing force calculation with linked cells
!        ilong = (1) do full force calculation
!                (0) skip node_config and finder when
!                    sx() haven't changed, like in box rescaling
      use sim_box
      use pot_module
      use IO

      implicit double precision (a-h,o-z)

      if(ilong.ne.0) then
       call node_config(nflag) ! Update ncell,ncell_per_node,
                                        ! natoms_alloc
       call finder                      ! Redistribute atoms to nodes !
      endif ! if(ilong.ne.0)... !

      select case(iPOT_file_ver)
       case(5) ! ANN potential !
        call get_ANN_nbrs
       case(6) ! PINN potential !
        call get_PINN_nbrs
       end select

      call force(ienergy)

! --- Sum Pot. Energy from all nodes ---

      PotEnrg_glb = ecoh
      PotEnrg_atm = PotEnrg_glb/natoms

!     call PROGRAM_END(1)  ! VVVV !

      return
      end   ! force_global !
!
!---------------------------------------------------------------------
! Performs loop over atoms in ALL of the neighboring cells
! Uses OMP, but prepares the system for ACC (GPU) implementation.
!---------------------------------------------------------------------
!
      subroutine get_ANN_nbrs

      use sim_box
      use atoms
      use pot_module

      implicit double precision (a-h,o-z)

      integer, dimension(natoms_per_cell3) :: ll_nbr
!
! do loop over all cells
!     
      h11 = h(1,1); h12 = h(1,2); h13 = h(1,3)
      h22 = h(2,2); h23 = h(2,3); h33 = h(3,3)

!$OMP  PARALLEL DO DEFAULT(NONE) PRIVATE (j)
!$OMP& SHARED(natoms, nbrs_per_atom, nbr_list)
!$OMP& SCHEDULE(STATIC)
      do i=1,natoms
        !DIR$ VECTOR UNALIGNED
        do j=0,nbrs_per_atom
          nbr_list(j,i) = i   ! Initial state: all nbrs are self-nbrs !
        enddo
      enddo

      max_nbrs = 0
      sz0_cut = r_cut_off/h33

!$OMP PARALLEL DO DEFAULT(NONE) REDUCTION(max:max_nbrs)
!$OMP& PRIVATE(ic,icell,jcell,iz_nr,iyx,iy_nr,ixx,nr_in_cell,n,ll_nbr)
!$OMP& PRIVATE(l_in_cell,izl,iyl,kzn,jyn,jyl,ns,nr,sxn,syn,szn,l)
!$OMP& PRIVATE(sx0,sy0,sz0,rx0,ry0,rz0,r2,k_all,k1_all,k2_all)
!$OMP& SHARED(natoms, nnx,nny,nnz, nXYlayer, ncells_all, nbr_list)
!$OMP& SHARED(h11,h22,h33,h12,h13,h23,r2_cut_off,sz0_cut)
!$OMP& SHARED(mpi_nodes, ncY_per_node, ncZ_per_node, nodes_on_Y)
!$OMP& SHARED(id_of_cell, natoms_in_cell, n_in_cell, sx,sy,sz)
!$OMP& SCHEDULE(STATIC)

      do_cells: do ic = 1, ncells_all  ! Each ic is independent !
       icell = id_of_cell(ic)
       iz_nr = icell/nXYlayer
       iyx = mod(icell,nXYlayer)
       iy_nr = iyx/nnx
       ixx = mod(iyx,nnx) + nnx

       nr_in_cell = natoms_in_cell(icell)    ! number atoms in icell !
       do n = 1,nr_in_cell              ! VECTORIZED: speedup: 4.760 !
        ll_nbr(n) = n_in_cell(n,icell)
       enddo

       l_in_cell = nr_in_cell
!
! Loop over atoms in all 27 neighboring cells
!
! 0 <= izl <= ncZ_per_node - 1 : range=ncZ_per_node
! 0 <= iyl <= ncY_per_node - 1 : range=ncY_per_node
! 
       do izl=-1,1
        kzn = mod(iz_nr+nnz+izl,nnz)*nXYlayer
       do iyl=-1,1
        jyn= kzn + mod(iy_nr+nny+iyl,nny)*nnx
       do i=-1,1
        jcell= jyn + mod(i+ixx,nnx)
        if(jcell.ne.icell) then
         ns = natoms_in_cell(jcell)     ! number atoms in jcell !
         do n = 1,ns                   ! VECTORIZED: speedup: 4.760 !
          ll_nbr(l_in_cell+n) = n_in_cell(n,jcell) ! atom n in jcell !
         enddo
         l_in_cell = l_in_cell+ns
        endif ! if(icell.ne.jcell)... !
       enddo ! do i=-1,1 !
       enddo ! do iyl = !
       enddo ! do izl = !
!
! Start: Find neighbors of atoms in icell.
!
       do n=1,nr_in_cell

        nr = ll_nbr(n)

        sxn=sx(nr)
        syn=sy(nr)
        szn=sz(nr)

        k_all = 0  ! no neighbors !

        do k = 1, l_in_cell
         l = ll_nbr(k)
         sz0 = sz(l) - szn
         if (sz0.ge. 0.5d0) then   ! make periodic along Z !
             sz0=sz0-1.0d0         ! 2x faster than sz0=sz0-dnint(sz0) !
         elseif (sz0.lt.-0.5d0) then
             sz0=sz0+1.0d0
         endif
         if (abs(sz0).lt.sz0_cut) then
          rz0 = h33*sz0
          sy0 = sy(l) - syn
          if (sy0.ge. 0.5d0) then   ! make periodic along Y !
              sy0=sy0-1.0d0
          elseif (sy0.lt.-0.5d0) then
              sy0=sy0+1.0d0
          endif
          ry0 = h22*sy0 + h23*sz0
          sx0 = sx(l) - sxn
          if (sx0.ge. 0.5d0) then   ! make periodic along X !
              sx0=sx0-1.0d0
          elseif (sx0.lt.-0.5d0) then
              sx0=sx0+1.0d0
          endif
          rx0 = h11*sx0 + h12*sy0 + h13*sz0
          r2 = rx0**2 + ry0**2 + rz0**2

          if ((r2.lt.r2_cut_off).AND.(l.ne.nr)) then
           k_all = k_all + 1
           nbr_list(k_all,nr) = l
          endif ! if (r2.lt.r2_cut_off)... 
         endif ! if (abs(sz0).lt.r_cut_off)... !
        enddo ! do do k = 1, l_in_cell !

        max_nbrs = max(k_all,max_nbrs)

       enddo ! do n=1,nr_in_cell

      enddo do_cells ! do ic = 1, ncells_all !

!$OMP END PARALLEL DO

      ! Ensure max_nbrs is a multiple of 8 to avoid remainder loops 
      ! after vectorization

      if (mod(max_nbrs,8) .ne. 0) max_nbrs = ((max_nbrs/8)+1)*8

      return
      end          ! get_ANN_nbrs !
!
!---------------------------------------------------------------------
! Performs loop over atoms in ALL of the neighboring cells
! dividing them into two gropus of:
! 1. neighbors in r < Rc range, in (1.. k1_all)
! 2. screening neighbors in (Rc < r < 1.5Rc) in nbrs_2nd_of(k1_all+1.. )
! Uses OMP, but prepares the system for ACC (GPU) implementation.
!---------------------------------------------------------------------
!
      subroutine get_PINN_nbrs

      use sim_box
      use atoms
      use pot_module
      use PINN

      implicit double precision (a-h,o-z)

      integer, dimension(natoms_per_cell3) :: ll_nbr
!
! do loop over all cells
!     
      h11 = h(1,1); h12 = h(1,2); h13 = h(1,3)
      h22 = h(2,2); h23 = h(2,3); h33 = h(3,3)

!$OMP  PARALLEL DO DEFAULT(NONE) PRIVATE (j)
!$OMP& SHARED(natoms, nbrs_per_atom, nbr_list)
!$OMP& SCHEDULE(STATIC)
      do i=1,natoms
        !DIR$ VECTOR UNALIGNED
        do j=0,nbrs_per_atom
          nbr_list(j,i) = i   ! Initial state: all nbrs are self-nbrs !
        enddo
      enddo

      max1_nbrs=0; max2_nbrs=0
      sz0_cut = r_cut_off/h33

!$OMP PARALLEL DO DEFAULT(NONE) REDUCTION(max:max1_nbrs)
!$OMP& PRIVATE(ic,icell,jcell,iz_nr,iyx,iy_nr,ixx,nr_in_cell,n,ll_nbr)
!$OMP& PRIVATE(l_in_cell,izl,iyl,kzn,jyn,jyl,ns,nr,sxn,syn,szn,l)
!$OMP& PRIVATE(sx0,sy0,sz0,rx0,ry0,rz0,r2,k_all,k1_all,k2_all)
!$OMP& SHARED(natoms, nnx,nny,nnz, nXYlayer, ncells_all, nbr_list)
!$OMP& SHARED(h11,h22,h33,h12,h13,h23,r2_cut_off,r2_cut_short,sz0_cut)
!$OMP& SHARED(nbrs_all_of, nbrs_2nd_of, nbrs_per_atom)
!$OMP& SHARED(mpi_nodes, ncY_per_node,ncZ_per_node,nodes_on_Y)
!$OMP& SHARED(id_of_cell, natoms_in_cell, n_in_cell, sx,sy,sz)
!$OMP& SCHEDULE(STATIC)

      do_cells: do ic = 1, ncells_all  ! Each ic is independent !
       icell = id_of_cell(ic)
       iz_nr = icell/nXYlayer
       iyx = mod(icell,nXYlayer)
       iy_nr = iyx/nnx
       ixx = mod(iyx,nnx) + nnx

       nr_in_cell = natoms_in_cell(icell)    ! number atoms in icell !
       do n = 1,nr_in_cell              ! VECTORIZED: speedup: 4.760 !
        ll_nbr(n) = n_in_cell(n,icell)
       enddo

       l_in_cell = nr_in_cell
!
! Loop over atoms in all 27 neighboring cells
!
! 0 <= izl <= ncZ_per_node - 1 : range=ncZ_per_node
! 0 <= iyl <= ncY_per_node - 1 : range=ncY_per_node
! 
       do izl=-1,1
        kzn = mod(iz_nr+nnz+izl,nnz)*nXYlayer
       do iyl=-1,1
        jyn= kzn + mod(iy_nr+nny+iyl,nny)*nnx
       do i=-1,1
        jcell= jyn + mod(i+ixx,nnx)
        if(jcell.ne.icell) then
         ns = natoms_in_cell(jcell)     ! number atoms in jcell !
         do n = 1,ns                   ! VECTORIZED: speedup: 4.760 !
          ll_nbr(l_in_cell+n) = n_in_cell(n,jcell) ! atom n in jcell !
         enddo
         l_in_cell = l_in_cell+ns
        endif ! if(icell.ne.jcell)... !
       enddo ! do i=-1,1 !
       enddo ! do iyl = !
       enddo ! do izl = !
!
! Start: Find neighbors of atoms in icell.
!
       do n=1,nr_in_cell

        nr = ll_nbr(n)

        sxn=sx(nr); syn=sy(nr); szn=sz(nr)

        k_all = 0; k1_all = 0; k2_all = nbrs_per_atom+1 ! no neighbors !

        do k = 1, l_in_cell
         l = ll_nbr(k)
         if (l.ne.nr) then
          sz0 = sz(l) - szn
          if (sz0.ge. 0.5d0) then  ! make periodic along Z !
              sz0=sz0-1.0d0        ! 2x faster than sz0=sz0-dnint(sz0) !
          elseif (sz0.lt.-0.5d0) then
              sz0=sz0+1.0d0
          endif
          if (abs(sz0).lt.sz0_cut) then
           rz0 = h33*sz0
           sy0 = sy(l) - syn
           if (sy0.ge. 0.5d0) then   ! make periodic along Y !
               sy0=sy0-1.0d0
           elseif (sy0.lt.-0.5d0) then
               sy0=sy0+1.0d0
           endif
           ry0 = h22*sy0 + h23*sz0
           sx0 = sx(l) - sxn
           if (sx0.ge. 0.5d0) then   ! make periodic along X !
               sx0=sx0-1.0d0
           elseif (sx0.lt.-0.5d0) then
               sx0=sx0+1.0d0
           endif
           rx0 = h11*sx0 + h12*sy0 + h13*sz0
           r2 = rx0**2 + ry0**2 + rz0**2
           if (r2.lt.r2_cut_off) then
            if (r2.lt.r2_cut_short) then
             k1_all = k1_all + 1
             k_all = k1_all
            else
             k2_all = k2_all - 1
             k_all = k2_all
            endif
            nbr_list(k_all,nr) = l       ! 0 < rij < 1.5Rc !
           endif ! if (r2.lt.r2_cut_off)... 
          endif ! if (abs(sz0).lt.r_cut_off)... !
         endif ! if (l.ne.nr)... !
        enddo ! do do k = 1, l_in_cell !

        nbrs_2nd_of(nr) = k2_all  ! k2_all,.. nbrs_per_atom : Rc <1.5Rc
        nbrs_all_of(nr) = k1_all  ! may be used in other subroutines !

        max1_nbrs = max(k1_all,max1_nbrs)

       enddo ! do n=1,nr_in_cell
      enddo do_cells ! do ic = 1, ncells_all !

!$OMP END PARALLEL DO

      if (mod(max1_nbrs,8).ne.0) max1_nbrs = ((max1_nbrs/8)+1)*8

c Rearrange screening nbrs to follow next after the interaction nbrs !

       if(I_have_GPU.gt.0) then
      
!$OMP PARALLEL DO DEFAULT(NONE) REDUCTION(max:max2_nbrs)
!$OMP& PRIVATE(new_place)
!$OMP& SHARED(nbrs_per_atom)
!$OMP& SHARED(natoms, max1_nbrs, nbr_list, nbrs_all_of, nbrs_2nd_of)
!$OMP& SCHEDULE(DYNAMIC,CHUNK)

      do nr = 1, natoms
       new_place = max1_nbrs  ! put screening nbrs at max1_nbrs+1,... !
       do k = nbrs_2nd_of(nr), nbrs_per_atom
        new_place = new_place + 1
        nbr_list(new_place,nr) = nbr_list(k,nr)
        nbr_list(k,nr) = nr
       enddo
       max2_nbrs = max(new_place,max2_nbrs) ! Rc < rij < 1.5Rc !
       nbrs_2nd_of(nr)=new_place  
      enddo ! do nr = 1, natoms !

!$OMP END PARALLEL DO

       else ! if(I_have_GPU.gt.0)... !
      
!$OMP PARALLEL DO DEFAULT(NONE) REDUCTION(max:max2_nbrs)
!$OMP& PRIVATE(new_place)
!$OMP& SHARED(nbrs_per_atom)
!$OMP& SHARED(natoms, nbr_list, nbrs_all_of, nbrs_2nd_of)
!$OMP& SCHEDULE(DYNAMIC,CHUNK)

      do nr = 1, natoms
       new_place = nbrs_all_of(nr) ! put screen nbrs after close nbrs !
       do k = nbrs_2nd_of(nr), nbrs_per_atom
        new_place = new_place + 1
        nbr_list(new_place,nr) = nbr_list(k,nr)
        nbr_list(k,nr) = nr
       enddo
       max2_nbrs = max(new_place,max2_nbrs) ! Rc < rij < 1.5Rc !
       nbrs_2nd_of(nr)=new_place  
      enddo ! do nr = 1, natoms !

!$OMP END PARALLEL DO

       endif ! if(I_have_GPU.gt.0)... !
      
      ! Ensure max_nbrs is a multiple of 8 to avoid remainder loops 
      ! after vectorization

      if (mod(max2_nbrs,8).ne.0) max2_nbrs = ((max2_nbrs/8)+1)*8

C     write(6,20) max1_nbrs, max2_nbrs
C 20  format('get_PINN_nbrs: max1_nbrs=',i4,' max2_nbrs=',i4)

      return
      end          ! get_PINN_nbrs !
c
c-----------------------------------------------------------------------
c Makes an invert list of neighbors nabors_num for each atom
c from nbr_list() to create nbr_list_inv()
c Called for PINN case only!
c-----------------------------------------------------------------------
c
      subroutine get_nbrs_inv

      use atoms
      use PINN

      implicit double precision (a-h,o-z)

!$OMP PARALLEL DO DEFAULT(NONE)
!$OMP& PRIVATE(l,ij_delta,nbrs_of_l,nbr_nr,n_nr,l2,nbrs2_of_l)
!$OMP& PRIVATE(nbrs1, nbrs2, nbrs1_l, nbrs2_l)
!$OMP& SHARED(natoms, nbrs_all_of,nbrs_2nd_of, nbr_list,nbr_list_inv)
!$OMP& SCHEDULE(STATIC)

       do nr=1,natoms

        nbrs1 = nbrs_all_of(nr)
        nbrs2 = nbrs_2nd_of(nr)

        do nbl = 1,nbrs1        ! loop over 1st nbrs of nr !
         l = nbr_list(nbl,nr)
         nbrs1_l = nbrs_all_of(l)

c loop over 1st nbrs of l to find the nr position as a neighbor of l

         do n_nr = 1,nbrs1_l       ! VECTORIZED: speedup: 2.390 !
          l2 = nbr_list(n_nr,l)     ! l2 is a 1st neighbor of l !
          if(l2.eq.nr) nbr_nr=n_nr
         enddo  ! do nbr_nr = 1,nbrs1_l !
         nbr_list_inv(nbl,nr) = nbr_nr
        enddo ! do nbl = 1,nbrs1 !

        do nbl = nbrs1+1,nbrs2  ! loop over 2nd nbrs of nr !
         l = nbr_list(nbl,nr)
         nbrs1_l = nbrs_all_of(l)
         nbrs2_l = nbrs_2nd_of(l)
         do n_nr = nbrs1_l+1,nbrs2_l  ! VECTORIZED: speedup: 2.390 !
          l2 = nbr_list(n_nr,l)        ! l2 is a 2nd neighbor of l !
          if(l2.eq.nr) nbr_nr=n_nr
         enddo  ! do nbr_nr = 1,nbrs2_l !
         nbr_list_inv(nbl,nr) = nbr_nr
        enddo ! do nbl = nbrs1_l+1,nbrs2 !

                  ! nr = nbr_list(nbr_nr,l)
                  ! nr is nbr_list_inv(nbl,nr)-th nbr of l !
                  ! or
                  ! l = nbr_list(nbl,nr)
                  ! nr = nbr_list(nbr_list_inv(nbl,nr),l)
                  
       enddo ! do nr = 1,natoms !

!$OMP END PARALLEL DO

      return
      end subroutine ! get_nbrs_inv !
!
!--------------------------------------------------------------------
!
      function nodeUp_of(node)
       use sim_box

       node_Z = node/nodes_on_Y
       node_Y = mod(node,nodes_on_Y)
       nodeUp_of = mod(node_Z+1,nodes_on_Z)*nodes_on_Y + node_Y

      return
      end  ! nodeUp_of(node) !
!
!--------------------------------------------------------------------
!
      function nodeDown_of(node)
       use sim_box

       node_Z = node/nodes_on_Y
       node_Y = mod(node,nodes_on_Y)
       nodeDown_of = 
     = mod(node_Z-1+nodes_on_Z,nodes_on_Z)*nodes_on_Y + node_Y

      return
      end  ! nodeDown_of(node) !
!
!--------------------------------------------------------------------
!
      function nodeLeft_of(node)
       use sim_box

       node_Z = node/nodes_on_Y
       node_Y = mod(node,nodes_on_Y)
       nodeLeft_of = 
     = node_Z*nodes_on_Y + mod(node_Y-1+nodes_on_Y,nodes_on_Y)

      return
      end  ! nodeLeft_of(node) !
!
!--------------------------------------------------------------------
!
      function nodeRight_of(node)
       use sim_box

       node_Z = node/nodes_on_Y
       node_Y = mod(node,nodes_on_Y)
       nodeRight_of = 
     = node_Z*nodes_on_Y + mod(node_Y+1,nodes_on_Y)

      return
      end  ! nodeRight_of(node) !
!
!-------------------------------------------------------------------
!
      subroutine link_cell_setup

      use sim_box

      implicit double precision (a-h,o-z)
!
! subroutine to set up a cell structure in which to assign atoms
!
      ncell_per_node_old = ncell_per_node

      ! Parallel values               ! 1D,   Serial values        !
      mynodZ = 0
      mynodY = 0

      ncZ_per_node = nnz
      lZstart = 0
      lZend = nnz-1 
      iZ_shift = nnz
      ncY_per_node = nny
      lYstart = 0
      lYend = nny-1
      iY_shift = nny
      nXYlayer = nnx*nny
      ncell_per_node = nnx*nny*nnz

!     write(6,10)nodes_on_Y,ncell_per_node
! 10  format('link_cell_setup: nodes_on_Y=',i2,' ncell_per_node=',i5)

      cellix=real(nnx)
      celliy=real(nny)
      celliz=real(nnz)

      ncell_per_node = (ncell_per_node/8 + 1)*8

! *** Neighbor nodes index ***

      !  kp1YZ | kp1Z  | kp1ZY !
      ! ---------------------- !
      !   km1Y | mynod | kp1Y  !
      ! ---------------------- !
      !  km1YZ | km1Z  | km1ZY !

      kp1Z = 0  ! Serial mode !
      km1Z = 0
      kp1Y = 0
      km1Y = 0

      kp1YZ = 0
      kp1ZY = 0
      km1YZ = 0
      km1ZY = 0

      if(ncell_per_node.gt.ncell_per_node_old) then
       if(NOT(allocated(id_of_cell))) then
        call alloc_cells(ierror) ! alloc_ates cells in aladyn_pi_mods !
        call error_check(ierror,'ERROR in alloc_cells...')
       endif 
      endif
!
! Collect cell ids and indices for use in get_ANN_nbrs, get_PINN_nbrs
!
      k = 0

      do iz=lZstart,lZend
       izz = iz+iZ_shift
       iz_nr = mod(izz,nnz)

       if(nodes_on_Y.eq.1) then  ! 1D node topology !
        do iy=0,nny-1
         iyy = iy+nny
         iy_nr = mod(iyy,nny)
        do ix=0,nnx-1
         k = k + 1
         ixx = ix+nnx
         id_of_cell(k) = iz_nr*nXYlayer + iy_nr*nnx + mod(ixx,nnx)
        enddo
        enddo
       else  ! 2D node topology !
        do iy=lYstart,lYend
         iyy = iy+iY_shift
         iy_nr = mod(iyy,nny)
        do ix=0,nnx-1
         k = k + 1
         ixx = ix+nnx
         id_of_cell(k) = iz_nr*nXYlayer + iy_nr*nnx + mod(ixx,nnx)
        enddo
        enddo
       endif  ! if(nodes_on_Y.eq.1)... !

      enddo

      ncells_all = k

      return   ! link_cell_setup !
      end
!
!**********************************************************************
!   ***  report istantaneous properties  ***
!**********************************************************************
!
      subroutine report(jstep)

      use pot_module
      use MD

      implicit double precision (a-h,o-z)

      integer, intent(in) :: jstep
!
!   ***  temperature and energy
!
      epot=PotEnrg_glb/natoms
      call get_T
      etott=epot+Ek_sys
 
      write(6,942) jstep,real_time,epot,Ek_sys,etott,T_sys

 942  format(i8,' t=',f10.2,' ps, Ep=',f12.8,' + Ek=',f12.8,
     1 ' = Etot=',f12.8,' eV/atom,   Tsys=',f8.2,' K')

      return
      end         ! report !
!
! ---------------------------------------------------------------------
!
      subroutine PROGRAM_END(ierr)

      use sys_OMP
      use sys_ACC
      use sim_box
      use pot_module
      use atoms
      use ANN
      use PINN

      implicit double precision (a-h,o-z)

      integer, intent(in) :: ierr

      close (unit = 5)

      call deall_atoms_sys(ierror)
      call deall_atoms_MD(ierror)

      call deall_types_ANN(ierror)
      call deall_atoms_ANN(ierror)
      call deall_atoms_PINN(ierror)

      call deall_buffers(ierror)
      call deall_pot_types(ierror)
      call deall_cells(ierror)

!$ACC  EXIT DATA
!$ACC& DELETE(Nodes_of_layer, r0_value, base_pot_param)
!$ACC& DELETE(W1_ann, W2_ann, W3_ann, B1_ann, B2_ann, B3_ann)

      if(ierr.ne.0) stop   ! PROGRAM_END was called due to an error !
                           ! STOP the execution, otherwise continue !
      return
      end       ! PROGRAM_END !
!
! ---------------------------------------------------------------------
!
