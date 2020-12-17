!------------------------------------------------------------------
! 11-20-2020
!
! Input - Output module subroutines for aladyn_pi.f code
!
! Vesselin Yamakov
! National Institute of Aerospace
! 100 Exploration Way,
! Hampton, VA 23666
! phone: (757)-864-2850
! fax:   (757)-864-8911
! e-mail: yamakov@nianet.org
!------------------------------------------------------------------
!
      MODULE IO

      use sim_box
      use pot_module
      use atoms
      use ANN
      !save

      CONTAINS
!
!---------------------------------------------------------------------
!  Command Options:
!  Version type:      aladyn_pi -v 1  (default 0)
!  Number MD steps:   aladyn_pi -n 10 (default 10)
!  Measure step:      aladyn_pi -m 5  (default 1)
!---------------------------------------------------------------------
!
      subroutine read_Args

      integer :: N_args

      character*8 str_opt,str_num
      character*4 chem

      ! Default values: !
       iver = 0
       nstep = 10
       measure_step = 1
       dt_step = 0.001d0  ! 1 fs = 0.001 ps !
       start_time = 0.d0
       T_set = 100.d0 ! K !

       k = 1
       N_args = IARGC()

       do while(k.le.N_args)
        call GETARG(k,str_opt); k = k+1

        if((str_opt(1:2).eq.'-v').OR.(str_opt(1:2).eq.'-V')) then
          call GETARG(k,str_num); k = k+1
          read(str_num,*) iver
          write(6,31)iver
        elseif((str_opt(1:2).eq.'-n').OR.(str_opt(1:2).eq.'-N')) then
          call GETARG(k,str_num); k = k+1
          read(str_num,*) nstep
          write(6,32) nstep
        elseif((str_opt(1:2).eq.'-m').OR.(str_opt(1:2).eq.'-M')) then
          call GETARG(k,str_num); k = k+1
          read(str_num,*) measure_step 
          write(6,33) measure_step
        endif ! if((str_opt(1:2)... !

       enddo  ! do while(k.le.N_args) !

       nodes_on_Y = 1
       nodes_on_Z = 1

  31   format(' Version:',i3)
  32   format(/,' Executing:    ',i8,' MD steps')
  33   format(' Measure at each ',i6,' MD steps')

      return
      end subroutine    ! read_Args !
!
!---------------------------------------------------------------------
!   ***  read the input parameters for md calculation  ***************
!        Called only once in Main after read_pot to read aladyn_pi.com
!        first initialization lines.
!---------------------------------------------------------------------
!
      subroutine read_com

      use constants

      implicit double precision (a-h,o-z)

      character*200 LINE0,LINE
      character*1 LeadChar
      character*4 chem
      character*8 str_format_inp,str_format_out, Astring

  10  format(A200)
  12  format(A2)
  16  format('LINE: ',A200)
  50  format(A)

      no_line = 0
      new_mu0 = 0
      new_alpha0 = 0
      iensemble = 0
      itriclinic = 0 ! assume ORTHORHOMBIC system !
      at_vol_inp = -1.d0 ! no preset atomic volume !
      iHEAD_WRITE = 1  ! write a header line in the *.dat file !
      iCrack = 0       ! No crack !

      ierror = 0
      start_time = 0.0

      do i=0,iatom_types
C      iZ_elem_in_com(i) = 13  ! Al !
       iZ_elem_in_com(i) = 14  ! Si !
      enddo

! --- Some Initializations --- 

      real_time = start_time    ! [fms] !
      T_sys = T_set

      r_cut_short = r_cut_off   ! pot. cut off !
      if(iPOT_file_ver.eq.6) r_cut_off = 1.5d0*r_cut_off

      HowBig = r_cut_off
      size = r_cut_off
      !write(6,*)'read_com: ',r_cut_off,HowBig

      r2_cut_off = r_cut_off**2  ! [Ang^2] !
      CutOffRadiusSquared = r2_cut_off
      r2_cut_short = r_cut_short**2

      return
      end subroutine           ! read_com !
!
!--------------------------------------------------------------------
! This subroutine reads in data from files set in pot.dat file 
! describing a trained Neural Network for a specific potential format
! iPOT_file_ver = 5: straight ANN; 6: BOP ANN;
!
! FORMAT:
! iflag_ann,range_min_ann,Rc_ann,d_ann,Gauss_ann
! n_set_ann, (r0_value(i), i=1,n_set_ann)
! w1(1,1), w1(1,2),…w1(1,20), w1(2,1),w1(2,2),…w1(60,20)
! b1(1), b1(2),…b1(20)
! w2(1,1), w2(1,2),…w2(1,20), w2(2,1),w2(2,2),…w2(20,20)
! b2(1), b2(2),…b2(20)
! w3(1,1), w3(2,1), w3(3,1),…w3(20,1)
! b3(1)
!--------------------------------------------------------------------
!
      subroutine input_pot_ANN(ierror)

      implicit double precision (a-h,o-z)

      integer, intent(out) :: ierror
      integer :: ierr
      double precision, dimension(:,:), allocatable :: WT_ann

      character*200 LINE, err_msg
      character*1024 LINE_base_pot    ! Long for multi-comp. potential !
      character*4 elem_symb_pot(112) ! maximum numb. elements possible !
      double precision gram_mol_pot(112)

      character*32 LINE32

      ierr = 0; net_in = 0; nBOP_params = 0
      net_atom_types = 1

      Nodes_of_layer(:) = 0
      nunit=40  ! PINN.dat file !

      err_msg = ' '
      ierror = 0
      nparam_base = 0

      open (nunit,file=filename(1),status='old')

      write(6,*)' '
      write(6,*)' READING pot file: ',filename(1),'...'

! *** Start reading the Artificial Neural Network file *.ann ***

      read(nunit,fmt=50,IOSTAT=ierr) LINE  ! 1st Line pot file version !
      read(LINE,*) iPOT_file_ver, gi_shift, nAct_func_type

      ActFunc_shift = 0.d0
      select case(nAct_func_type)   
       case(0) 
        ActFunc_shift = 0.d0
       case(1) 
        ActFunc_shift = -0.5d0 ! f(x) = 1/(1+exp(-x)) + ActFunc_shift !
      end select

      read(nunit,fmt=50,IOSTAT=ierr) LINE  ! 2nd Line !
      read(LINE,*) net_atom_types          ! number of chemical elem. !

      do i=1,net_atom_types
       read(nunit,*) elem_symb_pot(i), gram_mol_pot(i)
      enddo

      write(6,*)' '
      write(6,*)'Elements in ANN potential file:'
      write(6,*) (elem_symb_pot(n),n=1,net_atom_types)

      ipot(1) = 1
      write(6,*)' '

      read(nunit,fmt=50,IOSTAT=ierr) LINE  ! 4:Interaction Range Line !
  50  format(A)
      read(LINE,*) iflag_ann,range_min_ann,Rc_ann,d_ann,Gauss_ann

      read(nunit,fmt=50,IOSTAT=ierr) LINE32  ! 5:Legendre Pol. Line !

      read(nunit,fmt=50,IOSTAT=ierr) LINE  ! 6:Gaussians Line !
      if(ierr.eq.0) then
       read(LINE,*) n_set_ann
       if(n_set_ann.gt.0) then
        allocate(r0_value(n_set_ann), stat=ierr)
        if(ierr.eq.0) read(LINE,*) nn, (r0_value(i), i=1,n_set_ann)
       else
        ierr = 1; 
      write(6,*)'ERROR: No Gaussian positions in line 2 in ',filename(1)
       endif

       if(ierr.eq.0) then
      write(6,10)n_set_ann,(r0_value(i), i=1,n_set_ann) 
       endif
      endif ! if(ierr.eq.0)... 6:Gaussians Line !!

      read(nunit,fmt=50,IOSTAT=ierr) LINE  ! 7: ANN_nodes Line !
      if(ierr.eq.0) then
       read(LINE,*) net_layers
       if((0.lt.net_layers).AND.(net_layers.le.Max_net_layers)) then
        read(LINE,*)nn,(Nodes_of_layer(i),i=1,net_layers)
        net_in = Nodes_of_layer(1)
        nBOP_params = Nodes_of_layer(net_layers)
        net_in_check = 5*n_set_ann*iatom_types 

        if(net_in.ne.net_in_check) then
         ierr = 1
      write(6,13) net_in, net_in_check, n_set_ann, iatom_types
        endif
       else
        ierr = 1; 
      write(6,*)'ERROR: Incorect Net layers in line 3 in ',filename(1)
      write(6,22) net_layers, Max_net_layers
  22  format('Number of Net layers =',i3,' must be between 1 and ',i2)
       endif ! if((0.lt.net_layers)... !

      if(nBOP_params.eq.1) then
       iPOT_file_ver = 5
      else
       iPOT_file_ver = 6
      endif

      select case(iPOT_file_ver)
       case(5); write(6,*)'Straight ANN potential'
       case(6); write(6,*)'PINN: Physically Informed NN potential'
      end select

      if(ierr.eq.0) then
      write(6,12)net_layers,(Nodes_of_layer(i), i=1,net_layers) 
!     write(6,20) net_atom_types, net_in, nBOP_params, Rc_ann, d_ann

       nparam_base = Nodes_of_layer(net_layers) ! = last ANN layer !
       allocate(base_pot_param(nparam_base), stat=ierrb)
       base_pot_param(:) = 0.d0

      iBOP_base_exist = 0                   ! 8:Base BOP param. Line !
      read(nunit,fmt=50,IOSTAT=ierr) LINE_base_pot  ! < 1024 chars !
      read(LINE_base_pot,*) iBOP_base_exist

       if(iBOP_base_exist.gt.0) then   ! Read base potential param. !
        if(ierrb.eq.0) then
         read(LINE_base_pot,*,IOSTAT=ieol)
     1   ibase,(base_pot_param(i),i=1,nparam_base)
         if(ieol.eq.0) then
          write(6,23)(base_pot_param(i),i=1,min(nparam_base,8))
  23      format('Base Pot. Parameters: ',8(f12.8,1x),'...')
         else
          write(6,*)'ERROR reading base_pot_param on node 0...'
          write(6,24) nparam_base
  24      format('Number of base pot. parameters is less than the ANN',
     1    ' output (',i3,')...')
          ierr = 1
         endif
        else
         write(6,*)'ERROR allocating base_pot_param on node 0...'
         ierr = 1
        endif
       endif ! if(iBOP_base_exist.gt.0)... !
      endif ! if(ierr.eq.0)... !

      endif ! if(ierr.eq.0)... 8: Line !!

      call error_check(ierr,'ERROR reading ANN file in input_pot_ANN')

      call error_check(net_atom_types-iatom_types,
     1 'ERROR: Elements in ANN file and pot.dat do not match!')

      allocate(r0G_value(n_set_ann), stat=ierr)
      call error_check(ierr,'ERROR allocate r0G_value in input_pot_ANN')
      r0G_value(:) = r0_value(:)/Gauss_ann

      call alloc_types_ANN(ierr)
      call error_check(ierr,'ERROR alloc_types_ANN in input_pot_ANN')

      ierr = 0
      ww = 0.d0; bb = 0.d0

! --- Read Input Layer Parameters for atom of type itype ---
        Ncolumns = Nodes_of_layer(1)  ! 60 !
        Nraws = Nodes_of_layer(2)     ! 20 !
        do icol = 1, Ncolumns ! 1.. 60 !
         do iraw = 1, Nraws ! 1.. 20: w(1,1), w(1,2), w(1,3)... !
          read(nunit,fmt=50,IOSTAT=ierr) LINE  ! ANN Line !
          if(ierr.eq.0) read(LINE,*) W1_ann(icol,iraw)
         enddo
        enddo
        do iraw = 1, Nraws
         read(nunit,fmt=50,IOSTAT=ierr) LINE  ! ANN Line !
         if(ierr.eq.0) read(LINE,*) B1_ann(iraw)
        enddo

! --- Read Hidden Layers Parameters for atom of type itype ---
       do layer = 2,net_layers-2
        Ncolumns = Nodes_of_layer(layer)  ! 20 !
        Nraws = Nodes_of_layer(layer+1)       !  8 !
        do icol = 1, Ncolumns ! 1.. 20 !
         do iraw = 1, Nraws ! 1.. 20: w(1,1), w(1,2), w(1,3)... !
          read(nunit,fmt=50,IOSTAT=ierr) LINE  ! ANN Line !
          if(ierr.eq.0) read(LINE,*)ww
          W2_ann(icol,iraw,layer) = ww
         enddo
        enddo
        do iraw = 1, Nraws
         read(nunit,fmt=50,IOSTAT=ierr) LINE  ! ANN Line !
         if(ierr.eq.0) read(LINE,*)bb
         B2_ann(iraw,layer) = bb
        enddo
       enddo ! do layer = 2,net_layers-2

! --- Read Output Layer Parameters for atom of type itype ---
        Ncolumns = Nodes_of_layer(net_layers-1)  ! 20 !
        Nraws = Nodes_of_layer(net_layers)       !  8 !
        do icol = 1, Ncolumns ! 1.. 20 !
         do iraw = 1, Nraws ! 1.. 20: w(1,1), w(1,2), w(1,3)... !
          read(nunit,fmt=50,IOSTAT=ierr) LINE  ! ANN Line !
          if(ierr.eq.0) read(LINE,*) W3_ann(icol,iraw)
         enddo ! do iraw = 1, Nodes_of_layer(layer+1) !
        enddo ! do icol = 1, Nodes_of_layer(layer) !
        do iraw = 1, Nraws
         read(nunit,fmt=50,IOSTAT=ierr) LINE  ! ANN Line !
         if(ierr.eq.0) read(LINE,*) B3_ann(iraw)
        enddo ! do iraw = 1, Nodes_of_layer(layer+1) !

       close (nunit)

      call error_check(ierr,'ERROR reading W_ann and B_ann arrays')

! Swap index order from (n_set,l,...) to (l,n_set,...) of
! the FIRST layer weights - FORTRAN order.

      Ncolumns = Nodes_of_layer(1)  ! 60 !
      Nraws = Nodes_of_layer(2)     ! 16 !

      allocate(WT_ann(Ncolumns,Nraws), stat=ierr)

      WT_ann(:,:) = W1_ann(:,:)

         do l=0,4
          do n_set=1,n_set_ann
           ind = l*n_set_ann + n_set
           icol = (n_set-1)*5 + l+1
           W1_ann(icol,:) = WT_ann(ind,:)
          enddo ! n_set !
         enddo ! l=0,4 !

      if (allocated(WT_ann)) deallocate(WT_ann, stat=ierr)
     
      ierror = ierr

  10  format('Gaussian positions:',i3,':',32f6.2) 
  12  format('NN layers=:',i3,' of nodes:',10i5) 
  13   format(/, 
     1 'ERROR: Inconsistency b/n the number of input net nodes:',i4,
     1 ' and Structure Parameters:',i4,' = 5 *',i4,' * ',i2)
  15  format(/,'ERROR: Element ',A2,' from pot.dat not in ',A,/)
  20  format('Atom types in neural net:',i2,/,
     1 'Structure parameters (net_in) :',i4,/,
     1 'Potential parameters (nBOP_params):',i4,/,
     1 'Rc (R_cut-off) =',f8.4,',  d (in Fcut) =',f8.4,/)

      return
      end subroutine    ! input_pot_ANN !
!
!--------------------------------------------------------------------
! Called from calc_pot_param() in aladyn_pi.f
!--------------------------------------------------------------------
!
      subroutine init_param_ANN

      implicit double precision (a-h,o-z)

      ! Get the maximum cut-off distance, r_cut_off and r2_cut_off !
      ! and the maximum overlapping distance, rin                  !

      r_cut_off = 0.0d0

      mG_dim = (iatom_types*(iatom_types+1)*n_set_ann*5)/2
      max_tri_index = (iatom_types*(iatom_types+1))/2 ! Upper triang. !

      if(mG_dim.ne.Nodes_of_layer(1)) then
       call error_check(1,
     1 'ERROR dim. of Gis not equal to Nodes_of_layer(1)...')
      endif

      d4_ann = d_ann**4
      r_cut_off = Rc_ann
      r_max_pot = r_cut_off + 1.0 ! Add additional 1 Ang to pot array !
                          ! for consistency with tabulated potentials !
      rindr = 0.5d0       ! [Ang] !

!$ACC ENTER DATA   ! EXIT DATA at subroutine PROGRAM_END() !
!$ACC& COPYIN(Nodes_of_layer, r0_value, W1_ann, W2_ann, W3_ann)
!$ACC& COPYIN(B1_ann, B2_ann, B3_ann, base_pot_param)

      return
      end subroutine      ! init_param_ANN  !
!
!---------------------------------------------------------------------
! Read input structure files in plt or lammps format
! Note: read_structure must be called AFTER read_pot and read_com
!---------------------------------------------------------------------
!
      subroutine read_structure

      implicit double precision (a-h,o-z)

      select case(INP_STR_TYPE)
       case default
         call read_structure_plt     ! plt type !
      end select

      write(6,*)' '
      write(6,*)'h(i,j) matrix:'
      write(6,80) h(1,1),h(1,2),h(1,3)
      write(6,80) h(2,1),h(2,2),h(2,3)
      write(6,80) h(3,1),h(3,2),h(3,3)
      write(6,82) natoms

  80  format(3(f14.8,2x))
  82  format('Crystal structure has ',i9,' atoms',/)

      call structure_chem

      call matinv(h,hi,dh) ! h*hi=I
      hi11 = hi(1,1); hi12 = hi(1,2); hi13 = hi(1,3)
      hi22 = hi(2,2); hi23 = hi(2,3); hi33 = hi(3,3)

      do n=1,natoms
       sx(n) = hi11*rx(n) + hi12*ry(n) + hi13*rz(n)
       sy(n) = hi22*ry(n) + hi23*rz(n)
       sz(n) = hi33*rz(n)
      enddo

      return
      end subroutine           ! read_structure !
!
!---------------------------------------------------------------------
!  ***  read the input atomic structure for MD simulation          ***
!       using plt format from "sold"   
!---------------------------------------------------------------------
!
      subroutine read_structure_plt

      implicit double precision (a-h,o-z)

      character*40 file_in

      character*82 str82
      character*80 LINE
      character*85 LINE2
      character*2 DUM2

      rsmax0 = 1.0

      open (unit=11,status='old',file='structure.plt')

      itime=int(start_time)
!     write(6,15) itime
 15   format(' Start time:',i9,' [MCS]')

      write(6,*) 'Reading structure.plt ...'

      read(11,*)DUM2, xmn,  ymn,  zmn
      read(11,*)DUM2, xmx,  ymx,  zmx
      read(11,*)DUM2, xtmn, ytmn, ztmn
      read(11,*)DUM2, xtmx, ytmx, ztmx  
      read(11,*)DUM2, nbas,natoms_tot,natoms_buf,natoms_free
                       ! total, buffer, free atoms !
      read(11,*)DUM2, r_plt, mdx, mdy, mdz
      read(11,*)DUM2, ibcx, ibcy, ibcz
      read(11,*)DUM2, ipo, ipl
      read(11,*)DUM2, PotEnrg_atm 

! Check if the plt file has a dof column
        read (11,FMT=10) LINE
 10     format(A80)
        write(unit=LINE2,fmt=20) LINE,' -1' ! add "-1" to line end !
 20     format(A,A)
!       write(6,*) 'LINE2:',LINE2
        read (LINE2,*) id,xk,yk,zk,ktype,kdof
!       write(6,*) id,xk,yk,zk,ktype,kdof

        if(kdof.eq.-1) then
         kdof_exist = 0  ! old plt file (no constrains) !
        else
         kdof_exist = 1  ! new plt file with constrains !
        endif
        rewind 11

      read(11,*)DUM2, xmn,  ymn,  zmn
      read(11,*)DUM2, xmx,  ymx,  zmx
      read(11,*)DUM2, xtmn, ytmn, ztmn
      read(11,*)DUM2, xtmx, ytmx, ztmx  
      read(11,*)DUM2, nbas,natoms_tot,natoms_buf,natoms_free
                       ! total, buffer, free atoms !
      read(11,*)DUM2, r_plt, mdx, mdy, mdz
      read(11,*)DUM2, ibcx, ibcy, ibcz
      read(11,*)DUM2, ipo, ipl
      read(11,*)DUM2, PotEnrg_atm 

!     write(6,*)'kdof_exist=',kdof_exist

! First spread the atoms uniformly to all nodes

      natoms = natoms_tot 

      h(1,1) = (xtmx - xtmn)  ! system size from second 2 lines !
      h(2,2) = (ytmx - ytmn)
      h(3,3) = (ztmx - ztmn)
      !write(6,*)'read_structure_plt: ',h(1,1),h(2,2),h(3,3)

      h(1,2) = 0.0d0; h(1,3) = 0.0d0; h(2,3) = 0.0d0  ! ORT structure !
      h(2,1) = 0.0d0; h(3,1) = 0.0d0; h(3,2) = 0.0d0

      call matinv(h,hi,dh)

      write(6,*)'Input file Pot. Enrg=',PotEnrg_atm
      write(6,*)' '

      nxold=0;nyold=0;nzold=0  ! first node_config call !
      call node_config(nflag)  ! gets ncell, ncell_per_node  !
                               ! natoms_alloc, call alloc_atoms !
      ndof_flag = 0
 
        do i=1,natoms
         idof = 0
         if(kdof_exist.ne.0) then
          read(11,*) id,xk,yk,zk,ktype,idof
          ndof_flag = IOR(ndof_flag,idof) ! collect DOF constraints !
         else
          read(11,*) id,xk,yk,zk,ktype
         endif

         ident(i) = id
         rx(i) = xk
         ry(i) = yk
         rz(i) = zk
         ntype(i) = ktype
	enddo

       close(unit=11)

      return     
      end subroutine       ! read_structure_plt !
!
!---------------------------------------------------------------------
! Establish the initial chemical composition from structure file
! Called by read_structure (only once)
!---------------------------------------------------------------------
!
      subroutine structure_chem

      implicit double precision (a-h,o-z)

      ierror = 0

! *** Collect atom types ...
      natoms_of_type(:) = 0

      do n=1,natoms
       ntp = ntype(n)   ! corresponds to ann.dat type order !
       if((ntp.gt.iatom_types).OR.(ntp.eq.0)) then
        ierror = 1
        write(6,13)n,ident(n),ntp
        exit
       else
        natoms_of_type(ntp) = natoms_of_type(ntp) + 1
       endif
      enddo

      if(ierror.ne.0) then
       call PROGRAM_END(1)
      endif

       write(6,14) ncell_per_node, ncell, cell_volume, natoms_alloc,
     1 natoms_per_cell, natoms_per_cell3, nbrs_per_atom

 14   format(/,'ncell_per_node=',i5,'; ncell=',i6,
     1 '; cell_volume=',f10.2,' [Ang^3];',/,
     2 '    Atoms allocated per node:',i8,/,
     2 '    Atoms allocated per cell:',i8,/,
     2 '    Atoms allocated per 3x3x3 cells:',i8,/,
     2 'Neighbors allocated per atom:',i8)
 15   format('    Close neighbors per atom:',i8)

 13   format(/,'ERROR: atom n=',i8,' of id=',i8,
     1 ' is of unknown type:',i3,/)
 83   format(i2,': (',A2,'):  Atom %:',f9.4)

      return     
      end subroutine       ! structure_chem !
!
!--------------------------------------------------------------------
! This subroutine reads pot.dat file
!--------------------------------------------------------------------
!
      subroutine read_pot_dat

      use constants

      implicit double precision (a-h,o-z)

      character*200 filename0

      ierror = 0

      r_ADP_cut = 0.0d0
      iatom_types = 1   ! Number of chemical elements !

* ==== Reading chemical species and filenames of the potential ======

      ifile_numbs = iatom_types*(iatom_types+1)/2
      ipair_types = iatom_types**2

      call alloc_pot_types(ierror) ! in pot_module !
      ! allocates arrays that are common to ALL pot file types and
      ! formats 

      call error_check(ierror,
     1 'alloc_pot_types error in read_pot_dat...')

      ielement(0) = 0
!     elem_symb(1) = 'Al'
!     gram_mol(1) = 26.982
      elem_symb(1) = 'Si'
      gram_mol(1) = 28.085500
      ielem = numb_elem_Z(elem_symb(1))
      ielement(1) = ielem      ! Z - number !
      amass(1) = gram_mol(1)*atu  ! convert to atomic units !
      ! atu = 100.0d0/(cavog*evjoul) = 0.010364188 eV.ps^2/nm^2 !

!     filename(1) = './ANN.dat'
      filename(1) = './PINN.dat'

      write(6,103)
      do i=1,iatom_types
       write(6,104)i,elem_symb(i),ielement(i),gram_mol(i),iPOT_file_ver
      enddo

      nelem_in_com = iatom_types

      MC_rank = 3

 103  format(/' CHEMICAL ELEMENTS:','   TYPE   ELEMENT',
     1 '    Z    Atomic Mass  POT_func_type')
 104  format(20x,i5,6x,a4,3x,i3,2x,f10.4,6x,i4)
 110  format(A70)

      return     
      end subroutine       ! read_pot_dat !
!
!---------------------------------------------------------------------
!  Writes structure output in plt format
!---------------------------------------------------------------------
!
      subroutine write_structure_plt

      implicit double precision (a-h,o-z)

      character fname*80, ffname*200
      double precision h_out(3,3)

      isize = 7
      itime=int(real_time)

      dt1 = 1.d0/dt_step ! 1/[ps] !

        write(fname,80) 'structure.',itime,'.plt'

        n = len(fname)
        k = 0
        do i=1,n
         if(fname(i:i).eq.' ') then 
          k = k+1
         else
          fname(i-k:i-k) = fname(i:i)
         endif
        enddo
        do i=n-k+1,n 
         fname(i:i)=' '
        enddo

        write(ffname,90) file_path,fname
 15     open (unit=22,status='unknown',file=ffname) 

 70    format(A,A4)
 80    format(A10,I8.8,A4)
 90    format(A,A)

      xtmx = h(1,1)*0.5d0 ! half the size in X in [Ang.] !
      xtmn = -xtmx
      ytmx = h(2,2)*0.5d0 ! half the size in Y in [Ang.] !
      ytmn = -ytmx
      ztmx = h(3,3)*0.5d0 ! half the size in Z in [Ang.] !
      ztmn = -ztmx

      write(22,400) xtmn,ytmn,ztmn
      write(22,400) xtmx,ytmx,ztmx
      write(22,400) xtmn,ytmn,ztmn
      write(22,400) xtmx,ytmx,ztmx
      write(22,401) nbas, natoms, natoms_buf, natoms_free
      write(22,402) r_plt, mdx, mdy, mdz
      write(22,401) ibcx, ibcy, ibcz
      write(22,401) ipo, ipl
      write(22,405) PotEnrg_atm, T_sys

      do kk=1,natoms

       Cx = h(1,1)*sx(kk) + h(1,2)*sy(kk) + h(1,3)*sz(kk)
       Cy = h(2,2)*sy(kk) + h(2,3)*sz(kk)
       Cz = h(3,3)*sz(kk) 

       ntp = ntype(kk)  ! write type as in pot.dat file !
       write(22,127) ident(kk),Cx,Cy,Cz,ntp, 0

      enddo ! do kk=1,natoms !

      write(22,130) 0 ! 0 - no velocities; 1 - vel. !
      
      close(unit=22)

 127  format(I8,3E18.10, 2(1x,i2))
 128  format(I8,3E18.10,1x,i2,1x,f13.8)
 130  format(i2)
 137  format(I8,3E18.10)
 140  format(2f4.1)
 400  format('#',4E18.10)
 401  format('#',4I8)
 402  format('#',E18.8,3I8)
 405  format('#',E18.10,'  ',f8.1)

      return
      end subroutine     ! write_structure_plt !
!
!---------------------------------------------------------------------
!
      subroutine alloc_atoms

      use sys_OMP
      !use sys_CUDAFOR
      !use sys_ACC
      use ANN
      use PINN

      write(6,*)'ALLOCATE_ATOMS: ',natoms_alloc

      call alloc_atoms_sys(ierror)
      call error_check(ierror,
     1 'Memory alloc_atoms_sys error in alloc_atoms...')

      if(ierror.ne.0) return

       call alloc_atoms_ANN(ierror)
       call error_check(ierror,
     1 'Memory alloc_atoms_ANN error in alloc_atoms...')

       if(iPOT_file_ver.eq.6) then
        call alloc_atoms_PINN(ierror)
        call error_check(ierror,
     1 'Memory alloc_atoms_PINN error in alloc_atoms...')
       endif

      call alloc_atoms_MD(ierror) ! alloc x1,.. x5 !

      nbufsize = 16*natoms_alloc  ! max isize in finder_MD !
      if(ierror.eq.0) call alloc_buffers(ierror)

      call error_check(ierror,
     1 'MEMORY alloc_buffers error in alloc_atoms...')

      return
      end subroutine    ! alloc_atoms !
!
!--------------------------------------------------------------------
!
      subroutine finder

      implicit double precision (a-h,o-z)

      natoms_in_cell(0:ncell_per_node) = 0  ! speedup: 1.450 !

      cellix=dfloat(nnx)
      celliy=dfloat(nny)
      celliz=dfloat(nnz)

      do n=1,natoms          ! 104 vs 51: speedup: 2.010 !
        sxn = sx(n)
        syn = sy(n)
        szn = sz(n)

        if (sxn.ge. 0.5d0) then   ! make periodic along X !
            sxn=sxn-1.0d0         ! 2xfaster than sxn=sxn-dnint(sxn) !
        elseif (sxn.lt.-0.5d0) then
            sxn=sxn+1.0d0
        endif
        if (syn.ge. 0.5d0) then     ! make periodic along Y !
            syn=syn-1.0d0 
        elseif (syn.lt.-0.5d0) then
            syn=syn+1.0d0
        endif
        if (szn.ge. 0.5d0) then     ! make periodic along Z !
            szn=szn-1.0d0 
        elseif (szn.lt.-0.5d0) then
            szn=szn+1.0d0
        endif

        sx(n) = sxn
        sy(n) = syn
        sz(n) = szn
      enddo  ! do n=1,natoms !
!
! *** now update the link-cell arrays
!
      ierror = 0
      MC_rank_XY = MC_rank_X * MC_rank_Y

      do i=1,natoms

       ix=int((sx(i)+0.5d0)*cellix)
       iy=int((sy(i)+0.5d0)*celliy)
       iz=int((sz(i)+0.5d0)*celliz)

!  put atom "i" in cell "(ix,iy,iz)"

       icell = mod(iz+iZ_shift,nnz)*nXYlayer +
     +         mod(iy+iY_shift,nny)*nnx      + 
     +         mod(ix+nnx,nnx)

       ncell_of_atom(i) = icell

       if(natoms_in_cell(icell).lt.natoms_per_cell) then
        natoms_in_cell(icell) = natoms_in_cell(icell) + 1
        n_in_cell(natoms_in_cell(icell),icell) = i
       else
        ierror = 1
       endif

      enddo ! do i=1,natoms !
!     natomsInCell_d = natoms_in_cell
!     NumInCell_d = n_in_cell

       call error_check(ierror,
     1 'ERROR: FINDER: Too many atoms per cell...')

      return
      end subroutine  ! finder !
!
! ---------------------------------------------------------------------
!
      END MODULE  ! IO !
!
!=====================================================================
