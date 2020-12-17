!------------------------------------------------------------------
! 11-20-2020
!
! General Module Unit for aladyn_pi.f code
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
      MODULE constants
      save
!
! *** Set some physics constants as parameters here ***
!
       double precision, parameter :: Boltz_Kb = 8.6173324d-5 ! eV/K !
             ! = 8.6173324d-5*1.602189d-19 (J) = 1.380659518d-23 ! J/K !
       double precision, parameter :: comp=0.0001
       double precision, parameter :: erg=0.1602189d-11
       double precision, parameter :: hp=0.0041354d0
       double precision, parameter :: cavog=602.2045d0
       double precision, parameter :: evdyne=1.602189d0
       double precision, parameter :: evkbar=1.602189d-3
       double precision, parameter :: evjoul=16.02189d0
       double precision, parameter :: pi=4*ATAN(1.0d0)
       double precision, parameter :: pi2=pi/2.d0

       double precision, parameter :: eVAng2GPa=evdyne*100.0d0
      ! conversion from pressure units to kbar (pressure*eVAng2GPa=GPa)

      ! double precision, parameter :: atu = 1.03642696d-8 ! eV(s/m)^2 !
      ! [eV.(ps/A)^2]= [eV.(e-12s/e-10m)^2] = [eV.(s/m)^2]*e-4   !
      ! atu = 1/(N_Avog*e) = 1.03642696d-8 [eV.(s/m)^2] !
      ! atu = 1/(6.022140857e26*1.60217662e-19) = 1/9.6485333e7 = 
      ! = 1.03642696e-8 [eV*s^2/m^2] !
      ! Ek [eV] = 0.5 * m [gmol] * atu * v^2 [m/s]^2 !
      ! atomic mass units in eV.ps^2/nm^2 !
      ! The mass of the atomic elements is expressed in atu units !
      ! For example Al_mass = 26.98 atu= 26.98*1.036427d-8 eV.(ns/nm)^2
      ! For example Al_mass = 26.98 atu= 27*1.0364188   [eV.ps^2/Ang^2]

       double precision, parameter :: atu = 1.03642696d-4 ! eV(ps/A)^2 !
      ! 1 eV(s/m)^2 = eV(10^12ps/10^10A)^2 = 10^4 eV(ps/A)^2 !

       double precision, parameter :: tunit = sqrt(atu)

      END MODULE  ! constants !
!
! =====================================================================
!
      MODULE sim_box
      save
!
! *** Set some working parameters here ***
!
       type node_conf
        character(len=40) :: node_name
        integer :: devicetype 
        integer :: nACC_devices
        integer :: I_have_GPU
        integer :: My_GPU_id
        integer :: MP_procs
        integer :: MP_threads
       end type node_conf

       type group_conf
        character(len=40) :: node_name
        integer :: devicetype 
        integer :: nACC_devices
        integer :: nodes
       end type group_conf

       type(node_conf) :: mynode_conf

       character*2, parameter :: file_path="./"
       integer, parameter :: Max_Lines_INI = 1024
       integer, parameter :: CHUNK = 24

       character*64 :: str_filename

       INTEGER :: MP_threads = 0     ! OMP: Number of OMP threads      !
       INTEGER :: MP_max_threads = 0 ! OMP: MAx. Number of OMP threads !
       INTEGER :: MP_procs = 0   ! OMP: Number of OMP processors !

       integer :: nbufsize = 0

       integer :: MC_rank, MC_rank0, MC_rank_X, MC_rank_Y, MC_rank_Z
       integer, parameter :: MC_rank_max = 19 ! large simple number !

       integer :: nn_count

       integer :: natoms_alloc = 0
       integer :: nbrs_alloc = 0
       integer :: ncell_per_node = 0
       integer :: natoms
       integer :: natoms_tot,natoms_buf,natoms_free
       integer :: nbrs_per_atom, nbrs1_per_atom
       integer :: natoms_per_cell, natoms_per_cell3
       integer :: nbas, ibcx, ibcy, ibcz, ipo, ipl, mdx, mdy, mdz
       integer :: iend,new_mu,new_mu0,new_alpha,new_alpha0
       integer :: nstep,irigid,irigid0,iensemble,measure_step
       integer :: ncpny,ncpnz,ido_finder,nnx,nny,nnz,ncell,
     1 nxold,nyold,nzold,ncells_all
       integer :: mpi_nodes,mynod,ncZ_per_node,ncY_per_node,nXYlayer
       integer :: nodes_on_Y,lYstart,lYend,mynodY,kp1Y,km1Y,iY_shift
       integer :: nodes_on_Z,lZstart,lZend,mynodZ,kp1Z,km1Z,iZ_shift
       integer :: kp1YZ,kp1ZY,km1YZ,km1ZY
       integer :: itriclinic, iHEAD_WRITE = 0

       integer :: ihalt = 0
       integer :: nACC_devices
       integer :: I_have_GPU, My_GPU

       integer :: istep
       integer :: Line_Executed(Max_Lines_INI) 
                  ! = 0: line has not been executed
                  ! = 1: line has been executed
       integer :: n_of_all_moved ! all MC moved atoms up to now !

       integer, dimension(:), allocatable :: ncell_of_atom
       !dir$ attributes align:64 :: ncell_of_atom

       double precision :: BkT ! = 1.0d0/(Boltz_Kb*T_sys) !
       double precision :: real_time, start_time, start_run_time
       double precision :: T_sys, T_set, Ek_sys

       double precision :: h(3,3),hi(3,3)

       ! System shape matrices !
       double precision :: h1(3,3),h2(3,3),h3(3,3),h4(3,3),h5(3,3)

       double precision :: density,a0_nm,r_plt
       double precision :: cellix,celliy,celliz,size,
     1 cell_size_X,cell_size_Y,cell_size_Z,cell_volume

! --- Crack parameters ---

       integer :: iCrack, iLine_crack, iPlane_crack
       double precision :: Crack_tip1, Crack_tip2, Crack_plane
       double precision :: Crack_rAB
!
! Forms a crack as:
!                     A o-----------------------o B
!  where: A = CrackP_A(Crack_tip1,Crack_plane) 
!     and B = CrackP_B(Crack_tip2,Crack_plane)
!  coord_A(x,y) = 
!  coord(iLine_crack=(1,2,3 for x,y,z), iPlane_crack=(1,2,3 for x,y,z))
!

! --- Buffers ---

       double precision, dimension(:), allocatable :: buf, bufa,
     1 bufb, bufc
       !dir$ attributes align:64 :: buf
       !dir$ attributes align:64 :: bufa
       !dir$ attributes align:64 :: bufb
       !dir$ attributes align:64 :: bufc

       integer, dimension(:), allocatable :: nbuf,nbufa
       !dir$ attributes align:64 :: nbuf
       !dir$ attributes align:64 :: nbufa

       integer, dimension(:), allocatable :: kbuf,kbufa,kbufb,kbufc
       !dir$ attributes align:64 :: kbuf
       !dir$ attributes align:64 :: kbufa
       !dir$ attributes align:64 :: kbufb
       !dir$ attributes align:64 :: kbufc

       integer, dimension(:), allocatable :: natoms_in_cell
       !dir$ attributes align:64 :: natoms_in_cell

       integer, dimension(:,:), allocatable :: n_in_cell
       !dir$ attributes align:64 :: n_in_cell

       integer, dimension(:), allocatable :: ibufY1,ibufY2,ibufZ1,ibufZ2
       !dir$ attributes align:64 :: ibufY1
       !dir$ attributes align:64 :: ibufY2
       !dir$ attributes align:64 :: ibufZ1
       !dir$ attributes align:64 :: ibufZ2

       integer, dimension(:), allocatable :: id_of_cell
       !dir$ attributes align:64 :: id_of_cell

       integer :: ipack1Y1,ipack2Y2,jpack1Y1,jpack2Y2
       integer :: ipack1Z1,ipack2Z2,jpack1Z1,jpack2Z2

       double precision :: U(97), C, CD, CM
       double precision :: Us(97), Cs, CDs, CMs
       integer :: I97, J97, I97s, J97s

      CONTAINS
!
!********************************************************
!                                                       $
!  Calculates the argument of a complex number x +iy    $
!                                                       $
!********************************************************
!
        FUNCTION argument(x,y,r)

        double precision x,y,r,xr,yr,phi
        double precision PI = 3.141592654D0

        r = SQRT(x**2 + y**2)
        xr = ABS(x)/r
        yr = ABS(y)/r

        if(xr.gt.0.1D0) then
         phi = ATAN(yr/xr)
        else
         phi = PI/2.0D0 - ATAN(xr/yr)
        endif

        if((x.lt.0.0D0).AND.(y.ge.0.0D0)) phi = PI - phi
        if((x.lt.0.0D0).AND.(y.lt.0.0D0)) phi = PI + phi
        if((x.ge.0.0D0).AND.(y.lt.0.0D0)) phi = 2.0D0*PI - phi

        ARGUMENT = phi

        RETURN
        END FUNCTION
!       
!*************************************************
!						 
!       RANDOM NUMBER GENERATOR			
! 					
!       UNIVERSAL RANDOM NUMBER GENERATOR PROPOSED BY MARSAGLIA
!       AND ZAMAN IN REPORT FSU-SCRI-87-50
!       GENERATES VECTOR 'RVEC' OF LENGTH 'LEN' OF PSEUDORANDOM
!       NUMBERS; THE COMMON BLOCK INCLUDES EVERYTHING NEEDED TO
!       COMPLETELY SPECIFY THE STATE OF THE GENERATOR.
! 
!       Puts LEN random numbers in bufer buf(1.. LEN)
!						 
!*************************************************
!
        SUBROUTINE ranmar(LEN)

        double precision UNI

        DO IVEC=1,LEN
         UNI = U(I97) - U(J97)
         IF(UNI.LT.0.) UNI = UNI + 1.
         U(I97) = UNI
         I97 = I97 - 1
         IF(I97.EQ.0) I97 = 97
         J97 = J97 - 1
         IF(J97.EQ.0) J97 = 97
         C = C - CD
         IF(C.LT.0.) C = C + CM
         UNI = UNI - C
         IF(UNI.LT.0.) UNI = UNI + 1.
         buf(IVEC) = UNI
        ENDDO

        RETURN
        END SUBROUTINE

!*************************************************
!						 
!       RANDOM NUMBER INITIALIZER			
! 					
!       INITIALIZING ROUTINE FOR ranmar. THE INPUT VALUE SHOULD
!       BE IN THE RANGE:   0 <= IJKL <= 900 000 000
!       TO GET THE STANDARD VALUES IN THE MARSAGLIA - ZAMAN PAPER
!       (I=12, J=34, K=56, L=78) PUT  IJKL = 54217137
!						 
!*************************************************

        SUBROUTINE rmarin(IJKL)

        double precision T

        IJ = IJKL / 30082
        KL = IJKL - IJ * 30082
        I  = MOD(IJ/177,177) + 2
        J  = MOD(IJ,177) + 2
        K  = MOD(KL/169,178) + 1
        L  = MOD(KL,169)
!        WRITE(*,*) 'ranmar INITIALIZED: ',IJKL,I,J,K,L
        DO   2    II=1,97
         S = 0.
         T = 0.5
         DO   3    JJ=1,24
          M = MOD(MOD(I*J,179)*K,179)
          I = J
          J = K
          K = M
          L = MOD(53*L+1,169)
          IF(MOD(L*M,64).GE.32) S = S + T
3         T = 0.5 * T
2        U(II) = S
        C = 362436. / 16777216.
        CD = 7654321. / 16777216.
        CM = 16777213. / 16777216.
        I97 = 97
        J97 = 33
        RETURN
        END SUBROUTINE
!
!----------------------------------------------------------------------
!
      subroutine matinv(a,b,c)

      implicit double precision (a-h,o-z)
      double precision a(3,3),b(3,3),c

      d11=a(2,2)*a(3,3)-a(2,3)*a(3,2)
      d22=a(3,3)*a(1,1)-a(3,1)*a(1,3)
      d33=a(1,1)*a(2,2)-a(1,2)*a(2,1)
      d12=a(2,3)*a(3,1)-a(2,1)*a(3,3)
      d23=a(3,1)*a(1,2)-a(3,2)*a(1,1)
      d31=a(1,2)*a(2,3)-a(1,3)*a(2,2)
      d13=a(2,1)*a(3,2)-a(3,1)*a(2,2)
      d21=a(3,2)*a(1,3)-a(1,2)*a(3,3)
      d32=a(1,3)*a(2,1)-a(2,3)*a(1,1)
      c=a(1,1)*d11+a(1,2)*d12+a(1,3)*d13
      b(1,1)=d11/c
      b(2,2)=d22/c
      b(3,3)=d33/c
      b(1,2)=d21/c
      b(2,3)=d32/c
      b(3,1)=d13/c
      b(2,1)=d12/c
      b(3,2)=d23/c
      b(1,3)=d31/c

      return
      end
!
! ---------------------------------------------------------------------
!
      subroutine alloc_nodes(nodes,ierror)

      integer, intent(in) :: nodes
      integer, intent(out) :: ierror
      integer ialloc(2)

!     if(mynod.eq.0) write(6,*)'alloc_nodes(',nodes,')'

      ialloc(:) = 0

      allocate(nbuf(0:nodes), stat=ialloc(1))
      allocate(nbufa(0:nodes), stat=ialloc(2))

      ierror = 0
      do i=1,2
       ierror = ierror + ialloc(i)
      enddo

      return
      end subroutine        ! alloc_nodes !
!
! ---------------------------------------------------------------------
!
      subroutine alloc_buffers(ierror)

      integer, intent(out) :: ierror
      integer ialloc(14)

!     if(mynod.eq.0) write(6,*)'alloc_buffers: nbufsize=',nbufsize

      ialloc(:) = 0

      allocate(buf(nbufsize), stat=ialloc(1))
      allocate(bufa(nbufsize), stat=ialloc(2))
      allocate(bufb(nbufsize), stat=ialloc(3))
      allocate(bufc(nbufsize), stat=ialloc(4))

      allocate(ibufY1(natoms_alloc), stat=ialloc(5))
      allocate(ibufY2(natoms_alloc), stat=ialloc(6))
      allocate(ibufZ1(natoms_alloc), stat=ialloc(7))
      allocate(ibufZ2(natoms_alloc), stat=ialloc(8))

      allocate(kbuf(nbufsize), stat=ialloc(9))
      allocate(kbufa(nbufsize), stat=ialloc(10))
      allocate(kbufb(nbufsize), stat=ialloc(11))
      allocate(kbufc(nbufsize), stat=ialloc(12))
      allocate(ncell_of_atom(natoms_alloc), stat=ialloc(13))

      ierror = 0
      do i=1,14
       ierror = ierror + ialloc(i)
      enddo

      return
      end subroutine        ! alloc_buffers !
!
! ---------------------------------------------------------------------
!
      subroutine deall_buffers(ierror)

      integer, intent(out) :: ierror
      integer ialloc(14)

      ialloc(:) = 0

      if (allocated(buf)) deallocate(buf, stat=ialloc(1))
      if (allocated(bufa)) deallocate(bufa, stat=ialloc(2))
      if (allocated(bufb)) deallocate(bufb, stat=ialloc(3))
      if (allocated(bufc)) deallocate(bufc, stat=ialloc(4))

      if (allocated(ibufY1)) deallocate(ibufY1, stat=ialloc(5))
      if (allocated(ibufY2)) deallocate(ibufY2, stat=ialloc(6))
      if (allocated(ibufZ1)) deallocate(ibufZ1, stat=ialloc(7))
      if (allocated(ibufZ2)) deallocate(ibufZ2, stat=ialloc(8))

      if (allocated(kbuf)) deallocate(kbuf, stat=ialloc(9))
      if (allocated(kbufa)) deallocate(kbufa, stat=ialloc(10))
      if (allocated(kbufb)) deallocate(kbufb, stat=ialloc(11))
      if (allocated(kbufc)) deallocate(kbufc, stat=ialloc(12))
      if (allocated(ncell_of_atom)) deallocate(ncell_of_atom, 
     1    stat=ialloc(13))

      ierror = 0
      do i=1,14
       ierror = ierror + ialloc(i)
      enddo

      return
      end subroutine        ! deall_buffers !
!
! ---------------------------------------------------------------------
!
      subroutine alloc_cells(ierror)

      integer, intent(out) :: ierror
      integer ialloc(4)

      ialloc(:) = 0

      allocate(id_of_cell(ncell_per_node), stat=ialloc(1))
      allocate(natoms_in_cell(0:ncell_per_node), stat=ialloc(2))

      allocate(n_in_cell(natoms_per_cell,0:ncell_per_node),
     1 stat=ialloc(3))

      ierror = 0
      do i=1,4
       ierror = ierror + ialloc(i)
      enddo

      return
      end subroutine        ! alloc_cells !
!
! ---------------------------------------------------------------------
!
      subroutine deall_cells(ierror)

      integer, intent(out) :: ierror
      integer ialloc(4)

      ialloc(:) = 0

      if (allocated(id_of_cell)) deallocate(id_of_cell, stat=ialloc(1))
      if (allocated(natoms_in_cell)) deallocate(natoms_in_cell, 
     1 stat=ialloc(2))
      if (allocated(n_in_cell)) deallocate(n_in_cell, stat=ialloc(3))

      ierror = 0
      do i=1,4
       ierror = ierror + ialloc(i)
      enddo

      return
      end subroutine        ! deall_cells !
!
! -------------------------------------------------------------------
!
      subroutine error_check(ierror,message)
      implicit double precision (a-h,o-z)

      integer, intent(in) :: ierror
      character(*) message

      if(ierror.ne.0) then
       if(mynod.eq.0) write(6,10) message
       call PROGRAM_END(1)
      endif

 10   format(/,A)

      return
      end subroutine       ! error_check !


      END MODULE  ! sim_box !
!
! =====================================================================
!
      MODULE atoms

      use sim_box

      save
      integer :: nbr_tot, max_nbrs, max1_nbrs, max2_nbrs

      integer :: iatom_types, ipair_types, icomp_high, icomp_low

      double precision :: at_vol, at_vol_inp, sys_vol

      double precision, dimension(3) :: sum_dcm, Sys_dcmR,
     1 Sys_dcmS = 0.d0
      ! System center of mass deviation !

      integer, dimension(:), allocatable :: ident, ntype
      !dir$ attributes align:64 :: ident
      !dir$ attributes align:64 :: ntype

      ! ident() -> id number of an atom              !
      ! ntype() -> chem. type of an atom inside code !

      double precision, dimension(:), allocatable :: sx,sy,sz,
     1 rx,ry,rz
      !dir$ attributes align:64 :: sx
      !dir$ attributes align:64 :: sy
      !dir$ attributes align:64 :: sz
      !dir$ attributes align:64 :: rx
      !dir$ attributes align:64 :: ry
      !dir$ attributes align:64 :: rz

      integer, dimension(:,:), allocatable :: nbr_list
      !dir$ attributes align:64 :: nbr_list

!     MD atom arrays:

      double precision :: dt_step = 0.001d0 ! default MD step  !

      double precision :: A_fr=0.d0, Q_heat=0.d0
      double precision :: Ek_wall=0.d0

      double precision, dimension(:), allocatable :: x1,y1,z1,
     1 x2,y2,z2, x3,y3,z3, x4,y4,z4, x5,y5,z5,
     2 sumPx, sumPy, sumPz, Tscale
      !dir$ attributes align:64 :: x1
      !dir$ attributes align:64 :: x2
      !dir$ attributes align:64 :: x3
      !dir$ attributes align:64 :: x4
      !dir$ attributes align:64 :: x5
      !dir$ attributes align:64 :: y1
      !dir$ attributes align:64 :: y2
      !dir$ attributes align:64 :: y3
      !dir$ attributes align:64 :: y4
      !dir$ attributes align:64 :: y5
      !dir$ attributes align:64 :: z1
      !dir$ attributes align:64 :: z2
      !dir$ attributes align:64 :: z3
      !dir$ attributes align:64 :: z4
      !dir$ attributes align:64 :: z5

      double precision, dimension(:,:), allocatable :: frr
      !dir$ attributes align:64 :: frr

      double precision, dimension(:), allocatable :: Ep_of
      !dir$ attributes align:64 :: Ep_of


      CONTAINS
!
! ---------------------------------------------------------------------
!
      subroutine alloc_atoms_sys(ierror)

      integer, intent(out) :: ierror
      integer ialloc(12)

      ialloc(:) = 0

      allocate(ident(natoms_alloc), stat=ialloc(1))
      allocate(ntype(natoms_alloc), stat=ialloc(2))

      allocate(rx(natoms_alloc), stat=ialloc(3))
      allocate(ry(natoms_alloc), stat=ialloc(4))
      allocate(rz(natoms_alloc), stat=ialloc(5))

      allocate(sx(natoms_alloc), stat=ialloc(6))
      allocate(sy(natoms_alloc), stat=ialloc(7))
      allocate(sz(natoms_alloc), stat=ialloc(8))

      allocate(nbr_list(0:nbrs_per_atom,natoms_alloc), stat=ialloc(9))

      allocate(frr(3,natoms_alloc), stat=ialloc(10))
      allocate(Ep_of(natoms_alloc), stat=ialloc(11))

      ierror = 0
      do i=1,12
       ierror = ierror + ialloc(i)
      enddo

      return
      end subroutine        ! alloc_atoms_sys !
!
! ---------------------------------------------------------------------
!
      subroutine deall_atoms_sys(ierror)

      integer, intent(out) :: ierror
      integer ialloc(12)

      ialloc(:) = 0

      if (allocated(ident)) deallocate(ident, stat=ialloc(1))
      if (allocated(ntype)) deallocate(ntype, stat=ialloc(2))

      if (allocated(rx)) deallocate(rx, stat=ialloc(3))
      if (allocated(ry)) deallocate(ry, stat=ialloc(4))
      if (allocated(rz)) deallocate(rz, stat=ialloc(5))

      if (allocated(sx)) deallocate(sx, stat=ialloc(6))
      if (allocated(sy)) deallocate(sy, stat=ialloc(7))
      if (allocated(sz)) deallocate(sz, stat=ialloc(8))

      if (allocated(nbr_list)) deallocate(nbr_list, stat=ialloc(9))

      if (allocated(frr)) deallocate(frr, stat=ialloc(10))
      if (allocated(Ep_of))deallocate(Ep_of, stat=ialloc(11))

      ierror = 0
      do i=1,12
       ierror = ierror + ialloc(i)
      enddo

      return
      end subroutine        ! deall_atoms_sys !
!
! ---------------------------------------------------------------------
!
      subroutine alloc_atoms_MD(ierror)

      implicit double precision (a-h,o-z)

      integer, intent(out) :: ierror
      integer ialloc(20)

      ialloc(:) = 0

      allocate(x1(natoms_alloc), stat=ialloc(1))
      allocate(y1(natoms_alloc), stat=ialloc(2))
      allocate(z1(natoms_alloc), stat=ialloc(3))

      allocate(x2(natoms_alloc), stat=ialloc(4))
      allocate(y2(natoms_alloc), stat=ialloc(5))
      allocate(z2(natoms_alloc), stat=ialloc(6))

      allocate(x3(natoms_alloc), stat=ialloc(7))
      allocate(y3(natoms_alloc), stat=ialloc(8))
      allocate(z3(natoms_alloc), stat=ialloc(9))

      allocate(x4(natoms_alloc), stat=ialloc(10))
      allocate(y4(natoms_alloc), stat=ialloc(11))
      allocate(z4(natoms_alloc), stat=ialloc(12))

      allocate(x5(natoms_alloc), stat=ialloc(13))
      allocate(y5(natoms_alloc), stat=ialloc(14))
      allocate(z5(natoms_alloc), stat=ialloc(15))

      allocate(sumPx(iatom_types), stat=ialloc(16))
      allocate(sumPy(iatom_types), stat=ialloc(17))
      allocate(sumPz(iatom_types), stat=ialloc(18))
      allocate(Tscale(iatom_types), stat=ialloc(19))

      x1(:) = 0.d0; y1(:) = 0.d0; z1(:) = 0.d0
      sumPx(:) = 0.d0; sumPy(:) = 0.d0; sumPz(:) = 0.d0

      ierror = 0
      do i=1,20
       ierror = ierror + ialloc(i)
      enddo

      return
      end subroutine        ! alloc_atoms_MD !
!
! ---------------------------------------------------------------------
!
      subroutine deall_atoms_MD(ierror)

      implicit double precision (a-h,o-z)

      integer, intent(out) :: ierror
      integer ialloc(20)

      ialloc(:) = 0

      if (allocated(x1)) deallocate(x1, stat=ialloc(1))
      if (allocated(y1)) deallocate(y1, stat=ialloc(2))
      if (allocated(z1)) deallocate(z1, stat=ialloc(3))

      if (allocated(x2)) deallocate(x2, stat=ialloc(4))
      if (allocated(y2)) deallocate(y2, stat=ialloc(5))
      if (allocated(z2)) deallocate(z2, stat=ialloc(6))

      if (allocated(x3)) deallocate(x3, stat=ialloc(7))
      if (allocated(y3)) deallocate(y3, stat=ialloc(8))
      if (allocated(z3)) deallocate(z3, stat=ialloc(9))

      if (allocated(x4)) deallocate(x4, stat=ialloc(10))
      if (allocated(y4)) deallocate(y4, stat=ialloc(11))
      if (allocated(z4)) deallocate(z4, stat=ialloc(12))

      if (allocated(x5)) deallocate(x5, stat=ialloc(13))
      if (allocated(y5)) deallocate(y5, stat=ialloc(14))
      if (allocated(z5)) deallocate(z5, stat=ialloc(15))

      if (allocated(sumPx)) deallocate(sumPx, stat=ialloc(16))
      if (allocated(sumPy)) deallocate(sumPy, stat=ialloc(17))
      if (allocated(sumPz)) deallocate(sumPz, stat=ialloc(18))
      if (allocated(Tscale)) deallocate(Tscale, stat=ialloc(19))

      ierror = 0
      do i=1,20
       ierror = ierror + ialloc(i)
      enddo

      return
      end subroutine        ! deall_atoms_MD !

      END MODULE  ! atoms !
!
! =====================================================================
!
      MODULE pot_module

      use sim_box
      use atoms

      save
      integer, parameter :: maxconstr=63 ! Max number of defined    !
                                         ! constraints              !
      character*2 :: chem_symb(0:112)
      integer ipot(112)            ! maximum numb. elements possible !

      character*4, dimension(:), allocatable :: elem_symb, 
     1 Elem_in_com
      character*80, allocatable :: filename(:)
     
      double precision :: PotEnrg_atm,PotEnrg_glb,ecoh

      integer :: ifile_numbs,nelem_in_com,n_pot_files,
     1 INP_STR_TYPE, iOUT_STR_TYPE, N_points_max, Nrho_points_max

      integer :: maxf  ! Max number of potential files !
      integer :: iPOT_file_type ! 0 - plt, *.dat files !
                                ! 1 - plt, *.plt files !
				! 2 - LAMMPS files     !

      integer :: iPOT_file_ver = 0 ! default: old ANN: normal Gis !
                                   ! 1 - log Gis for ANN pot. !

      integer :: iPOT_func_type ! 0 - EAM !
                                ! 1 - ADP !

      integer :: inv_calls      ! for TEST only !

      double precision :: rindr, dr_min, dr_recp
      double precision :: drho_min, rho_max, drho_recp
      double precision :: r_cut_off, r2_cut_off, rin_max, avr_mass
      double precision :: r_cut_short, r2_cut_short, r_max_pot
      double precision :: elem_radius(0:112)

      integer, dimension(:), allocatable :: ielement,
     1 natoms_of_type, iZ_elem_in_com

      double precision, dimension(:), allocatable :: amass,
     1 sum_mass, gram_mol, chem_pot, Am_of_type, pTemp, E_kin

      double precision, dimension(:), allocatable :: 
     1 r_pair_cut, r2_pair_cut
      !dir$ attributes align:64 :: r_pair_cut
      !dir$ attributes align:64 :: r2_pair_cut

      CONTAINS
!
! -------------------------------------------------------------------
!   Changes a string to upper case
! -------------------------------------------------------------------
!
      subroutine conv_to_low (strIn,strOut,len_str)

      character(*) strIn, strOut
      integer len_str

      len_str = len_trim(strIn)

!     write(6,10)strIn,len_trim(strIn)
! 10  format('conv_to_low: len(',A,')=',i3)

      do i = 1, len_str
         j = iachar(strIn(i:i))
         if (j>= iachar("A") .and. j<=iachar("Z") ) then
                 strOut(i:i) = achar(iachar(strIn(i:i))+32)
         else
                 strOut(i:i) = strIn(i:i)
         endif
      enddo

      return
      end subroutine  ! conv_to_low !
!
!--------------------------------------------------------------------
!
      subroutine init_elements

      implicit double precision (a-h,o-z)

       chem_symb(0) = ''; elem_radius(0) = 0.0d0 ! Ang !
       chem_symb(1) = 'H'; elem_radius(1) = 0.53
       chem_symb(2) = 'He'; elem_radius(2) = 0.31
       chem_symb(3) = 'Li'; elem_radius(3) = 1.45
       chem_symb(4) = 'Be'; elem_radius(4) = 1.05
       chem_symb(5) = 'B'; elem_radius(5) = 0.85
       chem_symb(6) = 'C'; elem_radius(6) = 0.70
       chem_symb(7) = 'N'; elem_radius(7) = 0.65
       chem_symb(8) = 'O'; elem_radius(8) = 0.60
       chem_symb(9) = 'F'; elem_radius(9) = 0.50
       chem_symb(10) = 'Ne'; elem_radius(10) = 0.38
       chem_symb(11) = 'Na'; elem_radius(11) = 1.80
       chem_symb(12) = 'Mg'; elem_radius(12) = 1.50
       chem_symb(13) = 'Al'; elem_radius(13) = 1.18
       chem_symb(14) = 'Si'; elem_radius(14) = 1.10
       chem_symb(15) = 'P'; elem_radius(15) = 1.00
       chem_symb(16) = 'S'; elem_radius(16) = 1.00
       chem_symb(17) = 'Cl'; elem_radius(17) = 1.00
       chem_symb(18) = 'Ar'; elem_radius(18) = 0.71
       chem_symb(19) = 'K'; elem_radius(19) = 2.20
       chem_symb(20) = 'Ca'; elem_radius(20) = 1.80
       chem_symb(21) = 'Sc'; elem_radius(21) = 1.60
       chem_symb(22) = 'Ti'; elem_radius(22) = 1.40
       chem_symb(23) = 'V'; elem_radius(23) = 1.35
       chem_symb(24) = 'Cr'; elem_radius(24) = 1.40
       chem_symb(25) = 'Mn'; elem_radius(25) = 1.40
       chem_symb(26) = 'Fe'; elem_radius(26) = 1.40
       chem_symb(27) = 'Co'; elem_radius(27) = 1.35
       chem_symb(28) = 'Ni'; elem_radius(28) = 1.35
       chem_symb(29) = 'Cu'; elem_radius(29) = 1.35
       chem_symb(30) = 'Zn'; elem_radius(30) = 1.35
       chem_symb(31) = 'Ga'; elem_radius(31) = 1.30
       chem_symb(32) = 'Ge'; elem_radius(32) = 1.25
       chem_symb(33) = 'As'; elem_radius(33) = 1.15
       chem_symb(34) = 'Se'; elem_radius(34) = 1.15
       chem_symb(35) = 'Br'; elem_radius(35) = 1.15
       chem_symb(36) = 'Kr'; elem_radius(36) = 0.88
       chem_symb(37) = 'Rb'; elem_radius(37) = 2.35
       chem_symb(38) = 'Sr'; elem_radius(38) = 2.00
       chem_symb(39) = 'Y'; elem_radius(39) = 1.85
       chem_symb(40) = 'Zr'; elem_radius(40) = 1.55
       chem_symb(41) = 'Nb'; elem_radius(41) = 1.45
       chem_symb(42) = 'Mo'; elem_radius(42) = 1.45
       chem_symb(43) = 'Tc'; elem_radius(43) = 1.35
       chem_symb(44) = 'Ru'; elem_radius(44) = 1.30
       chem_symb(45) = 'Rh'; elem_radius(45) = 1.35
       chem_symb(46) = 'Pd'; elem_radius(46) = 1.40
       chem_symb(47) = 'Ag'; elem_radius(47) = 1.60
       chem_symb(48) = 'Cd'; elem_radius(48) = 1.55
       chem_symb(49) = 'In'; elem_radius(49) = 1.55
       chem_symb(50) = 'Sn'; elem_radius(50) = 1.45
       chem_symb(51) = 'Sb'; elem_radius(51) = 1.45
       chem_symb(52) = 'Te'; elem_radius(52) = 1.40
       chem_symb(53) = 'I'; elem_radius(53) = 1.40
       chem_symb(54) = 'Xe'; elem_radius(54) = 1.08
       chem_symb(55) = 'Cs'; elem_radius(55) = 2.60
       chem_symb(56) = 'Ba'; elem_radius(56) = 2.15
       chem_symb(57) = 'La'; elem_radius(57) = 1.95
       chem_symb(58) = 'Ce'; elem_radius(58) = 1.85
       chem_symb(59) = 'Pr'; elem_radius(59) = 1.85
       chem_symb(60) = 'Nd'; elem_radius(60) = 1.85
       chem_symb(61) = 'Pm'; elem_radius(61) = 1.85
       chem_symb(62) = 'Sm'; elem_radius(62) = 1.85
       chem_symb(63) = 'Eu'; elem_radius(63) = 1.85
       chem_symb(64) = 'Gd'; elem_radius(64) = 1.80
       chem_symb(65) = 'Tb'; elem_radius(65) = 1.75
       chem_symb(66) = 'Dy'; elem_radius(66) = 1.75
       chem_symb(67) = 'Ho'; elem_radius(67) = 1.75
       chem_symb(68) = 'Er'; elem_radius(68) = 1.75
       chem_symb(69) = 'Tm'; elem_radius(69) = 1.75
       chem_symb(70) = 'Yb'; elem_radius(70) = 1.75
       chem_symb(71) = 'Lu'; elem_radius(71) = 1.75
       chem_symb(72) = 'Hf'; elem_radius(72) = 1.55
       chem_symb(73) = 'Ta'; elem_radius(73) = 1.45
       chem_symb(74) = 'W'; elem_radius(74) = 1.35
       chem_symb(75) = 'Re'; elem_radius(75) = 1.35
       chem_symb(76) = 'Os'; elem_radius(76) = 1.30
       chem_symb(77) = 'Ir'; elem_radius(77) = 1.35
       chem_symb(78) = 'Pt'; elem_radius(78) = 1.35
       chem_symb(79) = 'Au'; elem_radius(79) = 1.35
       chem_symb(80) = 'Hg'; elem_radius(80) = 1.50
       chem_symb(81) = 'Tl'; elem_radius(81) = 1.90
       chem_symb(82) = 'Pb'; elem_radius(82) = 1.80
       chem_symb(83) = 'Bi'; elem_radius(83) = 1.60
       chem_symb(84) = 'Po'; elem_radius(84) = 1.90
       chem_symb(85) = 'At'; elem_radius(85) = 1.27
       chem_symb(86) = 'Rn'; elem_radius(86) = 1.20
       chem_symb(87) = 'Fr'; elem_radius(87) = 1.94
       chem_symb(88) = 'Ra'; elem_radius(88) = 2.15
       chem_symb(89) = 'Ac'; elem_radius(89) = 1.95
       chem_symb(90) = 'Th'; elem_radius(90) = 1.80
       chem_symb(91) = 'Pa'; elem_radius(91) = 1.80
       chem_symb(92) = 'U'; elem_radius(92) = 1.75
       chem_symb(93) = 'Np'; elem_radius(93) = 1.75
       chem_symb(94) = 'Pu'; elem_radius(94) = 1.75
       chem_symb(95) = 'Am'; elem_radius(95) = 1.75
       chem_symb(96) = 'Cm'; elem_radius(96) = 1.75
       chem_symb(97) = 'Bk'; elem_radius(97) = 1.75
       chem_symb(98) = 'Cf'; elem_radius(98) = 1.75
       chem_symb(99) = 'Es'; elem_radius(99) = 1.75
       chem_symb(100) = 'Fm'; elem_radius(100) = 1.75
       chem_symb(101) = 'Md'; elem_radius(101) = 1.75
       chem_symb(102) = 'No'; elem_radius(102) = 1.75
       chem_symb(103) = 'Lr'; elem_radius(103) = 1.75
       chem_symb(104) = 'Rf'; elem_radius(104) = 1.75
       chem_symb(105) = 'Db'; elem_radius(105) = 1.75
       chem_symb(106) = 'Sg'; elem_radius(106) = 1.75
       chem_symb(107) = 'Bh'; elem_radius(107) = 1.75
       chem_symb(108) = 'Hs'; elem_radius(108) = 1.75
       chem_symb(109) = 'Mt'; elem_radius(109) = 1.75
       chem_symb(110) = 'Ds'; elem_radius(110) = 1.75
       chem_symb(111) = 'Rg'; elem_radius(111) = 1.75
       chem_symb(112) = 'Cn'; elem_radius(112) = 1.75

      return     
      end subroutine       ! init_elements !
!
!--------------------------------------------------------------------
!
      function numb_elem_Z(chem)

      implicit double precision (a-h,o-z)

      character*4 chem,chem_low
      character*2 chem_symb_low

      numZ = 0
      call conv_to_low (chem,chem_low,len_chem)
      do iZ = 1,112
       call conv_to_low (chem_symb(iZ),chem_symb_low,len_symb)
       if(len_chem.eq.len_symb) then
        if(chem_low(1:len_symb).eq.chem_symb_low(1:len_symb)) then
         numZ = iZ
         exit
        endif
       endif
      enddo

!     write(6,*)'numb_elem_Z(',chem,')=',numZ

      numb_elem_Z = numZ

      return     
      end function       ! numb_elem_Z !
!
!---------------------------------------------------------------------
! *** Gets current chem. composition ***
! Includes all atoms that DO NOT have chemical constrain!
!---------------------------------------------------------------------
!
      subroutine get_chem
      
      implicit double precision (a-h,o-z)

! *** Collect atom types...
      natoms_of_type(:) = 0

      do n=1,natoms
       ntp = ntype(n)
       natoms_of_type(ntp) = natoms_of_type(ntp) + 1 
      enddo

      return
      end subroutine       ! get_chem !
!
! ---------------------------------------------------------------------
!
      subroutine alloc_pot_types(ierror)

      integer, intent(out) :: ierror
      integer ialloc(14)

      ialloc(:) = 0

      allocate(elem_symb(iatom_types), stat=ialloc(1))

      allocate(ielement(0:iatom_types), stat=ialloc(2))
      allocate(natoms_of_type(iatom_types), stat=ialloc(3))
      allocate(iZ_elem_in_com(0:iatom_types), stat=ialloc(4))

      allocate(amass(iatom_types), stat=ialloc(5))
      allocate(sum_mass(iatom_types), stat=ialloc(6))
      allocate(gram_mol(iatom_types), stat=ialloc(7))

      allocate(Am_of_type(iatom_types), stat=ialloc(8))
      allocate(Elem_in_com(iatom_types), stat=ialloc(9))

      allocate(r_pair_cut(ipair_types), stat=ialloc(10))
      allocate(r2_pair_cut(ipair_types), stat=ialloc(11))

      allocate(E_kin(iatom_types), stat=ialloc(12))
      allocate(pTemp(0:iatom_types), stat=ialloc(13))
      allocate(filename(1), stat=ialloc(14))

      ierror = 0
      do i=1,14
       ierror = ierror + ialloc(i)
      enddo

      return
      end subroutine        ! alloc_pot_types !
!
! ---------------------------------------------------------------------
!
      subroutine deall_pot_types(ierror)

      integer, intent(out) :: ierror
      integer ialloc(14)

      ialloc(:) = 0

      if (allocated(elem_symb)) deallocate(elem_symb, stat=ialloc(1))
      if (allocated(filename)) deallocate(filename, stat=ialloc(2))

      if (allocated(ielement)) deallocate(ielement, stat=ialloc(3))
      if (allocated(natoms_of_type)) deallocate(natoms_of_type,
     1    stat=ialloc(4))
      if (allocated(iZ_elem_in_com)) deallocate(iZ_elem_in_com,
     1    stat=ialloc(5))

      if (allocated(amass)) deallocate(amass, stat=ialloc(6))
      if (allocated(sum_mass)) deallocate(sum_mass, stat=ialloc(7))
      if (allocated(gram_mol)) deallocate(gram_mol, stat=ialloc(8))
      if (allocated(Am_of_type)) deallocate(Am_of_type, 
     1    stat=ialloc(9))
      if (allocated(Elem_in_com)) deallocate(Elem_in_com,
     1    stat=ialloc(10))

      if (allocated(r_pair_cut)) deallocate(r_pair_cut,stat=ialloc(11))
      if (allocated(r2_pair_cut))deallocate(r2_pair_cut,stat=ialloc(12))

      if (allocated(E_kin)) deallocate(E_kin, stat=ialloc(13))
      if (allocated(pTemp)) deallocate(pTemp, stat=ialloc(14))

      ierror = 0
      do i=1,14
       ierror = ierror + ialloc(i)
      enddo

      return
      end subroutine        ! deall_pot_types !
!
!---------------------------------------------------------------------
!
      END MODULE  ! pot_module !
!
!=====================================================================
!
