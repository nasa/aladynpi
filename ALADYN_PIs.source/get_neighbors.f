#define NBR_NILAST 0
!
!---------------------------------------------------------------------
! Performs loop over atoms in ALL of the neighboring cells
! Uses OMP, but prepares the system for ACC (GPU) implementation.
!---------------------------------------------------------------------
!
      subroutine get_neighbors

      use sim_box
      use atoms
      use pot_module

      implicit double precision (a-h,o-z)

      integer :: ierr
      integer, dimension(natoms_per_cell3) :: ll_nbr
!
! do loop over all cells
!     
      h11 = h(1,1); h12 = h(1,2); h13 = h(1,3)
      h22 = h(2,2); h23 = h(2,3); h33 = h(3,3)

!$OMP  PARALLEL DO DEFAULT(NONE) PRIVATE (j)
!$OMP& SHARED(natoms, nbrs_per_atom, nbr_list)
!$OMP& SCHEDULE(STATIC)
#if (NBR_NILAST==TRUE)
      do i=1,natoms
        do j=0,nbrs_per_atom
          nbr_list(j,i) = i   ! Initial state: all nbrs are self-nbrs !
        enddo
      enddo
#else
      do j=0,nbrs_per_atom
        do i=1,natoms
          nbr_list(i,j) = i   ! Initial state: all nbrs are self-nbrs !
        enddo
      enddo
#endif

      max_nbrs = 0
      sz0_cut = r_cut_off/h33

!$OMP PARALLEL DO DEFAULT(NONE) REDUCTION(max:max_nbrs)
!$OMP& PRIVATE(ic,icell,jcell,iz_nr,iyx,iy_nr,ixx,nr_in_cell,n,ll_nbr)
!$OMP& PRIVATE(l_in_cell,izl,iyl,kzn,jyn,jyl,ns,nr,sxn,syn,szn,l)
!$OMP& PRIVATE(sx0,sy0,sz0,rx0,ry0,rz0,r2,k_all,k1_all,k2_all)
!$OMP& SHARED(natoms, nnx,nny,nnz, nXYlayer, ncells_all, nbr_list)
!$OMP& SHARED(h11,h22,h33,h12,h13,h23,r2_cut_off,sz0_cut, nbrs_all_of)
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
#if (NBR_NILAST==TRUE)
           nbr_list(k_all,nr) = l
#else
           nbr_list(nr,k_all) = l
#endif
          endif ! if (r2.lt.r2_cut_off)... 
         endif ! if (abs(sz0).lt.r_cut_off)... !
        enddo ! do do k = 1, l_in_cell !

        max_nbrs = max(k_all,max_nbrs)

       enddo ! do n=1,nr_in_cell

      enddo do_cells ! do ic = 1, ncells_all !

!$OMP END PARALLEL DO

      ! ensure max_nbrs is a multiple of 8 to avoid remainder loops after vectorization
      !max_nbrs = 56
      if (mod(max_nbrs,8) .ne. 0) max_nbrs = ((max_nbrs/8)+1)*8

      return
      end          ! get_neighbors !
