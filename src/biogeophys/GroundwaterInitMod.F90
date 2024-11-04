module GroundwaterInitMod
   !-----------------------------------------------------------------------
   ! Module to find neighbors for GroundwaterMod
   ! Author: Aman Shrestha
   ! Created on: 2024/11/04
   !-----------------------------------------------------------------------

   ! USES
   use shr_kind_mod          , only : r8 => shr_kind_r8
   use decompMod             , only : get_proc_global, ldecomp, procinfo
   use domainMod             , only : ldomain
   use spmdMod               , only : MPI_REAL8, MPI_INTEGER, mpicom, npes, masterproc, iam
   use perf_mod              , only : t_startf, t_stopf

   ! CIME Globals
   use shr_infnan_mod            , only : nan => shr_infnan_nan, assignment(=)

   implicit none

   ! !PUBLIC MEMBER FUNCTIONS:
   public :: neighborInit              ! initializes neighboring grid cell indices
   public :: DetermineNeighbors_GW     ! determines neighbors

   ! Used to initialize and test unset integers
   integer, parameter, public :: gw_unset_int = -9999

   logical :: debug = .false.  ! for debugging this module

!-----------------------------------------------------------------------
contains
!-----------------------------------------------------------------------
   subroutine NeighborInit()
      !-----------------------------------------------------------------------
      ! Initialize neighbor grids
      !-----------------------------------------------------------------------
      ! Local variables:
      integer :: ier          ! error status

      integer :: numg               ! number of land gridcells
      call get_proc_global(ng=numg)

      ! Allocate ldecomp vars
      allocate(ldecomp%glat(numg), stat=ier)
      allocate(ldecomp%glon(numg), stat=ier)
   
      allocate(ldecomp%gtop(numg), stat=ier)
      allocate(ldecomp%gbot(numg), stat=ier)
      allocate(ldecomp%glft(numg), stat=ier)
      allocate(ldecomp%grgt(numg), stat=ier)

      allocate(ldecomp%gtoplft(numg), stat=ier)
      allocate(ldecomp%gtoprgt(numg), stat=ier)
      allocate(ldecomp%gbotlft(numg), stat=ier)
      allocate(ldecomp%gbotrgt(numg), stat=ier)
   
      allocate(ldecomp%gneighbors(numg), stat=ier)

      ! Initialize
      ldecomp%glat(:) = nan
      ldecomp%glon(:) = nan
  
      ldecomp%gneighbors(:) = 0
     
      ldecomp%gtop(:) = 0
      ldecomp%gbot(:) = 0
      ldecomp%glft(:) = 0
      ldecomp%grgt(:) = 0
  
      ldecomp%gtoplft(:) = 0
      ldecomp%gtoprgt(:) = 0
      ldecomp%gbotlft(:) = 0
      ldecomp%gbotrgt(:) = 0

      call DetermineNeighbors_GW()

   end subroutine NeighborInit

   !--------------------------------------------------------------------------

   subroutine DetermineNeighbors_GW()
      !-----------------------------------------------------------------------
      ! Determines neighbors for groundwater
      ! Most of the codes have been shamelessly copied from Fates module.
      ! Fates tag sci.1.77.2_api.36.0.0
      !
      ! Neighbor search uses:
      ! current lat +- delta lat and current lon +- delta lon
      ! to find the neighbors. delta lat and delta lon are namelist inputs.
      !-----------------------------------------------------------------------
      
      ! Local variables:
      real(r8) :: delta_lat         ! lat spacing from namelist
      real(r8) :: delta_lon         ! lon spacing from namelist
      integer :: i, g_out, g_in, ni ! indices
      integer :: numg               ! number of land gridcells
      integer :: ier, mpierr        ! error status

      integer,  allocatable :: ncells_array(:), begg_array(:) ! number of cells and starting global grid cell index per process 
      
      associate(&
         gclat             =>    ldecomp%glat         , & ! Output: [real(r8) (:) ] gridcell latitude global
         gclon             =>    ldecomp%glon           & ! Output: [real(r8) (:) ] gridcell longitude global
         )

      ! For test only at 0.9x1.25 resolution
      delta_lat = 0.9_r8 
      delta_lon = 1.25_r8

      ! Get total number of grid cells
      call get_proc_global(ng=numg)

      ! ! Allocate and initialize local lat and lon arrays
      ! allocate(gclat(numg), stat=ier)
      ! if(debug) write(*,*)'DGCN: gclat alloc: ', ier
   
      ! allocate(gclon(numg), stat=ier)
      ! if(debug) write(*,*)'DGCN: gclon alloc: ', ier
   
      ! gclon(:) = nan
      ! gclat(:) = nan

      ! Allocate and initialize MPI count and displacement values
      allocate(ncells_array(0:npes-1), stat=ier)
      if(debug) write(*,*)'DGCN: ncells alloc: ', ier

      allocate(begg_array(0:npes-1), stat=ier)
      if(debug) write(*,*)'DGCN: begg alloc: ', ier

      ncells_array(:) = gw_unset_int
      begg_array(:) = gw_unset_int

      call t_startf('gw-neighbor-allgather')

      if(debug) write(*,*)'DGCN: procinfo%begg: ', procinfo%begg
      if(debug) write(*,*)'DGCN: procinfo%ncells: ', procinfo%ncells
   
      ! Gather the sizes of the ldomain that each mpi rank is passing
      call MPI_Allgather(procinfo%ncells,1,MPI_INTEGER,ncells_array,1,MPI_INTEGER,mpicom,mpierr)
      if(debug) write(*,*)'DGCN: ncells mpierr: ', mpierr
   
      ! Gather the starting gridcell index for each ldomain 
      call MPI_Allgather(procinfo%begg,1,MPI_INTEGER,begg_array,1,MPI_INTEGER,mpicom,mpierr)
      if(debug) write(*,*)'DGCN: begg mpierr: ', mpierr
   
      ! reduce the begg_array displacements by one as MPI collectives expect zero indexed arrays
      begg_array = begg_array - 1
   
      if(debug .and. masterproc) write(*,*)'DGCN: ncells_array: ' , ncells_array
      if(debug .and. masterproc) write(*,*)'DGCN: begg_array: '   , begg_array
   
      ! Gather the domain information together into the neighbor type
      ! Note that MPI_Allgatherv is only gathering a subset of ldomain
      if(debug .and. masterproc) write(*,*)'DGCN: gathering latc'
      call MPI_Allgatherv(ldomain%latc,procinfo%ncells,MPI_REAL8,gclat,ncells_array,begg_array,MPI_REAL8,mpicom,mpierr)
   
      if(debug .and. masterproc) write(*,*)'DGCN: gathering lonc'
      call MPI_Allgatherv(ldomain%lonc,procinfo%ncells,MPI_REAL8,gclon,ncells_array,begg_array,MPI_REAL8,mpicom,mpierr)
   
      if (debug .and. masterproc) then
         write(*,*)'DGCN: sum(gclat):, sum(gclon): ', sum(gclat), sum(gclon)
      end if
   
      call t_stopf('gw-neighbor-allgather')

      ! Finding neighbors by comparing grid lat lon
      call t_startf('gw-neighbor-decomp')

      ! !------------------------------------------------------------------------------------------
      ! ! This loop has been modified from original FATES code. It now iterates through all grids.
      ! ! In Fates, the loop was gi = 1,numg-1 and gj=gi, numg

      ! ! Iterate through the grid cell indices and determine if any neighboring cells are in range
      ! gc_loop: do g_out = 1,numg ! outer loop

      !    if(debug) write(*, *) 'g_out: gclon, gclat: ', gclon(g_out), gclat(g_out)

      !    ! Seach all indices for neighbors to current grid cell index
      !    neighbor_search: do g_in = 1,numg ! inner loop

      !       if(debug) write(*, *) 'DGCN: g_out,g_in: ', g_out, g_in
            
      !       if(debug) write(*, *) 'g_in: gclon, gclat: ', gclon(g_in), gclat(g_in)

      !       ! identify neighbors with the ixy, jxy indices of grid cells
      !       latlon_check: if (gclon(g_out) == gclon(g_in) .and.  &
      !          gclat(g_out) == gclat(g_in) - delta_lat) then
               
      !          ldecomp%gneighbors(g_out) = ldecomp%gneighbors(g_out) + 1
      !          ldecomp%gtop(g_out)      = g_in

      !       else if (gclon(g_out) == gclon(g_in) + delta_lon .and.  &
      !                gclat(g_out) == gclat(g_in) - delta_lat) then
            
      !          ldecomp%gneighbors(g_out) = ldecomp%gneighbors(g_out) + 1
      !          ldecomp%gtoplft(g_out)      = g_in

      !       else if (gclon(g_out) == gclon(g_in) - delta_lon .and.  &
      !                gclat(g_out) == gclat(g_in) - delta_lat) then
            
      !          ldecomp%gneighbors(g_out) = ldecomp%gneighbors(g_out) + 1
      !          ldecomp%gtoprgt(g_out)      = g_in
         
      !       else if (gclon(g_out) == gclon(g_in)     .and.  &
      !                gclat(g_out) == gclat(g_in) + delta_lat) then
               
      !          ldecomp%gneighbors(g_out) = ldecomp%gneighbors(g_out) + 1
      !          ldecomp%gbot(g_out)      = g_in

      !       else if (gclon(g_out) == gclon(g_in) + delta_lon .and.  &
      !                gclat(g_out) == gclat(g_in) + delta_lat) then
               
      !          ldecomp%gneighbors(g_out) = ldecomp%gneighbors(g_out) + 1
      !          ldecomp%gbotlft(g_out)      = g_in

      !       else if (gclon(g_out) == gclon(g_in) - delta_lon .and.  &
      !                gclat(g_out) == gclat(g_in) + delta_lat) then
               
      !          ldecomp%gneighbors(g_out) = ldecomp%gneighbors(g_out) + 1
      !          ldecomp%gbotrgt(g_out)      = g_in
            
      !       else if (gclon(g_out) == gclon(g_in) + delta_lon .and.  &
      !                gclat(g_out) == gclat(g_in)) then

      !          ldecomp%gneighbors(g_out) = ldecomp%gneighbors(g_out) + 1
      !          ldecomp%glft(g_out)      = g_in

      !       else if (gclon(g_out) == gclon(g_in) - delta_lon .and.  &
      !                gclat(g_out) == gclat(g_in)) then
               
      !          ldecomp%gneighbors(g_out) = ldecomp%gneighbors(g_out) + 1
      !          ldecomp%grgt(g_out)      = g_in

      !       end if latlon_check
      !    end do neighbor_search
      ! end do gc_loop
      ! !-----------------------------------------------------------------------------------------

      ! This loop has been modified from Farshid's code.
      ! In Farshid's code, the loop was gi = 1,numg and gj=1, numg
      ! Fates code assigns neighbors both for inner and outer loop
      ! This should be faster than Farshid's code.

      ! Iterate through the grid cell indices and determine if any neighboring cells are in range
      gc_loop: do g_out = 1,numg-1 ! outer loop

         if(debug) write(*, *) 'g_out: gclon, gclat: ', gclon(g_out), gclat(g_out)
         
         ! Seach all indices for neighbors to current grid cell index
         neighbor_search: do g_in = g_out+1,numg ! inner loop

            if(debug) write(*, *) 'DGCN: g_out,g_in: ', g_out, g_in
            
            if(debug) write(*, *) 'g_in: gclon, gclat: ', gclon(g_in), gclat(g_in)

            ! identify neighbors with the ixy, jxy indices of grid cells
            latlon_check: if (gclon(g_out) == gclon(g_in) .and.  &
               gclat(g_out) == gclat(g_in) - delta_lat) then
               
               ldecomp%gneighbors(g_out) = ldecomp%gneighbors(g_out) + 1
               ldecomp%gtop(g_out)      = g_in

               ldecomp%gneighbors(g_in) = ldecomp%gneighbors(g_in) + 1
               ldecomp%gtop(g_in)      = g_out

            else if (gclon(g_out) == gclon(g_in) + delta_lon .and.  &
                     gclat(g_out) == gclat(g_in) - delta_lat) then
            
               ldecomp%gneighbors(g_out) = ldecomp%gneighbors(g_out) + 1
               ldecomp%gtoplft(g_out)      = g_in

               ldecomp%gneighbors(g_in) = ldecomp%gneighbors(g_in) + 1
               ldecomp%gbotrgt(g_in)      = g_out

            else if (gclon(g_out) == gclon(g_in) - delta_lon .and.  &
                     gclat(g_out) == gclat(g_in) - delta_lat) then
            
               ldecomp%gneighbors(g_out) = ldecomp%gneighbors(g_out) + 1
               ldecomp%gtoprgt(g_out)      = g_in

               ldecomp%gneighbors(g_in) = ldecomp%gneighbors(g_in) + 1
               ldecomp%gbotlft(g_in)      = g_out
         
            else if (gclon(g_out) == gclon(g_in)     .and.  &
                     gclat(g_out) == gclat(g_in) + delta_lat) then
               
               ldecomp%gneighbors(g_out) = ldecomp%gneighbors(g_out) + 1
               ldecomp%gbot(g_out)      = g_in

               ldecomp%gneighbors(g_in) = ldecomp%gneighbors(g_in) + 1
               ldecomp%gtop(g_in)      = g_out

            else if (gclon(g_out) == gclon(g_in) + delta_lon .and.  &
                     gclat(g_out) == gclat(g_in) + delta_lat) then
               
               ldecomp%gneighbors(g_out) = ldecomp%gneighbors(g_out) + 1
               ldecomp%gbotlft(g_out)      = g_in

               ldecomp%gneighbors(g_in) = ldecomp%gneighbors(g_in) + 1
               ldecomp%gtoprgt(g_in)      = g_out

            else if (gclon(g_out) == gclon(g_in) - delta_lon .and.  &
                     gclat(g_out) == gclat(g_in) + delta_lat) then
               
               ldecomp%gneighbors(g_out) = ldecomp%gneighbors(g_out) + 1
               ldecomp%gbotrgt(g_out)      = g_in

               ldecomp%gneighbors(g_in) = ldecomp%gneighbors(g_in) + 1
               ldecomp%gtoplft(g_in)      = g_out
            
            else if (gclon(g_out) == gclon(g_in) + delta_lon .and.  &
                     gclat(g_out) == gclat(g_in)) then

               ldecomp%gneighbors(g_out) = ldecomp%gneighbors(g_out) + 1
               ldecomp%glft(g_out)      = g_in

               ldecomp%gneighbors(g_in) = ldecomp%gneighbors(g_in) + 1
               ldecomp%grgt(g_in)      = g_out

            else if (gclon(g_out) == gclon(g_in) - delta_lon .and.  &
                     gclat(g_out) == gclat(g_in)) then
               
               ldecomp%gneighbors(g_out) = ldecomp%gneighbors(g_out) + 1
               ldecomp%grgt(g_out)      = g_in

               ldecomp%gneighbors(g_in) = ldecomp%gneighbors(g_in) + 1
               ldecomp%glft(g_in)      = g_out

            end if latlon_check
         end do neighbor_search
      end do gc_loop


      call t_stopf('gw-neighbor-decomp')

   end subroutine DetermineNeighbors_GW

end module GroundwaterInitMod