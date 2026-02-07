
!***********************************************************************
! MAIN PROGRAM
!***********************************************************************

program lycom
use lycom_par
use lycom_common
use lycom_interface
use lycom_local
use lycom_nova
use lycom_land
use lycom_global
use lycom_global2
use lycom_opt
implicit none

integer :: month0, day0

! initialisation

call init_mpi ! initialise MPI
!if (para) then
!  if (rank .eq. 0) call lycom_openG ! open forcing files (global)
!else
!  call lycom_openL ! local
!endif

if (para) then

  if (rank .eq. 0) call read_land ! read land mask file
  !write(*,*) "indvec value 1 readland", indvec(:,1)
  !write(*,*) "indvec value 2 readland", indvec(:,2)
  call bc_pp ! broadcast points per processor and other parameters

  nCPts = pp
    
  call alloc_global ! allocate global fields     # if (rank .eq. 0) 
endif

if (rank .eq. 0) call lycom_namelist ! read namelist

if (para) then
  call bc_namelist
else
  nCPts = nSites
endif

call lycom_alloc ! allocate variables

if (BSCtypes) call lycom_allocBSC

if (rank .eq. 0) then
  call lycom_readSpec ! read strategy parameters
endif
if (para) then
  call bc_specpar ! broadcast strategies

  if (rank .eq. 0) call lycom_readBCg ! read boundary conditions (global)

  call sc_BC ! scatter boundary conditions
else
  call inq_ntLoc

  call lycom_readBCl ! local
endif

!if (BSCtypes .and. NOHONO) then
!  if (rank .eq. 0) call lycom_readNOHONO ! read NO/HONO fluxes (optional)
  
!  if (para) call bc_NOHONO
!endif

if (para) then

  call lycom_init(ppvec(rank+1)) ! initialise variables (global)

else
  call lycom_init(nCPts) ! local

!  call lycom_init_speciesL
endif

call lycom_reset ! set accumulation variables to zero

if (para) then
!This is commented
!  if (rank .eq. 0) call lycom_openPreprocessing
  if (rank .eq. 0) call lycom_openG ! open forcing files (global)
  !write(*,*) "OPenG done"
else
  call lycom_openL ! local
endif

if (.not. para) open(kfile_survspec, file=sfile_survspec, status='replace', action='write') ! open survivors file

year = year0
cyear = cyear0
accts = accts0
tspd = 86400 / tsl

months30(:) = (/ 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12 /)

tsindata = tsindata0
tpos = 1
out1 = .true.

if (para .and. rank .eq. 0) open( kstatus, file=sstatus, status='replace', action='write' )

! main yearly loop

do while ( year .le. lastyear )

  if ( mod( cyear, 4 ) .eq. 0 ) then ! determine leap year
    leapyear = .true.
    if ( mod( cyear, 100 ) .eq. 0 .and. mod( cyear, 400 ) .ne. 0 ) leapyear = .false.
  else
    leapyear = .false.
  endif

  if ( year .ge. yearout1 .and. year .le. yearoutX ) then ! check for writing output
    writeout = .true.
  else
    writeout = .false.
  endif

! monthly loop

  do month0 = 1, 12

    month = month0

    if ( any( months30 .eq. month )) then ! determine days per month
      dpm = 30
    elseif ( month .eq. 2 ) then
      dpm = 28
      if ( leapyear ) dpm = 29
    else
      dpm = 31
    endif

    call lycom_readMon(month) ! read monthly forcing (from variable)

! daily loop

    do day0 = 1, dpm

      day = day0

! timestep loop

      do ts = 1, tspd

        if (para) then
          !write(*,*) "lycom_readHourG is called"
          !write(*,*) "rank is", rank

          !write(*,*) "indvec value 1 before Hourreadland", indvec(:,1)
          !write(*,*) "indvec value 2 before Houreadland", indvec(:,2)
          if (rank .eq. 0) call lycom_readHourG (tsindata) ! read hourly climate forcing (global)
          
          call sc_clim ! scatter climate forcing

          !fH2Ol_ad(:) = fH2Ol_ad(:) !* 0.001 ! adapt units for ERA5 or WATCH data (mm/s -> m/s)
          !fH2Os_ad(:) = fH2Os_ad(:) !* 0.001 
        else
          call lycom_readHourL ! local

          fH2Ol_ad(:) = fH2Ol_ad(:) !* 0.001
          fH2Os_ad(:) = fH2Os_ad(:) !* 0.001 
        endif

        if (para) then

          call lycom_step(ppvec(rank+1)) ! perform calculations
        else
          call lycom_step(nCPts)
        endif

        if ( writeout .and. mod( accts, outint ) .eq. 0 ) then ! output

          if (para) then

            call lycom_av_output(ppvec(rank+1)) ! average accumulated variables

            call lycom_outputG ! write global output

            tpos = tpos + 1
          else

            call lycom_av_output(nCPts) ! local

            if (out1) call lycom_openOutL

            call lycom_outputL
          endif

          if (out1) out1 = .false.

          call lycom_reset ! reset accumulated variables
        endif

!        if ( year .eq. lastyear .and. outdiurnal .and. .not. para) call lycom_outdnl ! write diurnal output

        accts = accts + 1

        tsindata = tsindata + 1

        if (tsindata .gt. tstep_total) then  ! shift to end of year ?
          
          tsindata = 1

          cyear = cyear0 -1
        endif

      enddo ! time step loop

      if (.not. para) call lycom_outsurv ! write number of surviving species

    enddo ! daily loop

  enddo ! monthly loop

  if (para) then
    if (rank .eq. 0) write( kstatus,* ) "year ",year," finished"
  else
    write(*,*) "year ",year," finished"
  endif

  year = year + 1 ! update simulation year

  cyear = cyear + 1 ! update calendar year

enddo ! yearly loop

! output at end of simulation

if (para) then

  call lycom_av_species(ppvec(rank+1)) ! average species properties

  tpos = 1 ! rewind

  call lycom_outspecG ! write averaged properties

  if (rank .eq. 0) call lycom_closeG ! close forcing
else

  call lycom_av_species(nCPts) ! local

!  call lycom_av_speciesL

  call lycom_outspecL ! write avg properties and weight of each species per habitat

  call lycom_closeL ! close forcing and output files

  close( kfile_survspec )
endif

call lycom_dealloc

if (BSCtypes) call lycom_deallocBSC

if (para) then

  call dealloc_global ! and stop MPI

  if (rank .eq. 0) then
    write( kstatus,* ) "simulation finished"
    close( kstatus )
  endif
endif

end

