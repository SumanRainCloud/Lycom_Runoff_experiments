
!###############################################################
! lycom_NC
!###############################################################

module lycom_nc
contains

! CHECK STATUS
!---------------------------------------------------------------

subroutine check(status)
use netcdf

  integer, intent(in) :: status
                     
  if(status /= nf90_noerr) then
    print *, trim(nf90_strerror(status))
    stop 2
  end if

end subroutine check

end module lycom_nc

!###############################################################
! lycom_GLOBAL2
!###############################################################

module lycom_global2
contains

! DEFINE GLOBAL OUTPUT VARIABLES
!---------------------------------------------------------------

subroutine def_varG (code0, ovID0)
use lycom_par
use netcdf
use lycom_nc
implicit none

integer, intent(in)                     :: code0
integer, intent(out)                    :: ovID0
character (len=8)                       :: vd0

! Define output variables
write(vd0,'(A4,I4)') "code",code0
call check(nf90_def_var(outID, vd0, NF90_REAL, dimIDs3, ovID0)) ! returns outvarID

return
end subroutine def_varG


! WRITE GLOBAL OUTPUT VARIABLES
!---------------------------------------------------------------
subroutine write_varG (var0, ovID1)
use lycom_par
use netcdf
use lycom_nc
use mpi
implicit none

integer                 :: i,k,l
integer                 :: mperr
integer                 :: ovID1
integer, dimension(1)   :: tcode
real                    :: test_val

real, dimension(pp) :: var0

real, allocatable, dimension(:) :: var1      !(ppnp)
real, allocatable, dimension(:) :: outvarL   !(nland)
real, allocatable, dimension(:,:) :: outvar    !(nx,ny)

allocate(var1(ppnp))
allocate(outvarL(nland))
allocate(outvar(nx,ny))

! Clean input data
do i = 1, pp
  test_val = var0(i)
  if (test_val .ne. test_val .or. abs(test_val) .gt. 1.0e30) then
    var0(i) = -9999.0
  endif
enddo

call MPI_BARRIER(MPI_COMM_WORLD, mperr)

! Gather output variable - FIXED: Use MPI_DOUBLE_PRECISION
call MPI_GATHER( var0, pp, MPI_DOUBLE_PRECISION, var1, pp, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mperr )

call MPI_BARRIER(MPI_COMM_WORLD, mperr)

if (rank .eq. 0) then

  ! assemble list of output points from all processors
  k=1
  do i = 1,numproc
    outvarL(k:k+ppvec(i)-1) = var1((i-1)*pp+1:(i-1)*pp+ppvec(i))
    k = k + ppvec(i)
  enddo

  ! Clean outvarL
  do l = 1, nland
    test_val = outvarL(l)
    if (test_val .ne. test_val .or. abs(test_val) .gt. 1.0e30) then
      outvarL(l) = -9999.0
    endif
  enddo

  ! Initialize outvar with fill value
  outvar(:,:) = -9999.0

  ! distribute list of output land points to output map
  do l=1,nland
    if (indvec1(l,1) .ge. 1 .and. indvec1(l,1) .le. nx .and. &
        indvec1(l,2) .ge. 1 .and. indvec1(l,2) .le. ny) then
      outvar(indvec1(l,1),indvec1(l,2)) = outvarL(l)
    endif
  enddo

  ! Final cleanup
  do i = 1, nx
    do k = 1, ny
      test_val = outvar(i,k)
      if (test_val .ne. test_val .or. abs(test_val) .gt. 1.0e30) then
        outvar(i,k) = -9999.0
      endif
    enddo
  enddo

  ! Write output variable
  tcode(1) = (year*100 + month)*100 + day

  call check(nf90_put_var(outID, t_varID, tcode, start = (/tpos/), count = (/1/) ))

  call check(nf90_put_var(outID, ovID1, outvar(:,:), start = (/ 1, 1, tpos /), count = (/ nx, ny, 1 /) ))
  
endif

deallocate(var1)
deallocate(outvarL)
deallocate(outvar)

return
end subroutine write_varG


end module lycom_global2


!###############################################################
! lycom_GLOBAL
!###############################################################

module lycom_global
contains

! INITIALISE MPI
!---------------------------------------------------------------

subroutine init_mpi ()
use lycom_par
use mpi
implicit none

integer :: mperr

call MPI_INIT(mperr)

call MPI_COMM_RANK(MPI_COMM_WORLD, rank, mperr)

call MPI_COMM_SIZE(MPI_COMM_WORLD, numproc, mperr)

para = .true.

allocate(ppvec(numproc))

call MPI_BARRIER(MPI_COMM_WORLD, mperr)

return
end subroutine init_mpi


! READ LAND FILE
!---------------------------------------------------------------

subroutine read_land ()
use lycom_par
use netcdf
use lycom_nc
implicit none

integer                 :: i, j, l
integer                 :: ncID1, ncID2, varID
integer                 :: ndims, xtype
integer                 :: pp0, addp, cntp
integer, dimension(2)   :: dimIDs2n
character (len=20)      :: dimName, varName0

! Open land mask file and determine resolution
call check(nf90_open(landmask, nf90_nowrite, ncID1) )
!!!DEBUG -- TEMPORARY FIX FOR PALEO !!!
!!call check(nf90_inquire_dimension(ncID1, 1, dimName, nx) )
!!call check(nf90_inquire_dimension(ncID1, 2, dimName, ny) )

call check(nf90_inquire_dimension(ncID1, 2, dimName, nx) )      !! TEMPORARY
call check(nf90_inquire_dimension(ncID1, 3, dimName, ny) )      !! TEMPORARY

!!!DEBUG!!!
write(*,*) "read_land"
write(*,*) "nx=",nx
write(*,*) "ny=",ny
write(*,*) "  "


allocate(lsdata(nx,ny))
allocate(xpos(nx))
allocate(ypos(ny))

! Get the values of the coordinates and put them in xpos & ypos
!!call check(nf90_inquire_variable(ncID1, 1, varName0, xtype, ndims, dimIDs2n))
call check(nf90_inquire_variable(ncID1, 1, varName0, xtype, ndims, dimIDs2n))
write(*,*) "Reading variable #1"
write(*,*) "  --> Name       :", trim(varName0)
write(*,*) "  --> Type       :", xtype
write(*,*) "  --> Num Dims   :", ndims
write(*,*) "  --> DimIDs     :", dimIDs2n(1:ndims)


!! TEMPORARY
call check(nf90_inq_varID(ncID1, varName0, varID))
call check(nf90_get_var(ncID1, varID, xpos))

write(*,*) "varID (xpos):", varID
write(*,*) "xpos sample:", xpos(1:min(5,nx))
write(*,*) "reached here 1"



!!call check(nf90_inquire_variable(ncID1, 2, varName0, xtype, ndims, dimIDs2n))
call check(nf90_inquire_variable(ncID1, 2, varName0, xtype, ndims, dimIDs2n))
write(*,*) "Reading variable #2"
write(*,*) "  --> Name       :", trim(varName0)
write(*,*) "  --> Type       :", xtype
write(*,*) "  --> Num Dims   :", ndims
write(*,*) "  --> DimIDs     :", dimIDs2n(1:ndims)



!! TEMPORARY
call check(nf90_inq_varID(ncID1, varName0, varID))
call check(nf90_get_var(ncID1, varID, ypos))
write(*,*) "varID (ypos):", varID
write(*,*) "ypos sample:", ypos(1:min(5,ny))
write(*,*) "reached here 2"


! Get the values of the land data and put them in lsdata
call check(nf90_inq_varID(ncID1, varName, varID))
call check(nf90_get_var(ncID1, varID, lsdata))
call check(nf90_close(ncID1))

! Open one climate data file (tair) and read length
call check(nf90_open(tairfileG, nf90_nowrite, ncID2) )
!!call check(nf90_inquire_dimension(ncID2, 3, dimName, nt) )
call check(nf90_inquire_dimension(ncID2, 1, dimName, nt) )   !! TEMPORARY
call check(nf90_close(ncID2))


!!!DEBUG!!!
write(*,*) "dimName=",dimName
write(*,*) "nt=",nt
write(*,*) "  "



nland = sum(lsdata(:,:))

allocate(indvec(nland,2))

allocate(indvec1(nland,2))
! Create list of land points
l = 0
do i = 1,nx
  do j = 1,ny

    if (lsdata(i,j) .gt. 0.0) then
      l = l + 1
      indvec(l,1) = i
      indvec(l,2) = j
      indvec1(l,1) = i
      indvec1(l,2) = j
        
      !write(*,*) "l", l
      !write(*,*) "indvec 1", indvec(l,1)
      !write(*,*) "indvec 2", indvec(l,2)
    endif
  enddo
enddo
!write (*,*) "list 1", indvec(:,1)

!write (*,*) "list 2", indvec(:,2)
! Compute land points per processor (pp)

pp0 = int(floor(real(nland / numproc)))

ppvec(:) = pp0

addp = modulo(nland, numproc)

if (addp .gt. 0) then
  pp = pp0 + 1

  cntp = 0
  do i = 1, numproc
    cntp = cntp + 1
    if (cntp .le. addp) ppvec(i) = ppvec(i) + 1
  enddo
else
  pp = pp0
endif

! total land points + extra points
ppnp = pp * numproc
write(*,*) "read land completed"
deallocate(lsdata)

return
end subroutine read_land


! BROADCAST POINTS PER PROCESSOR
!---------------------------------------------------------------

subroutine bc_pp ()
use lycom_par
use mpi
implicit none

integer :: mperr

call MPI_BCAST( pp, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mperr ) ! pp: dim=1, bcast process is root (0)

call MPI_BCAST( ppnp, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mperr )

call MPI_BCAST( nt, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mperr )

call MPI_BCAST( nx, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mperr )

call MPI_BCAST( ny, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mperr )

call MPI_BCAST( nland, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mperr )

call MPI_BCAST( ppvec, numproc, MPI_INTEGER, 0, MPI_COMM_WORLD, mperr )

return
end subroutine bc_pp


! ALLOCATE GLOBAL FIELDS
!---------------------------------------------------------------

subroutine alloc_global ()
use lycom_par
implicit none

if (rank .eq. 0) then

  allocate( indata(nx,ny), &
  indataL(nland), &
!  R_net(nland,tstep_total), &
!  ta_C(nland,tstep_total), &
!  desat_dT3(nland,tstep_total), &
!  ET_pot0(nland,tstep_total), &
!  E_Tdavg(nland), &
!  E_Trmax(nland), &
!  E_Tratio(nland,tstep_total), &
!  tair_for_avg(nland,tstep_total), &
!  srad_for_avg(nland,tstep_total), &
!  lrad_for_avg(nland,tstep_total), &
  bcdata(nx,ny,12), &
  bcdataL(nland,12), &
  bcdata2(nx,ny), &
  bcdata2L(nland) )
endif

allocate( laidata(ppnp,12), &
saidata(ppnp,12), &
ssadata(ppnp), &
biomedata(ppnp), &
ETrmaxdata(ppnp), &
ETdavgdata(ppnp), &
indata7(ppnp,7))
!indata_etmapdata(ppnp), &
!indata_etravg(ppnp))


! Initialize with fill value
laidata(:,:) = 1.0E12
saidata(:,:) = 1.0E12
ssadata(:)   = 1.0E12
biomedata(:) = 1.0E12
ETdavgdata(:) = 1.0E12
ETrmaxdata(:) = 1.0E12
!indata_etmapdata(:) = 1.0E12
!indata_etravg(:) = 1.0E12
indata7(:,:) = 1.0E12

return
end subroutine alloc_global


! BROADCAST NAMELIST
!---------------------------------------------------------------

subroutine bc_namelist ()
use lycom_par
use mpi
use netcdf
use lycom_nc

implicit none

integer :: mperr
!call MPI_BARRIER(MPI)
call MPI_BCAST( year0,         1, MPI_INTEGER,          0, MPI_COMM_WORLD, mperr )
call MPI_BCAST( cyear0,         1, MPI_INTEGER,          0, MPI_COMM_WORLD, mperr )
call MPI_BCAST( tsindata0,         1, MPI_INTEGER,          0, MPI_COMM_WORLD, mperr )
call MPI_BCAST( accts0,         1, MPI_INTEGER,          0, MPI_COMM_WORLD, mperr )
call MPI_BCAST( tpos0,         1, MPI_INTEGER,          0, MPI_COMM_WORLD, mperr )
call MPI_BCAST( lastyear,       1, MPI_INTEGER,          0, MPI_COMM_WORLD, mperr )
call MPI_BCAST( runperiod,         1, MPI_INTEGER,          0, MPI_COMM_WORLD, mperr )
call MPI_BCAST( tsl,            1, MPI_INTEGER,          0, MPI_COMM_WORLD, mperr )
call MPI_BCAST( yearout1,       1, MPI_INTEGER,          0, MPI_COMM_WORLD, mperr )
call MPI_BCAST( yearoutX,       1, MPI_INTEGER,          0, MPI_COMM_WORLD, mperr )
call MPI_BCAST( outint,         1, MPI_INTEGER,          0, MPI_COMM_WORLD, mperr )
call MPI_BCAST( nSites,         1, MPI_INTEGER,          0, MPI_COMM_WORLD, mperr )
call MPI_BCAST( p_nspec,        1, MPI_INTEGER,          0, MPI_COMM_WORLD, mperr )
call MPI_BCAST( frac_s_init,  1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mperr )
call MPI_BCAST( specout,        1, MPI_LOGICAL,          0, MPI_COMM_WORLD, mperr )
call MPI_BCAST( interCan,       1, MPI_LOGICAL,          0, MPI_COMM_WORLD, mperr )
call MPI_BCAST( lrestart,        1, MPI_LOGICAL,          0, MPI_COMM_WORLD, mperr )
call MPI_BCAST( NOHONO,         1, MPI_LOGICAL,          0, MPI_COMM_WORLD, mperr )
call MPI_BCAST( BSCtypes,       1, MPI_LOGICAL,          0, MPI_COMM_WORLD, mperr )
call MPI_BCAST( noVeg,          1, MPI_LOGICAL,          0, MPI_COMM_WORLD, mperr )
call MPI_BCAST( inDirect,       1, MPI_INTEGER,          0, MPI_COMM_WORLD, mperr )

return
end subroutine bc_namelist


! BROADCAST SPECIES PARAMETERS
!---------------------------------------------------------------

subroutine bc_specpar ()
use lycom_par
use mpi
use netcdf
use lycom_nc

implicit none

integer :: mperr, k

do k = 1,p_nspecpar

  call MPI_BCAST( vec_o(:,k),  p_nspec, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mperr )
enddo

return
end subroutine bc_specpar


! READ BOUNDARY CONDITIONS (GLOBAL)
!----------------------------------------
subroutine lycom_readBCg ()
use lycom_par
use netcdf
use lycom_nc
implicit none

integer :: ncID, varID, i,k,l,m
!!DEBUG
integer :: nx2,ny2,nt2,lv
character (len=20)      :: dimName



!!DEBUG
!write(*,*) "read BC started"

! LAI
call check(nf90_open(laifileG, nf90_nowrite, ncID) )
call check(nf90_inq_varID(ncID, varName, varID))

!!DEBUG
!!call check(nf90_inquire_dimension(ncID, 1, dimName, nx2) )      !! TEMPORARY
!!call check(nf90_inquire_dimension(ncID, 2, dimName, ny2) )      !! TEMPORARY
!!call check(nf90_inquire_dimension(ncID, 3, dimName, nt2) )      !! TEMPORARY
!!write(*,*) "dim1  :",nx2  !12
!!write(*,*) "dim2  :",ny2  !144
!!write(*,*) "dim3  :",nt2  !143
!!write(*,*) "  "
!!write(*,*) "calling get_var with nx :",11, " and ny :", 11, "time 1"
!!
 !!call check(nf90_get_var(ncID, varID, bcdata0, start = (/ 1, 1, 1 /), count =
 !(/ 1, 11, 11 /) ))    !!  TEMPORARY
!!
!!write(*,*) "calling get_var with nx :",12, " and ny :", 12, "time 1"
!!
 !!call check(nf90_get_var(ncID, varID, bcdata0, start = (/ 1, 1, 1 /), count =
 !(/ 1, 12, 12 /) ))    !!  TEMPORARY
!!
!!!!write(*,*) "calling get_var with nx :",nx-1, " and ny :", ny-1, "time 1"
!!!!
!!!! call check(nf90_get_var(ncID, varID, bcdata0, start = (/ 1, 1, 1 /), count
!= (/ 1, nx-1, ny-1 /) ))    !!  TEMPORARY
!!
!!write(*,*) "calling get_var with nx :",nx-1, " and ny :", ny-1, "time 1 at POS
!3"
!!
 !!call check(nf90_get_var(ncID, varID, bcdata0, start = (/ 1, 1, 1 /), count =
 !(/ nx-1, ny-1, 1 /) ))    !!  TEMPORARY
!!
!!write(*,*) "calling get_var with nx :",nx, " and ny :", ny, "time 1 at POS 3"
!!
 !!call check(nf90_get_var(ncID, varID, bcdata0, start = (/ 1, 1, 1 /), count =
 !(/ nx, ny, 1 /) ))    !!  TEMPORARY
!!
!!!!write(*,*) "calling get_var with nx :",nx, " and ny :", ny, "time 1"
!!!!
!!!! call check(nf90_get_var(ncID, varID, bcdata0, start = (/ 1, 1, 1 /), count
!= (/ 1, nx, ny /) ))    !!  TEMPORARY
!!
!!write(*,*) "bcdata0 finished"
!!


!!call check(nf90_get_var(ncID, varID, bcdata))
do lv = 1,12
  call check(nf90_get_var(ncID, varID, bcdata(:,:,lv), start = (/ 1, 1, lv /), count = (/ nx, ny, 1 /) ))    !!  TEMPORARY
enddo

call check(nf90_close(ncID))

! assign boundary condition data to list of land points
do l=1,nland
  bcdataL(l,:) = bcdata(indvec(l,1),indvec(l,2),:)
enddo

! distribute boundary condition data to processors
do m = 1,12
  k=1
  do i = 1,numproc
    laidata((i-1)*pp+1:(i-1)*pp+ppvec(i), m) = bcdataL(k:k+ppvec(i)-1, m)

    k = k + ppvec(i)
  enddo
enddo


!!DEBUG
write(*,*) "read SAI started "


! SAI
call check(nf90_open(saifileG, nf90_nowrite, ncID) )
call check(nf90_inq_varID(ncID, varName, varID))
!!call check(nf90_get_var(ncID, varID, bcdata))
do lv = 1,12
  call check(nf90_get_var(ncID, varID, bcdata(:,:,lv), start = (/ 1, 1, lv /), count = (/ nx, ny, 1 /) ))    !!  TEMPORARY
enddo
call check(nf90_close(ncID))

! assign boundary condition data to list of land points
do l=1,nland
 bcdataL(l,:) = bcdata(indvec(l,1),indvec(l,2),:)
enddo

! distribute boundary condition data to processors
do m = 1,12
  k=1
  do i = 1,numproc
    saidata((i-1)*pp+1:(i-1)*pp+ppvec(i), m) = bcdataL(k:k+ppvec(i)-1, m)

    k = k + ppvec(i)
  enddo
enddo

! SSA
call check(nf90_open(landmask, nf90_nowrite, ncID) )
call check(nf90_inq_varID(ncID, varName, varID))
call check(nf90_get_var(ncID, varID, bcdata2))
write(*,*) "start SSA"
!call check(nf90_get_var(ncID, varID, bcdata2, start = (/ 1, 1, 1 /), count = (/ nx, ny, 1 /) ))    !!  TEMPORARY

call check(nf90_close(ncID))

! assign boundary condition data to list of land points
do l=1,nland
  bcdata2L(l) = bcdata2(indvec(l,1),indvec(l,2))
enddo

! distribute boundary condition data to processors
k=1
do i = 1,numproc
  ssadata((i-1)*pp+1:(i-1)*pp+ppvec(i)) = bcdata2L(k:k+ppvec(i)-1)

  k = k + ppvec(i)
enddo

! biomes
!call check(nf90_open(ssafileG, nf90_nowrite, ncID) )
!call check(nf90_inq_varID(ncID, varName, varID))
!call check(nf90_get_var(ncID, varID, bcdata2))
!write(*,*) "start biome"
!call check(nf90_get_var(ncID, varID, bcdata2, start = (/ 1, 1, 1 /), count = (/ nx, ny, 1 /) ))    !!  TEMPORARY
!call check(nf90_close(ncID))

! assign boundary condition data to list of land points
!do l=1,nland
!  bcdata2L(l) = bcdata2(indvec(l,1),indvec(l,2))
!enddo

! distribute boundary condition data to processors
!k=1
!do i = 1,numproc
!  biomedata((i-1)*pp+1:(i-1)*pp+ppvec(i)) = bcdata2L(k:k+ppvec(i)-1)

!  k = k + ppvec(i)
!enddo

call check(nf90_open(ETdavgfileG, nf90_nowrite, ncID) )
call check(nf90_inq_varID(ncID, varName, varID))
call check(nf90_get_var(ncID, varID, bcdata2))
write(*,*) "start ETDavg"
!call check(nf90_get_var(ncID, varID, bcdata2, start = (/ 1, 1, 1 /), count = (/ nx, ny, 1 /) ))    !!  TEMPORARY
call check(nf90_close(ncID))

! assign boundary condition data to list of land points
do l=1,nland
  bcdata2L(l) = bcdata2(indvec(l,1),indvec(l,2))
enddo

! distribute boundary condition data to processors
k=1
do i = 1,numproc
  ETdavgdata((i-1)*pp+1:(i-1)*pp+ppvec(i)) = bcdata2L(k:k+ppvec(i)-1)

  k = k + ppvec(i)
enddo

call check(nf90_open(ETrmaxfileG, nf90_nowrite, ncID) )
call check(nf90_inq_varID(ncID, varName, varID))
call check(nf90_get_var(ncID, varID, bcdata2))
write(*,*) "start etrmax"
!call check(nf90_get_var(ncID, varID, bcdata2, start = (/ 1, 1, 1 /), count = (/ nx, ny, 1 /) ))    !!  TEMPORARY
call check(nf90_close(ncID))

! assign boundary condition data to list of land points
do l=1,nland
  bcdata2L(l) = bcdata2(indvec(l,1),indvec(l,2))
enddo

! distribute boundary condition data to processors
k=1
do i = 1,numproc
  ETrmaxdata((i-1)*pp+1:(i-1)*pp+ppvec(i)) = bcdata2L(k:k+ppvec(i)-1)

  k = k + ppvec(i)
enddo

! biomes
call check(nf90_open(biomefileG, nf90_nowrite, ncID) )
call check(nf90_inq_varID(ncID, varName, varID))
call check(nf90_get_var(ncID, varID, bcdata2))
!write(*,*) "start biome"
!call check(nf90_get_var(ncID, varID, bcdata2, start = (/ 1, 1, 1 /), count = (/ nx, ny, 1 /) ))    !!  TEMPORARY
call check(nf90_close(ncID))

! assign boundary condition data to list of land points
do l=1,nland
  bcdata2L(l) = bcdata2(indvec(l,1),indvec(l,2))
enddo

! distribute boundary condition data to processors
k=1
do i = 1,numproc
  biomedata((i-1)*pp+1:(i-1)*pp+ppvec(i)) = bcdata2L(k:k+ppvec(i)-1)

  k = k + ppvec(i)
enddo
!write (*,*) "indvec 1 inside Boundary c", indvec(:,1)
!write (*,*) "indvec 2 inside Boundary c", indvec(:,2)

!!DEBUG
write(*,*) "read BC finished"

return
end subroutine lycom_readBCg

! SCATTER BOUNDARY CONDITIONS
!---------------------------------------------------------------

subroutine sc_BC ()
use lycom_par
use mpi
implicit none

integer :: m, mperr

real, allocatable :: laimon(:)
real, allocatable :: saimon(:)

allocate(laimon(pp))
allocate(saimon(pp))

do m = 1,12

  call MPI_BARRIER(MPI_COMM_WORLD, mperr)

! parts of v0 have length pp*8760, root=0
  call MPI_SCATTER( laidata(:,m), pp, MPI_DOUBLE_PRECISION, &
                    laimon, pp, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mperr)

  call MPI_SCATTER( saidata(:,m), pp, MPI_DOUBLE_PRECISION, &
                    saimon, pp, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mperr)

  areaLEAF(:,m) = laimon(:)
  areaSTEM(:,m) = saimon(:)
enddo

call MPI_BARRIER(MPI_COMM_WORLD, mperr)

call MPI_SCATTER( ssadata, pp, MPI_DOUBLE_PRECISION, areaSOIL, pp, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mperr)

call MPI_SCATTER( biomedata, pp, MPI_DOUBLE_PRECISION, biome, pp, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mperr)

call MPI_SCATTER( ETrmaxdata, pp, MPI_DOUBLE_PRECISION, ETrmax, pp, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mperr)

call MPI_SCATTER( ETdavgdata, pp, MPI_DOUBLE_PRECISION, ETdavg, pp, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mperr)

deallocate(laimon)
deallocate(saimon)

call MPI_BARRIER(MPI_COMM_WORLD, mperr)

return
end subroutine sc_BC


! OPEN GLOBAL FILES
!---------------------------------------------------------------

subroutine lycom_openG ()
use lycom_par
use netcdf
use lycom_nc
implicit none

! Open climate data files
call check(nf90_open(tairfileG, nf90_nowrite, ncID0(1)) )
call check(nf90_open(rhumfileG, nf90_nowrite, ncID0(2)) )
call check(nf90_open(windfileG, nf90_nowrite, ncID0(3)) )
call check(nf90_open(rainfileG, nf90_nowrite, ncID0(4)) )
call check(nf90_open(snowfileG, nf90_nowrite, ncID0(5)) )
call check(nf90_open(sradfileG, nf90_nowrite, ncID0(6)) )
call check(nf90_open(lradfileG, nf90_nowrite, ncID0(7)) )

return
end subroutine lycom_openG

! PREPROCESSING
!---------------------------------------------------------------
subroutine lycom_openPreprocessing () 
use lycom_par
use netcdf
use lycom_nc
use mpi
implicit none

integer :: varID, c,i,h,k,l,y
real :: ETd
integer :: mperr
!integer :: tstep_total          !! NEW
ETd=0.0
!#############################################################Commented everything from here
! Open climate data files
!call check(nf90_open(tairfileG, nf90_nowrite, ncID0(8)) )
!call check(nf90_open(sradfileG, nf90_nowrite, ncID0(9)) )
!call check(nf90_open(lradfileG, nf90_nowrite, ncID0(10)) )


!! NEW: loop over all timesteps of the input data. Easiest solution is to determine it by hand from the netcdf files, say tstep_total = 1234567
!! More complicated, but flexible solution: Copy the calendar calculation loop from libry_main, which determines the number of total
!hours based on the start year of the data, then you still need to provide the final year.
!! Would be also useful to make tstep_total or the final year a namelist parameter, same as start year (cyear0).

!do h = 1, tstep_total


! loop over all needed climate forcing files

!  do c = 8,10
!    call check(nf90_inq_varID(ncID0(c), varName, varID))
!    call check(nf90_get_var(ncID0(c),varID,indata(:,:), start = (/ 1, 1, h /), count = (/ nx, ny, 1 /) ))

!    ! assign boundary condition data to list of land points
!    do l = 1,nland
      
!      indataL(l) = indata(indvec(l,1),indvec(l,2))
!    enddo
    
!    if (c .eq. 8) then
!      tair_for_avg(:,h) = indataL(:)   
!    elseif (c .eq. 9) then               
!      srad_for_avg(:,h) = indataL(:)
!    elseif (c .eq. 10) then
!      lrad_for_avg(:,h) = indataL(:)
!    endif
!  enddo
!  !!Now you can use indataL(l) for your calculations. If you need to sum up something, you have to create new variables where you
!  !can accumulate indataL(l) { note that indataL is overwritten for each new climate variable, so you need to copy it to another
!  !variable }

  
!  !Summing not required
!  !accumulated_var1(:) = accumulated_var1(:) + tair_for_avg(:)    
!  !accumulated_var1(:) = accumulated_var1(:) + srad_for_avg(:)
!  !accumulated_var1(:) = accumulated_var1(:) + lrad_for_avg(:)

!enddo




!! After you have processed all time steps, you can carry out the calculations where you need the accumulated variables...
!!!!!!!!!!!!!!!!!!!!!!!!!
!do l = 1,nland

!  do y =1,tstep_total


    !write (*,*) srad_for_avg(l,y)
    !write (*,*) p_albVEG
    !write (*,*) lrad_for_avg(l,y)
    !write (*,*) tair_for_avg(l,y)
    !write (*,*) p_eps
    !write (*,*) c_sigma
!    R_net(l,y) = srad_for_avg(l,y)*(1.0-p_albVEG) + p_eps*lrad_for_avg(l,y) - p_eps*c_sigma*tair_for_avg(l,y)*tair_for_avg(l,y)*tair_for_avg(l,y)*tair_for_avg(l,y)
    !write (*,*) R_net(l,y)  
!    ta_C(l,y) = tair_for_avg(l,y) - c_TH2Osl

!    if (ta_C(l,y) .gt. 0.0) then

!      desat_dT3(l,y) = exp(p_esatAIR1 *ta_C(l,y) /(p_esatAIR2 + ta_C(l,y))) *p_esatAIR3 *p_esatAIR1 * p_esatAIR2 /(p_esatAIR2 +ta_C(l,y))**2

!      ET_pot0(l,y) = 1.3 * R_net(l,y) * desat_dT3(l,y) / ( desat_dT3(l,y) + c_gamma) / c_HH2Olg / c_rhoH2Ol
!    else

!      ET_pot0(l,y) = 0.0
!    endif

!    ETd = ETd + ET_pot0(l,y)
!  enddo
!  E_Tdavg(l) = ETd / float(tstep_total)
  
!enddo

!write (*,*) "The ETdavg:preprocessing ", E_Tdavg(l)

!do l = 1,nland

!  do y = 1,tstep_total                              

!    R_net(l,y) = srad_for_avg(l,y)*(1.0-p_albVEG) + p_eps*lrad_for_avg(l,y) - p_eps*c_sigma*tair_for_avg(l,y)*tair_for_avg(l,y)*tair_for_avg(l,y)*tair_for_avg(l,y)

!    ta_C(l,y) = tair_for_avg(l,y) - c_TH2Osl

!    if (ta_C(l,y) .gt. 0.0) then

!      desat_dT3(l,y) = exp(p_esatAIR1 *ta_C(l,y) /(p_esatAIR2 + ta_C(l,y))) *p_esatAIR3 *p_esatAIR1 * p_esatAIR2 /(p_esatAIR2 +ta_C(l,y))**2

!      ET_pot0(l,y) = 1.3 * R_net(l,y) * desat_dT3(l,y) / ( desat_dT3(l,y) + c_gamma) / c_HH2Olg / c_rhoH2Ol
!    else

!      ET_pot0(l,y) = 0.0
!    endif
!    E_Tratio(l,y) = ET_pot0(l,y) / E_Tdavg(l)

!    E_Trmax(l)= max(E_Trmax(l), E_Tratio(l,y))
  
!  enddo
  
!enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! distribute boundary condition data to processors
!k=1
!do i = 1,numproc
!  indata_etmapdata((i-1)*pp+1:(i-1)*pp+ppvec(i)) = E_Trmax(k:k+ppvec(i)-1)  
!  k = k + ppvec(i)
  
!enddo
!k=1
!do i = 1,numproc
!  indata_etravg((i-1)*pp+1:(i-1)*pp+ppvec(i)) = E_Tdavg(k:k+ppvec(i)-1)  
!  k = k + ppvec(i)
  
!enddo
!write (*,*) "The preprocessing has been completed"

!call MPI_BARRIER(MPI_COMM_WORLD, mperr)


!call MPI_SCATTER( indata_etmapdata(:), pp, MPI_DOUBLE_PRECISION, Etr_max, pp, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mperr)
!call MPI_SCATTER( indata_etravg(:), pp, MPI_DOUBLE_PRECISION, Etr_avg, pp, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mperr)


!call check(nf90_close(ncID0(8)))
!call check(nf90_close(ncID0(9)))
!call check(nf90_close(ncID0(10)))

!call MPI_BARRIER(MPI_COMM_WORLD, mperr)

return
end subroutine lycom_openPreprocessing 



! READ GLOBAL INPUT DATA FILES
!---------------------------------------------------------------

subroutine lycom_readHourG (h)
use lycom_par
use netcdf
use lycom_nc
implicit none

integer :: varID, c,i,h,k,l
integer :: counter
counter=1
!write(*,*) "readHourG is entered"
! loop over all 7 climate forcing files
!write(*,*) "reading Hour G", h
!write(*,*) "Proc: in readHourG ",rank, " indv",size(indvec)
!write(*,*)"indvec 1 ", indvec(:,1)
!write(*,*)"indvec 2 ", indvec(:,2)
do c = 1,7
  !write(*,*) "c is",c
  call check(nf90_inq_varID(ncID0(c), varName, varID))
  call check(nf90_get_var(ncID0(c),varID,indata(:,:), start = (/ 1, 1, h /), count = (/ nx, ny, 1 /) ))

  !write(*,*) "indata is", indata(:,:)
  ! assign boundary condition data to list of land points
  do l = 1,nland
    !write(*,*)"indvec 1 ",indvec(l,1)
    !write(*,*)"indvec 2 ",indvec(l,2)
    !write(*,*)"l here is ",l
    indataL(l) = indata(indvec1(l,1),indvec1(l,2))
    !write(*,*) "counter is", indata(:,:)
    !if (indvec(l,1) .lt. 250) then
        !counter = counter +1
        !write(*,*) "counter is", counter
    !endif 
  enddo

  ! distribute boundary condition data to processors
  k=1
  do i = 1,numproc
    indata7((i-1)*pp+1:(i-1)*pp+ppvec(i), c) = indataL(k:k+ppvec(i)-1)

    k = k + ppvec(i)
  enddo
enddo
!write(*,*) "reading for h is done"
return
end subroutine lycom_readHourG


! SCATTER GLOBAL INPUT DATA
!---------------------------------------------------------------

subroutine sc_clim ()
use lycom_par
use mpi
implicit none

integer :: mperr

call MPI_BARRIER(MPI_COMM_WORLD, mperr)

call MPI_SCATTER( indata7(:,1), pp, MPI_DOUBLE_PRECISION, &
                  xT_a, pp, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mperr)

call MPI_SCATTER( indata7(:,2), pp, MPI_DOUBLE_PRECISION, &
                  rH2Og_RH, pp, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mperr)

call MPI_SCATTER( indata7(:,3), pp, MPI_DOUBLE_PRECISION, &
                  fAIR_s, pp, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mperr)

call MPI_SCATTER( indata7(:,4), pp, MPI_DOUBLE_PRECISION, &
                  fH2Ol_ad, pp, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mperr)

call MPI_SCATTER( indata7(:,5), pp, MPI_DOUBLE_PRECISION, &
                  fH2Os_ad, pp, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mperr)

call MPI_SCATTER( indata7(:,6), pp, MPI_DOUBLE_PRECISION, &
                  fRADs_ad, pp, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mperr)

call MPI_SCATTER( indata7(:,7), pp, MPI_DOUBLE_PRECISION, &
                  fRADl_ad, pp, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mperr)

call MPI_BARRIER(MPI_COMM_WORLD, mperr)

return
end subroutine sc_clim


! GLOBAL OUTPUT
!---------------------------------------------------------------

subroutine lycom_outputG ()
use lycom_par
use lycom_global2
use netcdf
use lycom_nc
implicit none

integer :: x_dimID, y_dimID, t_dimID, k
integer :: x_varID, y_varID

if (rank .eq. 0) then

  if (out1) then

    write(*,*) "Creating output file:", trim(sfile_outputG)
    ! Create global output file
    call check(nf90_create(sfile_outputG, NF90_CLOBBER, outID))
    write(*,*) "Output file created with ID:", outID   
    ! Define dimensions.
    call check(nf90_def_dim(outID, "lon", nx, x_dimID))    ! returns dimID
    call check(nf90_def_dim(outID, "lat", ny, y_dimID))
    call check(nf90_def_dim(outID, "time", NF90_UNLIMITED, t_dimID))
    dimIDs3 = (/ x_dimID, y_dimID, t_dimID /)
    
    !write(*,*) "nx Position is", nx
    !write(*,*) "ny Position is ", ny
    write(*,*) "Dimensions defined:"
    write(*,*) "  lon (x_dimID):", x_dimID, " size:", nx
    write(*,*) "  lat (y_dimID):", y_dimID, " size:", ny
    write(*,*) "  time (t_dimID):", t_dimID, " (unlimited)"
    ! Define coordinate variables
    call check(nf90_def_var(outID, "lon", NF90_REAL, x_dimID, x_varID))    ! returns varID
    call check(nf90_put_att(outID, x_varID, "long_name", "longitude"))
    call check(nf90_put_att(outID, x_varID, "units", "degrees_east"))
    call check(nf90_put_att(outID, x_varID, "standard_name", "longitude"))
    call check(nf90_def_var(outID, "lat", NF90_REAL, y_dimID, y_varID))
    call check(nf90_put_att(outID, y_varID, "long_name", "latitude"))
    call check(nf90_put_att(outID, y_varID, "units", "degrees_north"))
    call check(nf90_put_att(outID, y_varID, "standard_name", "latitude"))
   
    call check(nf90_def_var(outID, "time", NF90_INT, t_dimID, t_varID))
    call check(nf90_put_att(outID, t_varID, "units", "day as %Y%m%d.%f"))
    call check(nf90_put_att(outID, t_varID, "calendar", "proleptic_gregorian"))

    ! Define output variables

    do k = 1,1
      call def_varG(k*1000+kfile_rCO2d          ,  outvarID(1+25*(k-1)))
!      call def_varG(k*1000+kfile_sCO2d          ,  outvarID(2+25*(k-1)))
      call def_varG(k*1000+kfile_rCb            ,  outvarID(3+25*(k-1)))
      call def_varG(k*1000+kfile_rH2Ol          ,  outvarID(4+25*(k-1)))
      call def_varG(k*1000+kfile_rmaxH2Ol       ,  outvarID(5+25*(k-1)))
      call def_varG(k*1000+kfile_area           ,  outvarID(6+25*(k-1)))
      call def_varG(k*1000+kfile_act            ,  outvarID(7+25*(k-1)))
      call def_varG(k*1000+kfile_fCO2gc         ,  outvarID(8+25*(k-1)))
      call def_varG(k*1000+kfile_fCcg           ,  outvarID(9+25*(k-1)))
      call def_varG(k*1000+kfile_fCcb           ,  outvarID(10+25*(k-1)))
      call def_varG(k*1000+kfile_fCbo           , outvarID(11+25*(k-1)))
      call def_varG(k*1000+kfile_fH2Ol_ux       , outvarID(12+25*(k-1)))
      !call def_varG(k*1000+kfile_fH2Ol_lsat     , outvarID(13+25*(k-1)))
      !call def_varG(k*1000+kfile_fH2Ol_bsat     , outvarID(14+25*(k-1)))
      call def_varG(k*1000+kfile_fH2Ol_runoff_l , outvarID(15+25*(k-1)))
      call def_varG(k*1000+kfile_fCc_gpp        , outvarID(16+25*(k-1)))
      call def_varG(k*1000+kfile_fCc_npp        , outvarID(17+25*(k-1)))      
      call def_varG(k*1000+kfile_fH2Ol_xd       , outvarID(18+25*(k-1)))
      call def_varG(k*1000+kfile_fH2Ogl_ux      , outvarID(19+25*(k-1)))
      call def_varG(k*1000+kfile_fH2Olg_xu      , outvarID(20+25*(k-1)))
      call def_varG(k*1000+kfile_Ts             , outvarID(21+25*(k-1)))
!      call def_varG(k*1000+kfile_Ts_N           , outvarID(16+20*(k-1)))
      call def_varG(k*1000+kfile_H              , outvarID(22+25*(k-1)))
      call def_varG(k*1000+kfile_E              , outvarID(23+25*(k-1)))
      call def_varG(k*1000+kfile_C              , outvarID(24+25*(k-1)))
      call def_varG(k*1000+kfile_EB             , outvarID(25+25*(k-1)))
    enddo

    call def_varG(3000+kfile_rH2Ol_g1           , outvarID(100)) ! must be larger than highest index of outvarID 3 lines above
    call def_varG(3000+kfile_rH2Ol_g2           , outvarID(101))
    call def_varG(3000+kfile_rH2Os_g            , outvarID(102))
    call def_varG(3000+kfile_fH2Ol_ug           , outvarID(103))
    call def_varG(3000+kfile_fH2Ol_ug2          , outvarID(104))
    call def_varG(3000+kfile_fH2Ol_go           , outvarID(105))
    call def_varG(3000+kfile_fH2Ol_gb           , outvarID(106))
    call def_varG(3000+kfile_fH2Olg_ga          , outvarID(107))
    call def_varG(3000+kfile_fH2Osl_g           , outvarID(108))
    call def_varG(3000+kfile_fH2Os_ad           , outvarID(109))
    call def_varG(3000+kfile_fH2Ol_ad           , outvarID(110))
    call def_varG(3000+kfile_xT_a               , outvarID(111))
    call def_varG(3000+kfile_Tg                 , outvarID(112)) 
    call def_varG(3000+kfile_G                  , outvarID(113))
    call def_varG(3000+kfile_fRADs              , outvarID(114))
    call def_varG(3000+kfile_sCO2d              , outvarID(115))
    call def_varG(3000+kfile_Lai                , outvarID(116))
    ! Optional variables
    if (BSCtypes) then
      call def_varG(kfile_MareaLCgS     , outvarID2(1))                
      call def_varG(kfile_MareaDCgS     , outvarID2(2))                
      call def_varG(kfile_MareaCCgS     , outvarID2(3))
      call def_varG(kfile_MareaMCgS     , outvarID2(4))
      call def_varG(kfile_MrH2OlLC_S    , outvarID2(5))
      call def_varG(kfile_MrH2OlDC_S    , outvarID2(6))
      call def_varG(kfile_MrH2OlCC_S    , outvarID2(7))
      call def_varG(kfile_MrH2OlMC_S    , outvarID2(8))
      call def_varG(kfile_MTsLC_S       , outvarID2(9))
      call def_varG(kfile_MTsDC_S       , outvarID2(10))
      call def_varG(kfile_MTsCC_S       , outvarID2(11))
      call def_varG(kfile_MTsMC_S       , outvarID2(12))
      call def_varG(kfile_MfCcbLC_S     , outvarID2(13))
      call def_varG(kfile_MfCcbDC_S     , outvarID2(14))
      call def_varG(kfile_MfCcbCC_S     , outvarID2(15))
      call def_varG(kfile_MfCcbMC_S     , outvarID2(16))

      if (NOHONO) then

        call def_varG(kfile_MfNO_N      , outvarID2(20))
        call def_varG(kfile_MfHONO_N    , outvarID2(21))
      endif
    endif

    call check(nf90_enddef(outID)) !End Definitions  

    ! Write grid variables
    call check(nf90_put_var(outID, x_varID, xpos))
    call check(nf90_put_var(outID, y_varID, ypos))
  else

    ! Open output file for next write
    call check(nf90_open(sfile_outputG, NF90_WRITE, outID))
  endif
endif

! Write output (variable, code, outfileID, dimensions)

do k = 1,1
  call write_varG(ag_rCO2d(:,k)     ,  outvarID(1+25*(k-1)))
!  call write_varG(ag_sCO2d(:,k)     ,  outvarID(2+25*(k-1))) !commented out
  call write_varG(ag_rCb(:,k)       ,  outvarID(3+25*(k-1)))
  call write_varG(ag_rH2Ol(:,k)     ,  outvarID(4+25*(k-1)))
  call write_varG(ag_rmaxH2Ol(:,k)  ,  outvarID(5+25*(k-1)))
  call write_varG(ag_area_s(:,k)  ,  outvarID(6+25*(k-1)))
  call write_varG(ag_act(:,k)       ,  outvarID(7+25*(k-1)))
  call write_varG(ag_fCO2gc(:,k)    ,  outvarID(8+25*(k-1)))
  call write_varG(ag_fCcg(:,k)      ,  outvarID(9+25*(k-1)))
  call write_varG(ag_fCcb(:,k)      ,  outvarID(10+25*(k-1)))
  call write_varG(ag_fCbo(:,k)      , outvarID(11+25*(k-1)))
  call write_varG(ag_fH2Ol_ux(:,k)  , outvarID(12+25*(k-1)))
  !call write_varG(ag_fH2Ol_lsat(:,k), outvarID(13+25*(k-1)))
  !call write_varG(ag_fH2Ol_bsat(:,k), outvarID(14+25*(k-1)))
  call write_varG(ag_fH2Ol_runoff_l(:,k), outvarID(15+25*(k-1)))
  call write_varG(ag_fCc_gpp(:,k)   , outvarID(16+25*(k-1)))
  call write_varG(ag_fCc_npp(:,k)   , outvarID(17+25*(k-1)))
  call write_varG(ag_fH2Ol_xd(:,k)  , outvarID(18+25*(k-1)))
  
  call write_varG(ag_fH2Ogl_ux(:,k) , outvarID(19+25*(k-1)))
  call write_varG(ag_fH2Olg_xu(:,k) , outvarID(20+25*(k-1)))
  call write_varG(ag_Ts(:,k)        , outvarID(21+25*(k-1)))
!  call write_varG(ag_Ts_N(:,k)      , outvarID(16+20*(k-1)))
  call write_varG(ag_H(:,k)         , outvarID(22+25*(k-1)))
  call write_varG(ag_E(:,k)         , outvarID(23+25*(k-1)))
  call write_varG(ag_C(:,k)         , outvarID(24+25*(k-1)))
  call write_varG(ag_EB(:,k)        , outvarID(25+25*(k-1)))
enddo
call write_varG(ag_Lai(:)             , outvarID(116))
call write_varG(ag_sCO2d(:)           , outvarID(115))
call write_varG(ag_rH2Ol_g1(:)        , outvarID(100))
call write_varG(ag_rH2Ol_g2(:)        , outvarID(101))
call write_varG(ag_rH2Os_g(:)         , outvarID(102))
call write_varG(ag_fH2Ol_ug(:)        , outvarID(103))
call write_varG(ag_fH2Ol_ug2(:)       , outvarID(104))
call write_varG(ag_fH2Ol_go(:)        , outvarID(105))
call write_varG(ag_fH2Ol_gb(:)        , outvarID(106))
call write_varG(ag_fH2Olg_ga(:)       , outvarID(107))
call write_varG(ag_fH2Osl_g(:)        , outvarID(108))
call write_varG(ag_fH2Os_ad(:)        , outvarID(109))
call write_varG(ag_fH2Ol_ad(:)        , outvarID(110))
call write_varG(ag_xT_a(:)            , outvarID(111))
call write_varG(ag_Tg(:)              , outvarID(112)) 
call write_varG(ag_G(:)               , outvarID(113))
call write_varG(ag_fRADs(:)           , outvarID(114))

! Optional output variables
if (BSCtypes) then
  call write_varG(a_MareaTHLC_g_S(:), outvarID2(1))                 
  call write_varG(a_MareaTHDC_g_S(:), outvarID2(2))                 
  call write_varG(a_MareaTHCC_g_S(:), outvarID2(3))
  call write_varG(a_MareaTHMC_g_S(:), outvarID2(4))
  call write_varG(a_MrH2OlLC_S(:) , outvarID2(5))
  call write_varG(a_MrH2OlDC_S(:) , outvarID2(6))
  call write_varG(a_MrH2OlCC_S(:) , outvarID2(7))
  call write_varG(a_MrH2OlMC_S(:) , outvarID2(8))
  call write_varG(a_MTsLC_S(:)    , outvarID2(9))
  call write_varG(a_MTsDC_S(:)    , outvarID2(10))
  call write_varG(a_MTsCC_S(:)    , outvarID2(11))
  call write_varG(a_MTsMC_S(:)    , outvarID2(12))
  call write_varG(a_MfCcbLC_S(:)  , outvarID2(13))
  call write_varG(a_MfCcbDC_S(:)  , outvarID2(14))
  call write_varG(a_MfCcbCC_S(:)  , outvarID2(15))
  call write_varG(a_MfCcbMC_S(:)  , outvarID2(16))

  if (NOHONO) then

    call write_varG(a_MfNO_N(:)   , outvarID2(20))
    call write_varG(a_MfHONO_N(:) , outvarID2(21))
  endif
endif

if (rank .eq. 0) call check(nf90_close(outID))

return
end subroutine lycom_outputG


! WRITE SPECIES PROPERTIES
!---------------------------------------------------------------

subroutine lycom_outspecG ()
use lycom_par
use lycom_global2
use netcdf
use lycom_nc
implicit none

integer :: x_dimID, y_dimID, t_dimID
integer :: x_varID, y_varID

integer :: k,l,i0
integer :: kfile_o_spec, kfile_habitat

! Open output file
if (rank .eq. 0) then

!  if (out1) then

    ! Create global output file
    call check(nf90_create(sfile_outputGS, NF90_CLOBBER, outIDS))

    ! Define dimensions.
    call check(nf90_def_dim(outIDS, "lon", nx, x_dimID))    ! returns dimID
    call check(nf90_def_dim(outIDS, "lat", ny, y_dimID))
    call check(nf90_def_dim(outIDS, "time", NF90_UNLIMITED, t_dimID))
    dimIDs3 = (/ x_dimID, y_dimID, t_dimID /)
   
    ! Define coordinate variables
    call check(nf90_def_var(outIDS, "lon", NF90_REAL, x_dimID, x_varID))    ! returns varID
    call check(nf90_put_att(outIDS, x_varID, "long_name", "longitude"))
    call check(nf90_put_att(outIDS, x_varID, "units", "degrees_east"))
    call check(nf90_put_att(outIDS, x_varID, "standard_name", "longitude"))
    call check(nf90_def_var(outIDS, "lat", NF90_REAL, y_dimID, y_varID))
    call check(nf90_put_att(outIDS, y_varID, "long_name", "latitude"))
    call check(nf90_put_att(outIDS, y_varID, "units", "degrees_north"))
    call check(nf90_put_att(outIDS, y_varID, "standard_name", "latitude"))
   
    call check(nf90_def_var(outIDS, "time", NF90_INT, t_dimID, t_varID))
    call check(nf90_put_att(outIDS, t_varID, "units", "day as %Y%m%d.%f"))
    call check(nf90_put_att(outIDS, t_varID, "calendar", "proleptic_gregorian"))

    ! Define output variables
    call def_varG(kfile_count_spec, outvarIDS(1))

    i0 = 0

    do l = 1,p_nhabA

      kfile_o_spec = kfile_count_spec + l * 100

      call def_varG(kfile_o_spec, outvarIDS(1+l))

      do k = 1,p_nspecpar

        kfile_o_spec = kfile_count_spec + l * 100 + k

        i0 = i0 + 1

        call def_varG(kfile_o_spec, outvarIDS(1+p_nhabA+i0))
      enddo
    enddo

    call check(nf90_enddef(outIDS)) !End Definitions  

    ! Write grid variables
    call check(nf90_put_var(outIDS, x_varID, xpos))
    call check(nf90_put_var(outIDS, y_varID, ypos))

!  else
!
!    ! Open output file for next write
!    call check(nf90_open(sfile_outputGS, NF90_WRITE, outIDS))
!  endif
endif

! write species properties
call write_varG(count_spec(:), outvarIDS(1))

i0 = 0

do l = 1,p_nhabA

  call write_varG(count_spec_h(l,:), outvarIDS(1+l))

  do k = 1,p_nspecpar

    i0 = i0 + 1

    call write_varG(vec_o_avg(k,l,:), outvarIDS(1+p_nhabA+i0))
  enddo
enddo

if (rank .eq. 0) call check(nf90_close(outIDS))

return
end subroutine lycom_outspecG


! CLOSE GLOBAL FILES
!---------------------------------------------------------------

subroutine lycom_closeG ()
use lycom_par
use netcdf
use lycom_nc
implicit none

! Close climate data files
call check(nf90_close(ncID0(1)))
call check(nf90_close(ncID0(2)))
call check(nf90_close(ncID0(3)))
call check(nf90_close(ncID0(4)))
call check(nf90_close(ncID0(5)))
call check(nf90_close(ncID0(6)))
call check(nf90_close(ncID0(7)))

return
end subroutine lycom_closeG


! DEALLOCATE GLOBAL FIELDS
!---------------------------------------------------------------

subroutine dealloc_global ()
use lycom_par
use mpi
implicit none

integer                 :: mperr

if (rank .eq. 0) then
  deallocate( indata, &
  indataL, &
  bcdata, &
  bcdataL, &
  bcdata2, &
  bcdata2L, &
 ! R_net, &
 ! ta_C, &
 ! desat_dT3, &
 ! ET_pot0, &
 ! E_Tdavg, &
 ! E_Trmax, &
 ! E_Tratio, &
  xpos, &
  ypos, &
  indvec )
endif

write(*,*) "Proc:  ",rank, "  LAI",size(laidata)
write(*,*) "Proc:  ",rank, "  SAI",size(saidata)
write(*,*) "Proc:  ",rank, "  SSA",size(ssadata)
write(*,*) "Proc:  ",rank, "biome",size(biomedata)
!write(*,*) "Proc:  ",rank, "  IN7",size(indata7)
write(*,*) "Proc:  ",rank, "ppvec",size(ppvec)
!write(*,*) "Proc:  ",rank, " xpos",size(xpos)
!write(*,*) "Proc:  ",rank, " ypos",size(ypos)
!write(*,*) "Proc:  ",rank, " indv",size(indvec)
!write(*,*) "       "

deallocate( laidata, &
saidata, &
ssadata, &
biomedata, &
indata7, &
!indata_etmapdata, &
!indata_etravg, &
ETdavgdata, &
ETrmaxdata, &
ppvec )

!DEBUG
 write( kstatus,* ) "MPI_FINALIZE started"

! Stop MPI
call MPI_FINALIZE(mperr)

!DEBUG
 write( kstatus,* ) "MPI_FINALIZE status:", mperr

return
end subroutine dealloc_global

end module lycom_global

