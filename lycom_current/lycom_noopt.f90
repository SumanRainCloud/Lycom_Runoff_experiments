
module lycom_opt
contains

! DEFINE BSC TYPES
!---------------------------------------------------------------

subroutine lycom_defBSC (k)
integer :: k
end subroutine lycom_defBSC

! SUM UP BSC COVER
!---------------------------------------------------------------

subroutine lycom_covsBSC (k,k2)
integer :: k,k2
end subroutine lycom_covsBSC

! READ NOHONO FLUXES
!---------------------------------------------------------------
subroutine lycom_readNOHONO ()
end subroutine lycom_readNOHONO

! BROADCAST NOHONO DATA
!---------------------------------------------------------------
subroutine bc_NOHONO ()
end subroutine bc_NOHONO

! ALLOCATE BSC AND NO/HONO VARIABLES
!---------------------------------------------------------------

subroutine lycom_allocBSC ()
end subroutine lycom_allocBSC

! DEALLOCATE BSC AND NO/HONO VARIABLES
!---------------------------------------------------------------

subroutine lycom_deallocBSC ()
end subroutine lycom_deallocBSC

! CALCULATE NO AND HONO EMISSIONS
!---------------------------------------------------------------

subroutine lycom_fNOHONO (k,st,ac)
integer :: k
real    :: st,ac
end subroutine lycom_fNOHONO

! ACCUMULATE BSC PROPERTIES
!---------------------------------------------------------------

subroutine lycom_accBSC (i,i2,j)
integer :: i,i2,j
end subroutine lycom_accBSC

! AVERAGE BSC PROPERTIES
!---------------------------------------------------------------

subroutine lycom_avBSC (i)
integer :: i
end subroutine lycom_avBSC

! AVERAGE GRID CELL BSC PROPERTIES
!---------------------------------------------------------------

subroutine lycom_avgcBSC (i,i2,j)
integer :: i,i2,j
end subroutine lycom_avgcBSC

end module lycom_opt

