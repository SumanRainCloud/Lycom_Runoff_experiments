
! INITIALISE MPI
!---------------------------------------------------------------

subroutine init_mpi ()
use lycom_par
implicit none

rank = 0

para = .false.

return
end subroutine init_mpi

! Dummy lycom_GLOBAL*
!---------------------------------------------------------------

module lycom_global2

end module lycom_global2

!---------------------------------------------------------------

module lycom_global
contains

! other subroutines
!---------------------------------------------------------------
subroutine read_land ()
end subroutine read_land

!---------------------------------------------------------------
subroutine alloc_global ()
end subroutine alloc_global

!---------------------------------------------------------------
subroutine bc_pp ()
end subroutine bc_pp

!---------------------------------------------------------------
subroutine bc_namelist ()
end subroutine bc_namelist

!---------------------------------------------------------------
subroutine bc_specpar ()
end subroutine bc_specpar

!---------------------------------------------------------------
subroutine lycom_readBCg ()
end subroutine lycom_readBCg

!---------------------------------------------------------------
subroutine sc_BC ()
end subroutine sc_BC

!---------------------------------------------------------------
subroutine lycom_openG ()
end subroutine lycom_openG 

!---------------------------------------------------------------
subroutine lycom_openPreprocessing ()
end subroutine lycom_openPreprocessing 

!---------------------------------------------------------------
subroutine lycom_readHourG (dummy)
  integer :: dummy
end subroutine lycom_readHourG

!---------------------------------------------------------------
subroutine sc_clim ()
end subroutine sc_clim 

!---------------------------------------------------------------
subroutine lycom_outputG ()
end subroutine lycom_outputG 

!---------------------------------------------------------------
subroutine def_varG ()
end subroutine def_varG 

!---------------------------------------------------------------
subroutine write_varG ()
end subroutine write_varG 

!---------------------------------------------------------------
subroutine lycom_outspecG ()
end subroutine lycom_outspecG 

!---------------------------------------------------------------
subroutine lycom_closeG ()
end subroutine lycom_closeG 

!---------------------------------------------------------------
subroutine dealloc_global ()
end subroutine dealloc_global 

end module lycom_global

