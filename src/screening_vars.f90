!
!    author:   Jacopo Simoni
!    module file with definition of variables for dynamical calc.
!
! ---------------------------------------------------------------------------
!
!
MODULE dynvars
  !
  USE kinds,      ONLY : DP
  
  !
  !   ... dynamical variables used in GW
  !   and phonon calculations
  !
  
  SAVE
  !
  
  LOGICAL             :: dyncalc
  !   dynamic / static calculation
  real(DP)            :: freq_broadening
  !   i\eta freq. broadening
  integer             :: nfreq
  !   number of frequencies / dynamic calc
  real(DP)            :: en_max
  !   max. energy (eV) -> chi(w) / eps(w)
  !
END MODULE dynvars
!
!
MODULE screening_vars
  !
  USE kinds,          ONLY : DP
  
  SAVE
  !
  real(DP)                 :: t_ev                ! elec. temperature in eV
  real(DP)                 :: n_elec              ! elec. density (cm^-3)
  real(DP)                 :: n_carrier           ! carrier density (cm^-3)
  real(DP)                 :: Nc_n                ! Nc / n -> SC
  real(DP)                 :: Nv_p                ! Nv / p -> SC
  real(DP)                 :: n_elec_up           ! elec. density (cm^-3)
                                                  ! up electrons
  real(DP)                 :: n_elec_dn           ! dn electrons
  character(LEN=256)       :: scr_type            ! model of screening
  real(DP), parameter      :: relerr    = 1.e-12
  
  !
END MODULE screening_vars
!
