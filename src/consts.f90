!
!  author:  Jacopo Simoni
!  module file containing additional constants
!
!  ----------------------------------------------------------
!
MODULE consts
  !
  USE kinds,        ONLY : DP
  !
  SAVE
  !
  real(DP), parameter             :: AngtoBohr = 1.889725989
  real(DP), parameter             :: e_4pieps0 = 1.439803E-9   ! SI UNITS
  real(DP), parameter             :: elec_mass = 0.056856      ! (eV fs^2 / Ang^2)
  real(DP), parameter             :: hbar = 0.6582119514       ! (eV fs)
  real(DP), parameter             :: proton_mass = 29.23394467 ! (eV fs^2 / bohr^2)
  !
END MODULE consts
