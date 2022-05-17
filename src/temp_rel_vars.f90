!
!    Author:    Jacopo Simoni
!    module file containing variables for the
!    temperature relaxation calculation
!
! --------------------------------------------------------------------
!
MODULE temp_rel_vars
  
  !
  SAVE
  
  character (LEN=256)        :: temprel_data_dir
  !  dir. with 
  !  rel. matr. elem.
  character (LEN=256)        :: tensor_fil
  !             
  ! file with tensor data
  character (LEN=256)        :: temprel_fil
  !  file with 
  !  temp. rel. data
  character (LEN=256)        :: sprel_fil
  !  file with
  !  spin relax. data
  character (LEN=256)        :: filvscr
  !  file with
  !  scr. potential
  integer                    :: st_at
  !  starting atom 
  !  calculation
  integer                    :: fin_at
  !  final atom
  !  calculation
  
  
  
END MODULE temp_rel_vars
!
!
