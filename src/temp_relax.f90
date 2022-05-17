! ----------------------------------------------------------------------------
!
!      THIS PROGRAM COMPUTES THE TEMPERATURE RELAXATION
!      OF A PLASMA OR CRYSTALLINE SYSTEM AT FINITE TEMPERATURE
!
!      the program reads prefix.wfc written as a QE output, it computes
!      the screened potential gradient around each atom and
!      evaluate the temperature relaxation rate g(T)
!
!      AUTHOR:     Jacopo Simoni
!
!      atomic potential supported:
!      1) PAW
!      2) Norm-conserving - local
!      3) AE (never used)
!
!      input file:
!
!      &inputph
!         prefix='pwscf',
!         outdir='./tmp',
!         fil_chi='al.chi',                 # file with KS suscept. (input)
!         temprel_fil='al.temprel.dat',     # file with temp. rel. data
!         screening_model='lindhard',
!         elec_temperature=0.1,
!         temprel_data_dir='./tmp/mtxel/',
!         compute_forces=.true.,            # if true compute forces else jump to final section
!         filrhoc='al.rhoc'
!         filrhoc='al.rhov'
!         filvscr='al.vscr'
!         st_at=1,                          # initial at. to compute
!         fin_at=nat,                       # final at. to compute
!      /
! -----------------------------------------------------------------------------

PROGRAM temp_relax
  
  !
  !    MODULES TO LOAD
  !
  
  USE kinds,                  ONLY : DP
  USE constants,              ONLY : eps8, K_BOLTZMANN_SI, pi, RYTOEV
  USE consts,                 ONLY : hbar, proton_mass
  USE mp_images,              ONLY : intra_image_comm, me_image
  USE mp_pools,               ONLY : npool
  USE environment,            ONLY : environment_start, environment_end
  USE mp,                     ONLY : mp_bcast, mp_barrier, mp_sum
  USE mp_global,              ONLY : mp_startup
  USE mp_world,               ONLY : world_comm
  USE io_global,              ONLY : ionode, ionode_id
  USE io_files,               ONLY : prefix, tmp_dir
  USE fft_base,               ONLY : dffts
  USE ks_screening,           ONLY : fil_chi
  USE ions_base,              ONLY : amass, nat, ityp
  USE wvfct,                  ONLY : npwx, nbnd, et
  USE klist,                  ONLY : nelec, nks, wk
  USE temp_rel_vars,          ONLY : st_at, fin_at, temprel_fil,        &
       temprel_data_dir, filvscr
  USE ph_lt_vars,             ONLY : M, t, ni, ntd, t_max
  USE charge_density,         ONLY : calc_eff_charge, filrhoc, filrhov
  USE transitions,            ONLY : trans_list, trans_list_auxp
  USE uspp_param,             ONLY : upf
  USE force_mtxel,            ONLY : ngr, psi_gradv_psi, trans_groups, set_force_transitions,    &
       set_fmtxel_smoothg, read_matrxel_file, trans_list_full, set_fmtxel_rec_space, dviapsi
  USE screening_vars,         ONLY : t_ev, n_elec, scr_type
  USE apawpot,                ONLY : set_paw_dvie_dr, set_dvie_dr
  USE onsager_vars,           ONLY : pocc
  USE noncollin_module,       ONLY : npol
  
  implicit none
  
  !
  !    internal variables    
  !
  
  CHARACTER (LEN=10)               :: code = 'RELAX_TIME'
  character (LEN=100)              :: str
  character (LEN=256)              :: outdir                 ! out directory
  character (LEN=256)              :: fmt
  character (LEN=256)              :: makedirectory
  character (LEN=256)              :: out_file               ! output file
  character (LEN=256), external    :: trimcheck
  integer                          :: ios
  integer                          :: ia                     ! atom's index
  integer                          :: ik                     ! k-pt. index
  integer                          :: id                     ! direction
  LOGICAL                          :: compute_forces         ! compute force mtxel
  LOGICAL                          :: exist
  LOGICAL                          :: fully_local_pot        ! logical condition -> local pot.
  integer                          :: ib1, ib2               ! band index
  integer                          :: igr, itr, itr1, itr2   ! transitions
  integer                          :: it                     ! time index
  integer                          :: ierr                   ! error index
  real(DP), allocatable            :: gammat (:)             ! \gamma(t)
  real(DP), allocatable            :: gammast (:)            ! cum. sum
  real(DP), allocatable            :: gamma_aux (:)          ! aux. array
  real(DP), allocatable            :: correlfunc_oft (:)     ! correl. function
  real(DP), allocatable            :: temp_rel_oft (:)       ! temp. relax
  real(DP)                         :: C                      ! prefactor
  real(DP)                         :: de                     ! energy diff.
  real(DP)                         :: r                      ! param.
  real(DP)                         :: prefct                 ! prefactor
  real(DP)                         :: tanhfct                ! tanh factor
  
  !
  !    input file read directly here
  !    DO NOT USE phq_readin
  !
  
  NAMELIST / inputph / outdir, prefix, temprel_fil,                             &
       scr_type, t_ev, temprel_data_dir, st_at, fin_at, fil_chi, amass,         &
       calc_eff_charge, filrhoc, filrhov, filvscr, compute_forces, ntd, t_max,  &
       fully_local_pot
  
  !
  !
  !
#ifdef __MPI
  call mp_startup ( )
#endif
  call environment_start ( code )
  
  !
  prefix = 'pwscf'
  call get_environment_variable ( 'ESPRESSO_TMPDIR', outdir )
  IF ( TRIM( outdir ) == ' ' ) outdir = './'
  
  IF ( npool > 1 ) call errore ('temp_relax','pools not implemented',npool)
  !
  IF ( ionode ) THEN
     !
     call input_from_file ( )
     !
     READ (5, inputph, err=200, iostat=ios)
200  call errore ('temp_relax','reading inputph namelist', ABS (ios) )
     !
     tmp_dir = trimcheck (outdir)
     !
  END IF
  
  !
  !   ... broadcast variable
  !
  
  call mp_bcast (tmp_dir, ionode_id, intra_image_comm)
  call mp_bcast (prefix, ionode_id, intra_image_comm)
  call mp_bcast (t_ev, ionode_id, intra_image_comm)
  call mp_bcast (scr_type, ionode_id, intra_image_comm)
  call mp_bcast (fil_chi, ionode_id, intra_image_comm)
  call mp_bcast (temprel_data_dir, ionode_id, intra_image_comm)
  call mp_bcast (st_at, ionode_id, intra_image_comm)
  call mp_bcast (fin_at, ionode_id, intra_image_comm)
  call mp_bcast (calc_eff_charge, ionode_id, intra_image_comm)
  call mp_bcast (compute_forces, ionode_id, intra_image_comm)
  call mp_bcast (ntd, ionode_id, intra_image_comm)
  call mp_bcast (t_max, ionode_id, intra_image_comm)
  call mp_bcast (amass, ionode_id, intra_image_comm)
  call mp_bcast (fully_local_pot, ionode_id, intra_image_comm)
  
  !
  !   ... allocate space for pwscf variables
  !
  
  call read_file
  call openfil_pp
  
  call weights ( )
  !
  write(6,*) npwx
  write(6,*) dffts%nnr
  write(6,*) "screening type= ", trim(scr_type)
  write(6,*) "elec. temperature (eV)= ", t_ev
  write(6,*) "compute forces= ", compute_forces
  !
  call init_us_1
  
  !
  !   ... set up output directory
  !
  
  IF ( ionode ) THEN
     !
     makedirectory = "mkdir " // trim(temprel_data_dir)
     !
     inquire (file=TRIM (temprel_data_dir), exist=exist)
     IF ( .not. exist ) call system ( makedirectory )
     !
  END IF
  !
  call mp_barrier (intra_image_comm)
  
  !
  !   ... initialize temp. rel. variables
  !
  
  call initialize_trel ( )
  
  !
  write(6,*) "electron density (cm^-3)= ", n_elec
  write(6,*) "number of bands= ", nbnd
  write(6,*) "number of atoms= ", nat
  write(6,*) "M (eV fs^2/bohr^2)= ", M
  write(6,*) nelec
  
  !
  !   ... set up occupation probabilities
  !
  
  call set_pocc ( )
  
  !
  IF ( .not. compute_forces ) GOTO 100
  
  !
  !     calculation of electron-ion force fx(r - R0) for each ion
  !
  
  IF (.not. fully_local_pot) THEN
     
     !
     !   allocate non local arrays
     !
     
     call allocate_trel_nl ( )
     
     !
     !   initialize all arrays
     !
     
     call init_trel_arrays ( )
     
     !
  END IF
  
  !  ----------------------------------------------------------------------------------------------
  !     matrix elements definition:
  !     
  !     f_nk,mk(R0)^x = 1/\Omega \int_\Omega dr u_{n,k}*(r) fx(r - R0) u_{m,k}(r)
  !
  !   ---------------------------------------------------------------------------------------------
  
  !
  !     start iteration over atoms
  !
  
  DO ia= st_at, fin_at
     
     !
     !    compute external force -> atom ia
     !
     
     IF ( upf (ityp(ia))%tpawp ) THEN
        !
        !   compute potential inside core region
        !
        !
        !   compute potential over QE grid
        !
        call set_paw_dvie_dr ( ia )
        !
     ELSE
        !
        !   compute local potential derivative
        !
        call set_dvie_dr ( ia )
        !
     END IF
     
     !
     !  iterate over k pts
     !
     
     do ik= 1, nks
        
        !
        !    set up set of transitions:
        !    force matrix elements
        !
        
        call set_force_transitions ( ik )
        
        !
        !    set up the matrix elements:  f_nk,mk(R0)^x 
        !
        
        IF ( upf (ityp(ia))%tpawp ) THEN
           !
           !   compute matrix elements inside the core region
           !
           
           !
           !   TO BE IMPLEMENTED
           !
        ELSE
           
           !
           !   perform matrix elements calculation over space grid
           !
           
           IF ( fully_local_pot ) THEN
              call set_fmtxel_smoothg ( ia, ik )
           ELSE
              call set_fmtxel_rec_space ( ia, ik )
           END IF
           !
           
        END IF
        
        !
        !   write matrix elements data on external file
        !
        
        IF ( ionode ) THEN
           
           !
           !   write file name
           !
           
           IF (ia < 10) THEN
              write (str, '(I1.1)') ia
           ELSE IF (ia < 100) THEN
              write (str, '(I2)') ia
           ELSE
              write (str, '(I3)') ia
           END IF
           
           !
           out_file = TRIM(temprel_data_dir) // "/force_mtxel_at" // TRIM( adjustl (str) ) // ".txt"
           !
           inquire (file=TRIM (out_file), exist=exist)
           IF ( exist ) THEN
              open (100, file=TRIM (out_file), status="old", position="append", action="write")
           ELSE
              open (100, file=TRIM (out_file), status="new", action="write")
           END IF
           !
           write (fmt,*) '(I5.1,4X,I5.1,4X,I5.1,4X,I3.1,4X,E14.7,4X,E14.7)'
           !
           
           !
           !   run over trans. groups 
           !
           
           do igr= 1, ngr
              
              !
              !    transition indexes
              !
              
              itr1 = trans_groups (igr,1)
              itr2 = trans_groups (igr,2)
              
              !
              !   band index
              !
              
              ib1 = trans_list (itr1,1)
              !
              
              do itr= itr1, itr2
                 
                 !
                 ib2 = trans_list (itr,2)
                 !
                 do id= 1, 3
                    !
                    
                    !   write matrix elements
                    !
                    write (100,fmt) ib1, ib2, ik, id, real ( psi_gradv_psi (itr,id) ),        &
                         aimag ( psi_gradv_psi (itr,id) )
                    !
                 end do
                 !
              end do
              !
           end do
           !
           close (100, iostat=ios)
           !
        END IF
        !
        call mp_barrier ( intra_image_comm )
        !
        
     end do
     
     !
     !    end of calculation for atom ia
     !
     
     !
     write(6,*)  "atom ", ia, " computed"
     !
  END DO
  !
  GOTO 101
  
100 CONTINUE
  !
  !   ... allocate space for arrays
  !
  
  allocate ( gammat (1:ntd), STAT=ierr )
  if (ierr/=0) call errore ('temp_relax','allocating gammat', ABS(ierr))
  !
  allocate ( gammast (1:ntd), STAT=ierr )
  if (ierr/=0) call errore ('temp_relax','allocating gammast', ABS(ierr))
  !
  allocate ( gamma_aux (1:ntd), STAT=ierr )
  if (ierr/=0) call errore ('temp_relax','allocating gamma_aux', ABS(ierr))
  !
  
  !
  !   ... allocate temp. relax arrays
  !
  
  allocate ( correlfunc_oft (1:ntd), STAT=ierr )
  if (ierr/=0) call errore ('temp_relax','allocating correlfunc_oft', ABS(ierr))
  correlfunc_oft (:) = 0._dp
  !
  allocate ( temp_rel_oft (1:ntd), STAT=ierr )
  if (ierr/=0) call errore ('temp_relax','allocating temp_rel_oft', ABS(ierr))
  temp_rel_oft (:) = 0._dp
  
  !
  !  ----------------------------------------------------------------------------------------------
  !
  !       starting full calculation
  !       g_{ie}(T) = -\pi\hbar/M \sum_{R0}\sum_x \sum_{n,m}'\sum_k Wk (Pn(k)-Pm(k))/(En(k)-Em(k))
  !       f_nk,mk(R0)^x f_mk,nk(R0)^x \delta(En(k)-Em(k))
  !
  !  ----------------------------------------------------------------------------------------------
  !
  
  DO ia= 1, nat
     
     !
     !    acquire matrix elements from file
     !
     
     call read_matrxel_file ( ia )
     
     !
     !    atom's mass
     !
     
     M = 0._dp
     M = amass (ityp (ia)) * proton_mass
     
     !
     !    compute relaxation rate
     !
     
     DO id= 1, 3
        
        !
        !    initialize arrays
        !
        
        gammat (:) = 0._dp
        gammast (:) = 0._dp
        gamma_aux (:) = 0._dp
        
        !
        !    run over transitions
        !
        
        do itr= trans_list_auxp (me_image+1, 1), trans_list_auxp (me_image+1, 2)
           
           !
           !   select band indexes
           !
           
           ib1 = trans_list_full (itr,1)
           ib2 = trans_list_full (itr,2)
           ik = trans_list_full (itr,3)
           
           !
           de = et (ib1,ik) - et (ib2,ik)
           de = de * RYTOEV                ! eV
           
           !
           !   prefactor  
           !
           
           prefct = pocc (ib1,ik) + pocc (ib2,ik) - 2._dp * pocc (ib1,ik) * pocc (ib2,ik)
           
           !
           r = real ( psi_gradv_psi (itr,id) * conjg ( psi_gradv_psi (itr,id) ) +              &
                conjg ( psi_gradv_psi (itr,id) ) * psi_gradv_psi (itr,id) )
           !
           
           IF (ib1 .NE. ib2) THEN
              !
              IF (ABS (de) == 0._dp) THEN
                 !
                 gammat (:) = wk (ik) * prefct * r / pi + gammat (:)
                 !
                 !    cumulative sum
                 !
                 gammast (:) = wk (ik) * prefct * r * t (:) / pi + gammast (:)
                 !
              ELSE
                 !
                 tanhfct = TANH (de / (2._dp * t_ev)) / (de / (2._dp * t_ev))
                 !
                 gammat (:) = wk (ik) * prefct * r * tanhfct * cos (de * t (:)) / pi +           &
                      gammat (:)
                 !   eV^2 / bohr^2
                 !
                 !   cumulative sum
                 !
                 gammast (:) = wk (ik) * prefct * r * tanhfct * sin (de * t (:)) / (pi * de) +   &
                      gammast (:)
                 !
                 !   eV / bohr^2
                 !
              END IF
              !
           END IF
           
           !
        end do         ! end transition cycle
        !
        call mp_barrier ( intra_image_comm )
        
        !
        !     include prefactor
        !
        
        gammat (:) = pi / (2._dp * M * t_ev) * gammat (:)
        !
        !   bohr^2 / (eV^2 fs^2) * eV^2 / bohr^2 = fs^-2
        !
        
        gammast (:) = pi * hbar / (2._dp * M * t_ev) * gammast (:)
        !
        !   eV/bohr^2 * fs/(eV fs^2)*bohr^2=
        !   = fs^-1
        !
        
        !
        !   collect data into single proc.
        !
        
        gamma_aux (:) = gammat (:)
        call mp_sum ( gamma_aux (:), world_comm )
        gammat (:) = gamma_aux (:)
        !
        gamma_aux (:) = 0._dp
        gamma_aux (:) = gammast (:)
        call mp_sum ( gamma_aux (:), world_comm )
        gammast (:) = gamma_aux (:)
        
        !
        call mp_barrier ( intra_image_comm )
        
        !
        !   add to energy relaxation
        !
        
        correlfunc_oft (:) = correlfunc_oft (:) + gammat (:)
        !
        temp_rel_oft (:) = temp_rel_oft (:) + gammast (:)
        !
        !
     END DO
     !
     
     write(6,*) "atom= ", ia, " completed"
     !
  END DO
  
  !
  !   write data on file
  !
  
  IF ( ionode ) THEN
     !
     open (100, file=temprel_fil, status="new", action="write")
     !
     write (fmt,*) '(E14.7,4X,E14.7,4X,E14.7)'
     !
     !    prefactor
     !
     
     C = ni * K_BOLTZMANN_SI / nat
     !
     !    m^-3 * J / K
     correlfunc_oft (:) = correlfunc_oft (:) * C * 1.E15
     !
     !    W/(m^3 K fs)
     
     temp_rel_oft (:) = temp_rel_oft (:) * C * 1.E15
     !
     !    W/(m^3 K)
     
     !
     !    run over time array
     !
     
     do it= 1, ntd
        !
        write (100,fmt) (hbar * t (it)), correlfunc_oft (it), temp_rel_oft (it)
        !
     end do
     !
     close (100, iostat=ios)
     
     !
  END IF
  !
101 CONTINUE
  call mp_barrier ( intra_image_comm )
  call environment_end ( code )
  STOP
  !
END PROGRAM temp_relax
