! --------------------------------------------------------------------
!
!       MODULE
!       calculation of screened potential
! --------------------------------------------------------------------
!       1)  PAW
!       2)  local pseudo potentials
!
! --------------------------------------------------------------------
!       type of screening:
!       1)  bare
!       2)  lindhard
!       3)  Thomas-Fermi
!       4)  full RPA
! --------------------------------------------------------------------
!
!       Author:  Jacopo Simoni
!
! --------------------------------------------------------------------
!
!

MODULE apawpot
  !
  USE kinds,                   ONLY : DP
  
  !
  SAVE
  
  !
  !     global variables
  !
  
  real(DP), allocatable             :: vsloc_ofG (:,:)
  !     screened potential in rec. space
  real(DP), allocatable             :: vsloc_ofG_nc (:,:,:)
  !     screened nc potential
  complex(DP), allocatable          :: vsloc_ofGw (:,:,:)
  !     screened potential of w
  real(DP), allocatable             :: vsqloc_ofG (:,:)
  !     screened potential (q)
  real(DP), allocatable             :: vsqloc_ofG_nc (:,:,:)
  !     screened nc potential (q)
  real(DP), allocatable             :: dr_vie (:,:)
  !     \grad_r Vie_scr(r-Ri) on spatial grid
  real(DP), allocatable             :: ddvie_dr (:,:)
  !     grad_r grad_Ix Vie(r-RI)
  real(DP), allocatable             :: dr_view (:,:,:)
  !     \grad_r Vie_scr(r-Ri,w)
  real(DP), allocatable             :: dr_vss (:,:,:)
  !     \grad_r v_ss(r-Ri)
  logical                           :: finiteq
  !     calc. at q point
  !     -> .TRUE.
  
CONTAINS
  !
  ! ------------------------------------------------------------
  SUBROUTINE set_potentialofG_qegrid ( )
    ! ----------------------------------------------------------
    !
    !    F.T. of coulomb/PAW/pseudo local potential
    !    on Q.E. G grid
    !
    ! ----------------------------------------------------------
    !
    USE dynvars,                 ONLY : dyncalc, nfreq
    USE screening_vars,          ONLY : scr_type
    USE uspp_param,              ONLY : upf
    USE gvect,                   ONLY : ngl, gl, ngm
    USE ions_base,               ONLY : ntyp => nsp
    USE tf_screening,            ONLY : set_tf_screening
    USE lindhard_screening,      ONLY : set_lindhard_screening, set_dyn_screening, set_lindhard_screening_nc
    USE ks_screening,            ONLY : set_ks_screening, set_gq_grid
    USE paw_grids,               ONLY : ngp
    USE potentials,              ONLY : set_vpaw_ofg, vspaw_ofg, vslpaw_ofg, set_vpaw_ofg_qegrid
    USE noncollin_module,        ONLY : noncolin, npol
    USE local_field_corrections, ONLY : set_local_field_factors, set_localq_field_factors
    !
    IMPLICIT NONE
    
    !
    !    internal variables
    !
    integer                     :: ierr
    !    error index
    integer                     :: it
    !    atom type index
    logical                     :: is_any_paw
    !    check if paw potential used
    
    !
    !    ... prepare dielectric function
    !
    
    SELECT CASE ( scr_type )
    CASE ( "tf" )
       !
       !   screening length calculation
       !
       call set_tf_screening ( )
       !
    CASE ( "lindhard" )
       !
       !   set variables for Lindhard calculation
       !
       call set_lindhard_screening ( )
       !
       IF ( dyncalc ) call set_dyn_screening ( )
       !
       IF ( noncolin ) THEN
          !
          !   NC lindhard
          !
          call set_lindhard_screening_nc ( )
          !
          !   local field factors
          !
          IF (finiteq) THEN
             call set_localq_field_factors ( )
          ELSE
             call set_local_field_factors ( )
          END IF
          !
       END IF
       !
    CASE ( "RPA" )
       !
       !   set Gq grid
       !
       call set_gq_grid ( )
       !
       !   static KS screening calculation
       !
       call set_KS_screening ( )
       !
    CASE DEFAULT
       !
       write(6,*) "BARE calculation"
       !
    END SELECT
    
    !
    !    check if paw potential is used
    !
    
    is_any_paw = .FALSE.
    !
    do it= 1, ntyp
       IF ( upf(it)%tpawp ) is_any_paw = .TRUE.
    end do
    
    IF ( scr_type == "RPA" .and. is_any_paw )            &
         call errore ('set_potentialofG_qegrid','PAW + RPA not implemented', 1)
    !
    
    !
    !    allocate final screened potential
    !
    IF ( is_any_paw ) THEN
       !
       IF ( .not. ALLOCATED ( vspaw_ofg ) ) THEN
          ALLOCATE ( vspaw_ofg (1:ngp,1:ntyp), STAT=ierr )
          if (ierr/=0) call errore ('set_potentialofG_qegrid','allocating vspaw_ofg',ABS(ierr))
       END IF
       vspaw_ofg (:,:) = 0._dp
       !
       IF ( .not. ALLOCATED ( vslpaw_ofg ) ) THEN
          ALLOCATE ( vslpaw_ofg (1:ngl,1:ntyp), STAT=ierr )
          if (ierr/=0) call errore ('set_potentialofG_qegrid','allocating vslpaw_ofg',ABS(ierr))
       END IF
       vslpaw_ofg (:,:) = 0._dp
       !
    ELSE
       !      
       IF ( finiteq ) THEN
          !
          IF (noncolin) THEN
             IF (.not. ALLOCATED ( vsqloc_ofG_nc ) ) THEN
                ALLOCATE ( vsqloc_ofG_nc (1:ngm,1:npol,1:ntyp), STAT=ierr )
                if (ierr/=0) call errore ('set_potentialofG_qegrid','allocating vsqloc_ofG_nc',ABS(ierr))
             END IF
             vsqloc_ofG_nc (:,:,:) = 0._dp
             !
          ELSE
             !
             IF (.not. ALLOCATED ( vsqloc_ofG ) ) THEN
                ALLOCATE ( vsqloc_ofG (1:ngm,1:ntyp), STAT=ierr )
                if (ierr/=0) call errore ('set_potentialofG_qegrid','allocating vsqloc_ofG',ABS(ierr))
             END IF
             vsqloc_ofG (:,:) = 0._dp
             !
          END IF
          !
       ELSE
          !
          IF ( dyncalc ) THEN
             !
             IF ( scr_type == "RPA" ) THEN
                !
                ALLOCATE ( vsloc_ofGw (1:ngm,1:nfreq,1:ntyp), STAT=ierr )
                if (ierr/=0) call errore ('set_potentialofG_qegrid','allocating vsloc_ofGw',ABS(ierr))
                vsloc_ofGw (:,:,:) = cmplx (0._dp, 0._dp, kind=DP)
                !
             ELSE
                !
                ALLOCATE ( vsloc_ofGw (1:ngl,1:nfreq,1:ntyp), STAT=ierr )
                if (ierr/=0) call errore ('set_potentialofG_qegrid','allocating vsloc_ofGw',ABS(ierr))
                vsloc_ofGw (:,:,:) = cmplx (0._dp, 0._dp, kind=DP)
                !
             END IF
             !
          ELSE
             !
             IF ( scr_type == "RPA" ) THEN
                !
                ALLOCATE ( vsloc_ofG (1:ngm,1:ntyp), STAT=ierr )
                if (ierr/=0) call errore ('set_potentialofG_qegrid','allocating vsloc_ofG',ABS(ierr))
                vsloc_ofG (:,:) = 0._dp
                !
             ELSE
                !
                ALLOCATE ( vsloc_ofG (1:ngl,1:ntyp), STAT=ierr )
                if (ierr/=0) call errore ('set_potentialofG_qegrid','allocating vsloc_ofG',ABS(ierr))
                vsloc_ofG (:,:) = 0._dp
                !
                IF ( noncolin ) THEN
                   !
                   ALLOCATE ( vsloc_ofG_nc (1:ngl,1:npol,1:ntyp), STAT=ierr )
                   if (ierr/=0) call errore ('set_potentialofG_qegrid','allocating vsloc_ofG_nc',ABS(ierr))
                   vsloc_ofG_nc (:,:,:) = 0._dp
                   !
                END IF
                !
             END IF
             !
          END IF
          !
       ENDIF
       !
    END IF
    
    !
    !    starting cycle iteration
    !
    
    DO it= 1, ntyp
       !
       !    define vsloc_ofG
       !
       IF ( upf(it)%tcoulombp ) THEN
          !
          IF ( finiteq ) THEN
             call vscrq_coul_ofg ( it )
          ELSE
             call vscr_coul_ofg ( it )
          END IF
          !
       ELSE IF ( upf(it)%tpawp ) THEN
          !
          !  v(G) -> PAW more G vectors
          !
          call set_vpaw_ofg ( it )
          !
          !  v(G) -> PAW QE grid
          !
          call set_vpaw_ofg_qegrid ( it )
          !
       ELSE
          !
          IF ( finiteq ) THEN
             IF ( noncolin ) THEN
                call vscrq_loc_ofg_nc ( it )
             ELSE
                call vscrq_loc_ofg ( it )
             END IF
          ELSE
             IF ( dyncalc ) THEN
                call vscr_loc_ofgw ( it, gl, ngl )
             ELSE
                call vscr_loc_ofg ( it, gl, ngl )
                !
                IF ( noncolin ) call vscr_loc_ofg_nc ( it, gl, ngl )
                !
             END IF
          END IF
          !
       END IF
       !
    END DO
    !
    return
    !
  END SUBROUTINE set_potentialofG_qegrid
  !
  ! -------------------------------------------------------------------
  SUBROUTINE set_dvie_dr ( at_ind )
    ! -----------------------------------------------------------------
    !
    !    calculation of periodic potential gradient
    !    on the full spatial grid
    !    \grad_{Ra} v_ie ( r - Ra )
    !
    !    input1:   atom index
    !
    ! -----------------------------------------------------------------
    USE constants,         ONLY : tpi, RYTOEV
    USE cell_base,         ONLY : tpiba
    USE ions_base,         ONLY : tau, ityp
    USE fft_base,          ONLY : dffts, dfftp
    USE control_flags,     ONLY : gamma_only
    USE gvect,             ONLY : ngm, g, igtongl
    USE mp_world,          ONLY : world_comm
    USE Coul_cut_2D,       ONLY : do_cutoff_2D, lr_Vloc
    USE fft_interfaces,    ONLY : invfft
    USE screening_vars,    ONLY : scr_type
    USE noncollin_module,  ONLY : noncolin
    
    !
    IMPLICIT NONE
    !
    !    input variables
    !
    integer, intent(in)               :: at_ind
    !    atomic index to compute
    !
    !    internal variables
    !
    complex(DP), allocatable          :: aux (:)
    !    auxiliary array
    integer                           :: ierr
    !    error index
    integer                           :: id
    !    direction index
    integer                           :: ig
    !    G vector index
    real(DP)                          :: arg
    !    auxiliary variable
    complex(DP)                       :: caux
    !    e^{-i G Ri}
    
    !
    !    dr_vie allocation
    !
    
    IF ( .not. allocated ( dr_vie ) ) THEN
       
       !
       !    allocate space for periodic
       !    potential array
       !
       
       allocate ( dr_vie ( 1:dffts%nnr, 1:3 ), STAT=ierr )
       if (ierr/=0) call errore('apaw_pot:set_dvie_dr','allocating dr_vie',ABS(ierr))
       !
       
    END IF
    
    !
    dr_vie (:,:) = 0._dp
    
    !
    !    allocate auxiliary fields
    !
    
    allocate ( aux ( 1:dffts%nnr ), STAT=ierr )
    if (ierr/=0) call errore('apaw_pot:set_dvie_dr','allocating aux',ABS(ierr))
    
    ! --------------------------------------------------------------------------
    !
    !    potential in reciprocal space grid
    !    already computed
    !
    ! 
    !    compute -Grad_i V_scr (G) for given atom
    !    V(r) = \Sum_i \Sum_n v_Ie(r - Ri + nL)
    !
    ! --------------------------------------------------------------------------
    
    do id= 1, 3                    ! iteration over directions
       !
       aux (:) = cmplx (0._dp, 0._dp, kind=DP)
       !
       do ig= 1, ngm
          !
          !    arg= i Gx e^{-i G Ri} vloc(G)
          !
          arg = ( g (1,ig) * tau (1,at_ind) + g (2,ig) * tau (2,at_ind) +      &
               g (3,ig) * tau (3,at_ind) ) * tpi
          !
          caux = cmplx ( sin(arg), cos(arg), kind=DP )
          !
          !    grad_r v(r-Ri) (Ry/bohr)
          !
          IF ( scr_type == "RPA" ) THEN
             !
             aux ( dfftp%nl(ig) ) = aux ( dfftp%nl(ig) ) + caux * g (id,ig) *                &
                  vsloc_ofG (ig, ityp (at_ind)) * tpiba
             !
          ELSE
             !
             aux ( dfftp%nl(ig) ) = aux ( dfftp%nl(ig) ) + caux * g (id,ig) *                &
                  vsloc_ofG (igtongl (ig), ityp (at_ind)) * tpiba
             !
          END IF
          !
       end do
       !
       IF ( gamma_only ) THEN
          !
          do ig= 1, ngm
             aux ( dfftp%nlm (ig) ) = CONJG ( aux ( dfftp%nl (ig) ) )
          end do
          !
       END IF
       !
       IF ( do_cutoff_2D ) THEN
          !
          !  ... re-add the cutoff FT of erf
          !
          do ig= 1, ngm
             !
             aux ( dfftp%nl(ig) ) = aux ( dfftp%nl(ig) ) + caux * g (id,ig) *              &
                  lr_Vloc (ig,ityp (at_ind)) * tpiba
             !
          end do
          !
       END IF
       !
       !  ...  aux = potential gradient in G space -> FFT to real space
       !
       call invfft ( 'Rho', aux, dfftp )
       !
       dr_vie (:,id) = DBLE ( aux (:) )
       !
    end do
    
    !
    !    potential derivative (eV/bohr) units
    !
    dr_vie (:,:) = dr_vie (:,:) * RYTOEV
    !
    !    for non collinear calculation compute
    !    dvss_dr
    !
    IF (noncolin) call set_dvss_dr ( at_ind )
    !
    return
    !
  END SUBROUTINE set_dvie_dr
  !
  ! ----------------------------------------------------------------------
  SUBROUTINE set_dvss_dr ( at_ind )
    ! -----------------------------------------------------------------
    !
    !    calculation of periodic potential gradient
    !    on the full spatial grid
    !    \grad_{Ra} v_ie ( r - Ra )
    !
    !    input1:   atom index
    !
    ! -----------------------------------------------------------------
    USE constants,         ONLY : tpi, RYTOEV
    USE cell_base,         ONLY : tpiba
    USE ions_base,         ONLY : tau, ityp
    USE fft_base,          ONLY : dffts, dfftp
    USE control_flags,     ONLY : gamma_only
    USE gvect,             ONLY : ngm, g, igtongl
    USE mp_world,          ONLY : world_comm
    USE Coul_cut_2D,       ONLY : do_cutoff_2D, lr_Vloc
    USE fft_interfaces,    ONLY : invfft
    USE screening_vars,    ONLY : scr_type
    USE noncollin_module,  ONLY : npol
    
    !
    IMPLICIT NONE
    !
    !    input variables
    !
    integer, intent(in)               :: at_ind
    !    atomic index to compute
    !
    !    internal variables
    !
    complex(DP), allocatable          :: aux (:)
    !    auxiliary array
    integer                           :: ierr
    !    error index
    integer                           :: id
    !    direction index
    integer                           :: isp
    !    spin index
    integer                           :: ig
    !    G vector index
    real(DP)                          :: arg
    !    auxiliary variable
    complex(DP)                       :: caux
    !    e^{-i G Ri}
    
    !
    !    dr_vie allocation
    !
    
    IF ( .not. allocated ( dr_vss ) ) THEN
       
       !
       !    allocate space for periodic
       !    potential array
       !
       
       allocate ( dr_vss ( 1:dffts%nnr, 1:npol, 1:3 ), STAT=ierr )
       if (ierr/=0) call errore('apaw_pot:set_dvss_dr','allocating dr_vss',ABS(ierr))
       !
       
    END IF
    
    !
    dr_vss (:,:,:) = 0._dp
    
    !
    !    allocate auxiliary fields
    !
    
    allocate ( aux ( 1:dffts%nnr ), STAT=ierr )
    if (ierr/=0) call errore('apaw_pot:set_dvie_dr','allocating aux',ABS(ierr))
    
    ! --------------------------------------------------------------------------
    !
    !    potential in reciprocal space grid
    !    already computed
    !
    ! 
    !    compute -Grad_i V_scr (G) for given atom
    !    V(r) = \Sum_i \Sum_n v_Ie(r - Ri + nL)
    !
    ! --------------------------------------------------------------------------
    
    do id= 1, 3                    ! iteration over directions
       !
       do isp= 1, npol             ! iteration over spin
          !
          aux (:) = cmplx (0._dp, 0._dp, kind=DP)
          !
          do ig= 1, ngm
             !
             !    arg= i Gx e^{-i G Ri} vloc(G)
             !
             arg = ( g (1,ig) * tau (1,at_ind) + g (2,ig) * tau (2,at_ind) +      &
                  g (3,ig) * tau (3,at_ind) ) * tpi
             !
             caux = cmplx ( sin(arg), cos(arg), kind=DP )
             !
             !    grad_r vss(r-Ri) (Ry/bohr)
             !
             IF ( scr_type == "lindhard" ) THEN
                !
                aux ( dfftp%nl(ig) ) = aux ( dfftp%nl(ig) ) + caux * g (id,ig) *                &
                     vsloc_ofG_nc (igtongl (ig), isp, ityp (at_ind)) * tpiba
                !
             ELSE
                !
                call errore('apaw_pot:set_dvss_dr','only lindhard screening implemented', 1)
                !
             END IF
             !
          end do
          !
          IF ( gamma_only ) THEN
             !
             do ig= 1, ngm
                aux ( dfftp%nlm (ig) ) = CONJG ( aux ( dfftp%nl (ig) ) )
             end do
             !
          END IF
          !
          IF ( do_cutoff_2D ) THEN
             !
             !  ... re-add the cutoff FT of erf
             !
             do ig= 1, ngm
                !
                aux ( dfftp%nl(ig) ) = aux ( dfftp%nl(ig) ) + caux * g (id,ig) *              &
                     lr_Vloc (ig,ityp (at_ind)) * tpiba
                !
             end do
             !
          END IF
          !
          !  ...  aux = potential gradient in G space -> FFT to real space
          !
          call invfft ( 'Rho', aux, dfftp )
          !
          dr_vss (:,isp,id) = DBLE ( aux (:) )
          !
       end do
       !
    end do
    !
    !    potential derivative (eV/bohr) units
    !
    dr_vss (:,:,:) = dr_vss (:,:,:) * RYTOEV
    !
    return
    !
  END SUBROUTINE set_dvss_dr
  !
  ! ----------------------------------------------------------------------
  SUBROUTINE set_dview_dr ( at_ind )
    ! --------------------------------------------------------------------
    !
    !    calculation of periodic potential gradient (freq. dependent)
    !    on the full spatial grid
    !    \grad_{Ra} v_ie ( r - Ra, w )
    !
    !    input1:   atom index
    !
    ! --------------------------------------------------------------------
    
    USE constants,        ONLY : RYTOEV, tpi
    USE fft_base,         ONLY : dffts, dfftp
    USE fft_interfaces,   ONLY : invfft
    USE Coul_cut_2D,      ONLY : do_cutoff_2D, lr_Vloc
    USE control_flags,    ONLY : gamma_only
    USE cell_base,        ONLY : tpiba
    USE ions_base,        ONLY : tau, ityp
    USE gvect,            ONLY : ngm, g, igtongl
    USE screening_vars,   ONLY : scr_type
    USE dynvars,          ONLY : nfreq
    
    !
    implicit none
    !
    !    input variables
    !
    integer, intent(in)         :: at_ind
    !    atom index
    
    !
    !    internal variables
    !
    
    complex(DP), allocatable    :: aux (:)
    !    aux array
    real(DP)                    :: arg
    complex(DP)                 :: caux
    !    aux. variable
    integer                     :: id
    !    dir. index
    integer                     :: ierr
    integer                     :: ig
    integer                     :: iw
    !    freq. index
    
    !
    IF ( .not. allocated ( dr_view ) ) THEN
       !
       !    allocate space for
       !    periodic array potential
       !
       allocate ( dr_view (1:dffts%nnr, 1:nfreq, 1:3), STAT=ierr )
       if (ierr/=0) call errore('set_dview_dr','allocating dr_view',ABS(ierr))
       !
    END IF
    
    !
    dr_view (:,:,:) = 0._dp
    
    !
    !    allocate auxiliary fields
    !
    
    allocate ( aux (1:dffts%nnr), STAT=ierr )
    if (ierr/=0) call errore ('set_dview_dr','allocating aux',ABS(ierr))
    
    ! --------------------------------------------------------------------------
    !
    !    potential in reciprocal space grid
    !    already computed
    !
    ! 
    !    compute -Grad_i V_scr (G) for given atom
    !    V(r) = \Sum_i \Sum_n v_Ie(r - Ri + nL)
    !
    ! --------------------------------------------------------------------------
    
    do id= 1, 3               ! iterate over directions
       !
       do iw= 1, nfreq        ! iterate over frequencies
          !
          aux (:) = cmplx (0.d0, 0.d0, kind=DP)
          !
          do ig= 1, ngm
             !
             !    arg= i Gx e^{-i G Ri} vloc(G)
             !
             arg = ( g (1,ig) * tau (1,at_ind) + g (2,ig) * tau (2,at_ind) +          &
                  g (3,ig) * tau (3,at_ind) ) * tpi
             !
             caux = cmplx ( sin(arg), cos(arg), kind=DP )
             !
             !    grad_r v(r-Ri, w)   (Ry/bohr)
             !
             IF ( scr_type == "RPA" ) THEN
                !
                aux ( dfftp%nl (ig) ) = aux ( dfftp%nl (ig) ) + caux * g (id,ig) *    &
                     vsloc_ofGw (ig,iw,ityp (at_ind)) * tpiba
                !
             ELSE
                !
                aux ( dfftp%nl (ig) ) = aux ( dfftp%nl (ig) ) + caux * g (id,ig) *    &
                     vsloc_ofGw (igtongl (ig),iw,ityp (at_ind)) * tpiba
                !
             END IF
             !
          end do
          !
          IF ( gamma_only ) THEN
             !
             do ig= 1, ngm
                aux ( dfftp%nlm (ig) ) = CONJG ( aux ( dfftp%nl (ig) ) )
             end do
             !
          END IF
          !
          IF ( do_cutoff_2D ) THEN
             !
             !  ... re-add the cutoff FT of erf
             !
             do ig= 1, ngm
                !
                aux ( dfftp%nl (ig) ) = aux ( dfftp%nl (ig) ) + caux * g (id,ig) *    &
                     lr_Vloc (ig,ityp (at_ind)) * tpiba
                !
             end do
             !
          END IF
          !
          !  ...  aux = potential gradient in G space -> FFT to real space
          !
          call invfft ( 'Rho', aux, dfftp )
          !
          dr_view (:,iw,id) = DBLE ( aux (:) )
          !
       end do
       !
    end do
    !
    !    potential derivative (eV/bohr) units
    !
    dr_view (:,:,:) = dr_view (:,:,:) * RYTOEV
    !
    return
    !
  END SUBROUTINE set_dview_dr
  !
  !
  ! ------------------------------------------------------------
  SUBROUTINE set_ddvie_dr ( idx )
    ! ----------------------------------------------------------
    !
    !      this subroutine computes the force
    !      \nabla_Ix \nabla_r Vscr^ei(r - RI)
    !
    !      to be used in the calculation of f_Ix^SO(r)
    !
    ! ----------------------------------------------------------
    
    USE control_flags,         ONLY : gamma_only
    USE cell_base,             ONLY : tpiba
    USE fft_base,              ONLY : dffts, dfftp
    USE fft_interfaces,        ONLY : fwfft, invfft
    USE gvect,                 ONLY : ngm, g
    
    !
    implicit none
    
    !
    !      internal variables
    !
    
    integer, intent(in)                :: idx
    !       x dir. input
    integer                            :: id, ig
    !       dir.
    integer                            :: ierr
    real(DP), allocatable              :: f_r (:)
    !      aux. field
    complex(DP), allocatable           :: f_G (:)
    !      F.T.
    complex(DP), allocatable           :: aux (:)
    !      aux. field
    
    !
    !      1) dvie_dr around atom I
    !         already computed : main program
    !         extract idx component
    !
    
    !
    !      allocate spatial field
    !
    
    allocate ( f_r (1:dffts%nnr), STAT=ierr )
    if (ierr/=0) call errore ('set_ddvie_dr','allocating f_r', ABS(ierr))
    f_r (:) = 0._dp
    
    !
    f_r (:) = dr_vie (:,idx)      ! eV / bohr
    !
    
    !
    !      allocate field
    !
    
    IF ( .not. allocated ( ddvie_dr ) ) THEN
       !
       allocate ( ddvie_dr (1:dffts%nnr, 1:3), STAT=ierr )
       if (ierr/=0) call errore ('set_ddvie_dr','allocating ddvie_dr', ABS(ierr))
       !
    END IF
    !
    ddvie_dr (:,:) = 0._dp
    !
    
    !
    !     allocate auxiliary field
    !
    
    allocate ( f_G (1:dffts%nnr), STAT=ierr )
    if (ierr/=0) call errore ('set_ddvie_dr','allocating f_G',ABS(ierr))
    f_G (:) = cmplx (0.d0, 0.d0, kind=DP)
    !
    allocate ( aux (1:dffts%nnr), STAT=ierr )
    if (ierr/=0) call errore ('set_ddvie_dr','allocating aux',ABS(ierr))
    
    !
    !     perform f_r FT
    !
    
    f_G (:) = cmplx ( f_r (:), 0.d0, kind=DP )
    !
    call fwfft ('Rho', f_G, dffts)
    
    !
    !     compute potential gradient
    !     in reciprocal space
    !
    
    do id= 1, 3         ! iterate over dir
       !
       aux (:) = cmplx (0.d0, 0.d0, kind=DP)
       aux (:) = f_G (:)
       
       !
       do ig= 1, ngm
          !
          aux (dfftp%nl (ig)) = aux (dfftp%nl (ig)) +                 &
               cmplx (0.d0, g (id,ig), kind=DP) * f_G (dfftp%nl (ig)) * tpiba
          !
       end do
       
       !
       IF ( gamma_only ) THEN
          !
          do ig= 1, ngm
             aux ( dfftp%nlm (ig) ) = CONJG ( aux ( dfftp%nl (ig) ) )
          end do
          !
       END IF
       !
       
       !
       !     FFT back in real space
       !
       
       call invfft ( 'Rho', aux, dfftp )
       
       !
       ddvie_dr (:,id) = DBLE ( aux (:) )      ! eV / bohr^2
       !
    end do
    
    !
    return
    !
  END SUBROUTINE set_ddvie_dr
  !
  !  
  ! ------------------------------------------------------------
  SUBROUTINE set_paw_dvie_dr ( at_ind )
    ! ----------------------------------------------------------
    !
    !        This subroutine set up the screened dvie_dr
    !        1) inside the atomic core region: |r-Ra| < rc^a
    !           ... for atom a
    !        2) outside the core region
    !
    !        in region 1: dvie_dr = dvhxc_dr + dvae_dr
    !
    ! ----------------------------------------------------------
    
    USE constants,              ONLY : RYTOEV
    USE paw_variables,          ONLY : paw_info
    USE uspp_param,             ONLY : upf
    USE ions_base,              ONLY : ityp, tau
    USE potentials,             ONLY : clean_potentials, set_wfgrid_dviedr, set_radgrid_dvcdr, potential_wfr
    USE potentials,             ONLY : vslpaw_ofg
    USE fft_base,               ONLY : dffts, dfftp
    USE gvect,                  ONLY : ngm, g, igtongl
    USE cell_base,              ONLY : tpiba
    USE control_flags,          ONLY : gamma_only
    USE Coul_cut_2D,            ONLY : do_cutoff_2D, lr_Vloc
    !
    
    implicit none
    
    !
    !    internal variables 
    !
    
    integer, intent(in)             :: at_ind
    !    atom's index
    TYPE( paw_info )                :: i
    !    minimal atom's info
    integer                         :: id
    !    counter on dim.
    integer                         :: ig
    !    G vector count
    integer                         :: ierr
    complex(DP), allocatable        :: aux (:)
    !    aux. array
    real(DP)                        :: arg
    complex(DP)                     :: caux
    
    !
    !    this subroutine must be called only
    !    if PAW pseudo potential is used 
    !
        
    !
    !    atom's info
    !
    
    i%a = at_ind
    i%t = ityp(i%a)
    
    !
    !    check paw atom
    !
    
    IF ( .not. upf(i%t)%tpawp ) call errore ('set_paw_dvie_dr','paw potential needed',1)
    
    !
    !    deallocate old potentials
    !
    
    call clean_potentials ( )
    
    !
    !    set up the screened potential in reciprocal space
    !
    
    call set_wfgrid_dviedr (i%a)
    
    !
    !    compute all core contributions 
    !
    
    call set_radgrid_dvcdr (i%a)
    stop
    !
    !    compute potential over QE real grid
    !
    
    !
    !    dr_vie allocation
    !
    
    IF ( .not. allocated ( dr_vie ) ) THEN
       !
       !   allocate space for periodic
       !   potential array
       !
       allocate ( dr_vie (1:dffts%nnr, 1:3), STAT=ierr )
       if (ierr/=0) call errore ('apaw_pot:set_paw_dvie_dr','allocating dr_vie',ABS(ierr))
       !
    END IF
    
    !
    dr_vie (:,:) = 0._dp
    
    !
    !    allocate auxiliary fields
    !
    
    allocate ( aux (1:dffts%nnr), STAT=ierr )
    if (ierr/=0) call errore ('apaw_pot:set_paw_dvie_dr','allocating aux',ABS(ierr))
    
    ! -------------------------------------------------
    !
    !    potential in reciprocal space
    !    was already computed
    !
    !    evaluate -Grad_i vscr (G) for atom a
    !    V(r) = \Sum_i \Sum_n v_Ie(r - Ri + nL)
    !
    ! -------------------------------------------------
    
    do id= 1, 3
       !
       aux (:) = cmplx (0._dp, 0._dp, kind=DP)
       !
       do ig= 1, ngm
          !
          !   arg= i Gx e^{-i G Ri} vpaw(G)
          !
          arg = ( g (1,ig) * tau (1,i%a) + g (2,ig) * tau (2,i%a) +          &
               g (3,ig) * tau (3,i%a) )
          !
          caux = cmplx ( sin(arg), cos(arg), kind=DP )
          !
          !   grad_r Vpaw(r - Ri)     (Ry/bohr)
          !
          aux ( dfftp%nl (ig) ) = aux ( dfftp%nl (ig) ) +                    &
               caux * g (id,ig) * vslpaw_ofg ( igtongl (ig), i%t ) * tpiba
          !
       end do
       !
       IF ( gamma_only ) THEN
          !
          do ig= 1, ngm
             aux ( dfftp%nlm (ig) ) = CONJG ( aux ( dfftp%nl (ig) ) )
          end do
          !
       END IF
       !
       IF ( do_cutoff_2D ) THEN
          !
          !   ... re-add the cutoff FT of erf
          !
          do ig= 1, ngm
             !
             aux ( dfftp%nl (ig) ) = aux ( dfftp%nl (ig) ) +                 &
                  caux * g (id,ig) * lr_Vloc (ig,i%t) * tpiba
             !
          end do
          !
       END IF
       !
       !  ... aux = potential gradient in G space -> FFT to real space
       !
       call invfft ( 'Rho', aux, dfftp )
       !
       dr_vie (:,id) = DBLE ( aux (:) )
       !
    end do
    !
    !   eV / bohr units
    !
    dr_vie (:,:) = dr_vie (:,:) * RYTOEV
    !
    return
    !
  END SUBROUTINE set_paw_dvie_dr
  !
  !
  ! 
END MODULE apawpot
