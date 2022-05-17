!
!   author:  Jacopo Simoni
!   this file includes modules needed for the calculation of the
!   screened electron-ion potential
!
!   module : system_vars -> system's variables
!   module : atomic_charge -> set the atom's charges
!   module : tf_screening -> set the variables for tf screening
!   module : ld_screening -> set the variables for lindhard screening
!   module : rpa_screening -> set variables for RPA calculation
!   module : rad_mesh -> set the radial mesh variables
!
! -----------------------------------------------------------------------------
!
!
MODULE tf_screening
  ! ------------------------------------------------------------------------------
  !  
  !          in this module we define all the variables
  !          and subroutines needed to compute the tf screening
  !
  ! ------------------------------------------------------------------------------
  USE kinds,           ONLY : DP
  !
  SAVE
  !
  !          TF screening variables         
  !
  real(DP)                  :: ktf          ! screening length (bohr^-1)
  !
CONTAINS
  !
  ! -------------------------------------------------------------------------------
  SUBROUTINE set_tf_screening ( )
    ! -----------------------------------------------------------------------------
    !
    !      this subroutine initialize the tf screening length
    !      ktf
    !
    ! -----------------------------------------------------------------------------
    USE consts,          ONLY : AngtoBohr, e_4pieps0
    USE constants,       ONLY : BOHR_RADIUS_SI, fpi, pi
    USE screening_vars,  ONLY : n_elec, relerr, t_ev
    !
    implicit none
    !
    !      variables
    !
    real(DP)            :: at
    real(DP)            :: eta
    real(DP)            :: betamu     ! mu/kBT
    real(DP)            :: rs         ! screening parameter
    real(DP)            :: rs_cm      ! screening radius (cm)
    real(DP)            :: rs_m       ! screening radius (m)
    real(DP)            :: De         ! screening length
    real(DP)            :: f_p12      ! f{+1/2}(betamu)
    real(DP)            :: f_m12      ! f{-1/2}(betamu)
    real(DP)            :: gamma      ! coupling parameter
    integer             :: ierr       ! error index
    real*8, external    :: fermi
    !
    !      compute ktf : TF screening length
    !
    !      rs in cm
    rs_cm = (3._dp / (fpi * n_elec)) ** (1._dp / 3)
    !      conversion in m
    rs_m = rs_cm / 100._dp
    !      rs parameter
    rs = rs_m / BOHR_RADIUS_SI
    !      coupling parameter
    gamma = e_4pieps0 / (rs_m * t_ev)
    !
    at = (2._dp * rs / gamma) ** (3._dp/2) / (2._dp * pi ** 2)
    !
    eta = 3._dp / fpi * 1._dp / at
    !      chemical potential divided by kBT
    !      beta_mu = mu / kBT
    betamu = fermi (0.5D0, eta, -1)
    call fermid (0.5D0, betamu, relerr, f_p12, ierr)    ! f_{+1/2}(betamu)
    !
    call fermid (-0.5D0, betamu, relerr, f_m12, ierr)   ! f_{-1/2}(betamu)
    !      Thomas-Fermi screening length in m
    De = 1._dp / SQRT (fpi * 1.e6 * n_elec * e_4pieps0 / t_ev * f_m12 / f_p12)
    !
    De = De * 1.e10 * AngtoBohr                         ! scr. length in bohr
    !
    ktf = 1._dp / De                                    ! bohr^-1
    !
    return
    !
  END SUBROUTINE set_tf_screening
  !
END MODULE tf_screening
!
!
MODULE lindhard_screening
  ! ------------------------------------------------------------------------------
  !  
  !          in this module we define all the variables
  !          and subroutines needed to compute the lindhard (FEG) screening
  !
  ! ------------------------------------------------------------------------------
  !
  USE kinds,          ONLY : DP
  !
  SAVE
  !
  !          variables
  !
  complex(DP), parameter   :: zim = cmplx( 0._dp, 0._dp, kind=DP )
  real(DP)                 :: kf            ! fermi moment (bohr^-1)
  real(DP)                 :: N0            ! feg DOS(ef)
  real(DP)                 :: theta         ! kBT / ef
  real(DP)                 :: mu_ef         ! mu(T) / ef
  real(DP), allocatable    :: w (:)         ! freq. (eV)
  real(DP), allocatable    :: w_2ef (:)     ! w / (2*ef)
  complex(DP), allocatable :: zlist (:)     ! cmplx. freq. list: z / (2*ef)
  real(DP)                 :: ef_feg        ! FEG fermi level (eV)
  !
  !   non collinear quantities
  !
  real(DP)                 :: ef_sp (2)     ! fermi energy spin pol.
  real(DP)                 :: kf_sp (2)     ! fermi moment - spin pol.
  real(DP)                 :: theta_sp (2)  ! kT / ef
  real(DP)                 :: muef_sp (2)   ! mu(T) / ef
  real(DP)                 :: N0_sp (2)     ! DOS
  !
CONTAINS
  !
  ! ------------------------------------------------------------------------------
  SUBROUTINE set_lindhard_screening ( )
    ! ----------------------------------------------------------------------------
    !
    !         This subroutine computes the Lindhard (FEG)
    !         model of screening
    !
    ! ----------------------------------------------------------------------------
    USE screening_vars,   ONLY : t_ev
    !
    implicit none
    !
    !         variables
    !
    real(DP)                  :: ef      ! fermi energy
    real(DP)                  :: mu_T    ! f.e.g. chemical potential
    !
    !         FEG kf (bohr^-1)
    !
    call feg_fermi_moment ( )
    !         FEG Fermi energy
    ef = feg_fermi_energy ( )            ! eV
    !
    theta = t_ev / ef                    ! theta = kB T / ef
    !         FEG chemical pot.
    mu_T = feg_chemical_pot ( ef )       ! eV
    !
    mu_ef = mu_T / ef
    !         FEG DOS (fermi level)
    call feg_fermi_dos ( )
    !
    return
    !
  END SUBROUTINE set_lindhard_screening
  !
  ! ---------------------------------------------------------------------------
  SUBROUTINE set_dyn_screening ( )
    ! -------------------------------------------------------------------------
    !
    !        here we set the variables needed for the dynamical
    !        Lindhard dielectric function \epsilon(w)
    !
    ! -------------------------------------------------------------------------
    
    USE constants,      ONLY : RYTOEV
    USE dynvars,        ONLY : dyncalc, en_max, freq_broadening, nfreq
    USE ener,           ONLY : ef
    !
    implicit none
    !
    !   internal variables
    !
    integer                  :: ierr        ! error index
    integer                  :: iw          ! freq. index
    complex(DP)              :: wc
    real(DP)                 :: ddw         ! freq. interval
    
    !
    IF (dyncalc .eqv. .false.)                 &
         call errore ('set_dyn_screening','dyncalc must be .TRUE.',1)
    
    !
    !   define w list frequencies
    !
    IF ( .NOT. ALLOCATED ( w ) ) THEN
       allocate ( w (1:nfreq), STAT=ierr )
       if (ierr/=0) call errore ('set_dyn_screening','w not allocatable', ABS(ierr))
    END IF
    w (:) = 0._dp
    !
    IF ( .NOT. ALLOCATED ( w_2ef ) ) THEN
       allocate ( w_2ef (1:nfreq), STAT=ierr )
       if (ierr/=0) call errore ('set_dyn_screening','w_2ef not allocatable', ABS(ierr))
    END IF
    w_2ef (:) = 0._dp
    !
    !
    !   compute ef_feg
    !
    ef_feg = feg_fermi_energy ( )                      ! eV
    !
    ddw = en_max / (2._dp * ef * RYTOEV) / nfreq       ! no units
    !
    do iw= 1, nfreq
       w_2ef (iw) = (iw - 1) * ddw
    end do
    
    !
    !   full freq. list w (eV)
    !
    w (:) = w_2ef (:) * 2._dp * ef_feg                 ! eV
    !
    !   complex freq. list
    !
    IF ( .NOT. ALLOCATED ( zlist ) ) THEN
       allocate ( zlist (1:nfreq), STAT=ierr )
       if (ierr/=0) call errore ('set_dyn_screening','zlist not allocatable', ABS(ierr))
    END IF
    zlist (:) = cmplx (0.d0, 0.d0, kind=DP)
    !
    do iw= 1, nfreq
       wc = cmplx (w (iw), freq_broadening, kind=DP)
       wc = wc / (2._dp * ef_feg) 
       zlist (iw) = wc                   ! no units
    end do
    !
    return
    !
  END SUBROUTINE set_dyn_screening
  !
  ! ---------------------------------------------------------------------------
  SUBROUTINE set_lindhard_screening_nc ( )
    ! -------------------------------------------------------------------------
    USE noncollin_module,      ONLY : npol
    USE screening_vars,        ONLY : t_ev
    !
    implicit none
    !
    !    internal variables
    !
    integer                         :: is                !    spin ind.
    real(DP)                        :: muT_sp (2)        !    chemical pot.
    !
    !         FEG kf (bohr^-1)
    !
    call feg_fermi_moment_spol ( )
    !
    !         Fermi energy (eV)
    !
    call set_spol_feg_fermi_energy ( )
    !
    theta_sp (:) = t_ev / ef_sp (:)
    muT_sp (:) = 0._dp
    muef_sp (:) = 0._dp
    !
    do is= 1, npol
       !         FEG chemical potential
       muT_sp (is) = feg_chemical_pot ( ef_sp (is) )
       muef_sp (is) = muT_sp (is) / ef_sp (is)
       !
    end do
    !
    !         Fermi DOS
    !
    call feg_fermi_dos_spol ( )
    !
    return
    !
  END SUBROUTINE set_lindhard_screening_nc 
  !
  ! ---------------------------------------------------------------------------
  real(DP) function feg_chemical_pot ( ef )
    ! -------------------------------------------------------------------------
    ! 
    !            computation of F.E.G. chemical potential
    !            \mu(T) = ef * (1 - 1/3 * (\pi/2 * kT / ef)^2)
    !
    ! -------------------------------------------------------------------------
    USE constants,            ONLY : pi
    USE screening_vars,       ONLY : t_ev
    !
    implicit none
    !
    !     internal variables
    !
    real(DP), intent(in)           :: ef     ! fermi energy
    real(DP)                       :: mu_T   ! chemical potential
    !
    mu_T = ef * ( 1._dp - 1._dp/3 * (0.5D0 * pi * t_ev / ef) ** 2 )
    feg_chemical_pot = mu_T                  ! eV
    !
    return
    !
  end function feg_chemical_pot
  !
  ! ----------------------------------------------------------------------------
  SUBROUTINE feg_fermi_dos ( )
    ! --------------------------------------------------------------------------
    !
    !            computation F.E.G. DOS at the Fermi level
    !            N0 = d n / (2 ef)      (d=3)
    !
    ! --------------------------------------------------------------------------
    USE constants,            ONLY : pi
    USE consts,               ONLY : AngtoBohr, elec_mass, hbar
    !
    N0 = elec_mass * kf * AngtoBohr / (pi ** 2 * hbar ** 2)
    !
    return
    !
  END SUBROUTINE feg_fermi_dos
  !
  ! ----------------------------------------------------------------------------
  SUBROUTINE feg_fermi_dos_spol ( )
    ! --------------------------------------------------------------------------
    USE constants,            ONLY : pi
    USE consts,               ONLY : AngtoBohr, elec_mass, hbar
    !
    implicit none
    !
    N0_sp (:) = 0._dp
    N0_sp (:) = elec_mass * kf_sp (:) * AngtoBohr / (pi ** 2 * hbar ** 2)
    !
    return
    !
  END SUBROUTINE feg_fermi_dos_spol
  !
  ! ----------------------------------------------------------------------------
  complex(DP) function chiT(q,z,theta,mu_eF,dim)
    ! --------------------------------------------------------------------------
    ! Evaluate the finite-temperature Lindhard function
    ! at wave vector "q" and complex frequency "z"
    ! devided by the density of states "N0"
    !
    ! "q" and "z" are in reduced units such that
    !             "q==q/kF" where kF is the Fermi momentum
    !             "z==hbar z/2eF" where eF=hbar**2kF**2/2m is the Fermi energy
    !
    ! "theta"=(kB*T)/eF
    ! "mu_eF"=mu(T)/eF where mu(T) is the chemichal potential at temperature T
    ! --------------------------------------------------------------------------
    implicit none
    !
    ! argument variables
    !
    real(DP)    :: q
    complex(DP) :: z
    real(DP)    :: theta,mu_eF
    integer     :: dim
    !
    ! local variables
    !
    complex(DP) :: zloc,temp
    real(DP)    :: qloc,y,cosh2,arg
    real(DP)    :: dy
    real(DP)    :: fsumrule
    !
    dy=0.0001
    zloc=0.
    !
    if (theta.gt.0) then  ! finite temperature
       !
       chiT=0.
       fsumrule=0.
       y=0.
       !
       do while (abs(fsumrule-1.).gt.4.E-3)
          !  
          y=y+dy
          !
          arg=0.5*(y-mu_eF)/theta
          !
          cosh2=dcosh(arg)
          cosh2=cosh2**2
          !
          qloc=q/sqrt(y)
          zloc=z/y
          !
          temp=chi0(qloc,zloc,dim)
          !
          temp=temp/cosh2
          !
          fsumrule=fsumrule+dy*0.25*y*sqrt(y)/(theta*cosh2)
          !
          chiT=chiT+dy*temp*sqrt(y)
          !
          !        write(6,*) 'arg',arg,fsumrule-1.,y
          !
       end do
       !
       chiT=0.25*chiT/theta
       !
    else ! zero temperature
       !
       chiT=chi0(q,z,dim)
       !
    end if
    !
    return
    !
  end function chiT
  !
  ! ---------------------------------------------------------------------------
  complex(DP) function chi0(q,z,dim)
    ! -------------------------------------------------------------------------
    ! Evaluate the zero-temperature Lindhard function
    ! at wave vector "q" and complex frequency "z"
    ! devided by the density of states "N0"
    !
    ! "q" and "z" are in reduced units such that
    !             "q==q/kF" where kF is th eFermi momentum
    !             "z==hbar z/2eF" where eF=hbar**2kF**2/2m is the Fermi energy
    !
    ! "dim" is the system's dimension
    !
    ! REFERENCE :: "Quantum Theory of the Electron Liquid"
    !              G.F. Giuliani and G. Vignale
    !              Eq. 4.21 p. 161 and Table 4.1 p. 162 for the response function
    !              Table 4.2 p. 163 for the response function at "z"=0.
    !              Eq. 4.22 for the density of states "N0"
    ! -------------------------------------------------------------------------
    !
    !   argument variable
    !
    real(DP)    :: q
    complex(DP) :: z
    integer     :: dim   
    !
    !   local variables
    !
    complex(DP) :: zm,zp,psim,psip
    complex(DP), parameter :: eta=(0.,1.e-16)
    !
    z=z+eta
    !
    if (dim.eq.3) then
       !  
       if (z.eq.0.) then
          !  
          if (q.eq.2.) then
             chi0=-0.5
          else if (q.eq.0.) then
             chi0=-1.
          else
             chi0=-(0.5+(q**2-4.)/(8.*q)*log(abs((q-2.)/(q+2.))))
          end if
          !
       else if (q.ne.0) then
          !
          zm=z/q-0.5*q
          zp=z/q+0.5*q
          !
          psim=psi(zm)
          psip=psi(zp)
          !
          chi0=psim-psip
          chi0=chi0/q
          !
       else if (q.eq.0.) then
          !
          chi0=-1.
          !
       end if
       !
    end if
    !
    return
    !
  end function chi0
  !
  ! ----------------------------------------------------------------------
  complex(DP) function psi(z)
    ! --------------------------------------------------------------------
    implicit none
    !
    !    argument variable
    !
    complex(DP) :: z
    !
    !    local variables
    !
    complex(DP) :: ztemp
    !
    if (z.eq.(1.,0.)) then
       write(6,*) 'ERROR: z=1. in function psi'
       call exit(0)
    end if
    !
    ztemp=(z+1.)/(z-1.)
    psi=0.5*z+0.25*(1-z**2)*log(ztemp)
    !
    return
    !
  end function psi
  !
  ! ---------------------------------------------------------------------
  SUBROUTINE feg_fermi_moment ( )
    ! -------------------------------------------------------------------
    ! 
    !          computation F.E.G. Fermi momentum
    !          kf = (3 * pi^2 * n_elec)^(1/3)
    !
    ! -------------------------------------------------------------------
    USE constants,       ONLY : pi
    USE consts,          ONLY : AngtoBohr
    USE screening_vars,  ONLY : n_elec
    !
    kf = (3._dp * pi ** 2 * n_elec) ** (1._dp/3)
    !
    kf = kf * 1.e-8 / AngtoBohr       ! bohr^-1
    !
  END SUBROUTINE feg_fermi_moment
  !
  ! -----------------------------------------------------------------------
  SUBROUTINE feg_fermi_moment_spol ( )
    ! ---------------------------------------------------------------------
    !
    !      F.E.G. Fermi momentum
    !      kf(is) = (3 * pi^2 * n_elec(is))^(1./3)
    !
    ! ---------------------------------------------------------------------
    USE constants,           ONLY : pi
    USE consts,              ONLY : AngtoBohr
    USE screening_vars,      ONLY : n_elec_up, n_elec_dn
    !
    implicit none
    !
    kf_sp (:) = 0._dp
    kf_sp (1) = (3._dp * pi ** 2 * n_elec_up) ** (1.d0/3)
    kf_sp (2) = (3._dp * pi ** 2 * n_elec_dn) ** (1.d0/3)
    !
    kf_sp (:) = kf_sp (:) * 1.e-8 / AngtoBohr    ! bohr^-1
    !
  END SUBROUTINE feg_fermi_moment_spol
  !
  ! -----------------------------------------------------------------------
  SUBROUTINE set_spol_feg_fermi_energy ( )
    ! ---------------------------------------------------------------------
    !
    !          ef(is) = \hbar^2 / (2 * m_e) * kf(is)^2
    !
    ! ---------------------------------------------------------------------
    USE consts,           ONLY : AngtoBohr, elec_mass, hbar
    !
    implicit none
    !
    ef_sp (:) = 0._dp
    ef_sp (:) = (kf_sp (:) * AngtoBohr * hbar) ** 2 / (2._dp * elec_mass)     ! eV
    !
  END SUBROUTINE set_spol_feg_fermi_energy
  !
  ! -----------------------------------------------------------------------
  real(DP) function feg_fermi_energy ( )
    ! ---------------------------------------------------------------------
    !
    !          computation of F.E.G. Fermi energy
    !          ef = \hbar^2 / (2 * m_e) * kf^2
    !
    ! ---------------------------------------------------------------------
    USE consts,               ONLY : AngtoBohr, elec_mass, hbar
    !
    implicit none
    !
    !          internal variables
    !
    real(DP)                       :: ef       ! fermi energy
    !
    ef = (kf * AngtoBohr * hbar) ** 2 / (2._dp * elec_mass)       ! eV
    !
    feg_fermi_energy = ef
    !
  end function feg_fermi_energy
  !
  !
END MODULE lindhard_screening
!
!
MODULE ks_screening
  ! --------------------------------------------------------------------
  ! 
  !    This module defines all the variables needed
  !    for the KS model of screening
  !
  ! --------------------------------------------------------------------
  
  USE kinds,               ONLY : DP
  
  !
  SAVE
  
  !
  !    global variables
  !
  
  character (LEN=256)           :: fil_chi        ! suscept. input file
  !
  character (LEN=256)           :: fil_chi0       ! suscept. input file
  
  TYPE KS_susceptibility_object
     !
     real (DP)                  :: qp (1:3)       ! q point
     !
     real (DP)                  :: eta            ! broadening
     !
     integer                    :: ng             ! n. G vectors in matrix
     integer                    :: nw             ! n. freq.
     !
     real(DP), allocatable      :: w (:)          ! freq. list
     complex (DP), allocatable  :: chi_stat (:)   ! static suscept.
     complex (DP), allocatable  :: chi_dyn (:,:)  ! dynamic suscept.
     !
  END TYPE KS_susceptibility_object
  !
  TYPE ( KS_susceptibility_object )      :: chi   ! KS susceptibility 
  !
  integer, allocatable          :: igl_to_igg (:)
  !   map from local to global ig
  
  !
  !    G -  q vectors arrays
  !
  INTEGER, ALLOCATABLE           :: ngqp (:)
  !    ngq proc. distr.
  INTEGER, ALLOCATABLE           :: displs (:)
  !    displ. array
  INTEGER                        :: ngq
  !    PWs for q
  INTEGER                        :: ngq_g
  !    global ngq
  
  !
  INTEGER, ALLOCATABLE           :: igq_q (:)
  !!   index G vector corresponding
  !!   to q+G index
  INTEGER, ALLOCATABLE           :: igq_g (:)
  !!   index G vector -> G+q
  !!   global
  real(DP), allocatable          :: gq2 (:)
  !!   local |G+q|^2 array
  real(DP), allocatable          :: gq2_g (:)
  !!   global array
  real(DP), allocatable          :: gq (:,:)
  !!   G+q
  real(DP), allocatable          :: gq_g (:,:)
  !!   G+q global
  
CONTAINS
  !
  !
  SUBROUTINE set_gq_grid ( )
    ! -----------------------------------------------------------------
    !
    !    This subroutine defines the G+q grid
    !    see the GW routine this one is equivalent
    !
    ! ------------------------------------------------------------------
    
    !
    implicit none
    
    !    deallocate igq_q
    
    call deallocate_igq ( )
    
    !
    !    call init_igq
    !
    
    call set_igq ( )
    
    !
    !    set G+q array
    !
    
    call set_gg_plus_q ( )
    
    !
    RETURN
    !
  END SUBROUTINE set_gq_grid
  !
  !
  SUBROUTINE deallocate_igq ( )
    !
    IF ( allocated ( igq_q ) ) DEALLOCATE ( igq_q )
    !
  END SUBROUTINE deallocate_igq
  !
  !
  SUBROUTINE set_igq ( )
    ! ------------------------------------------------------------------
    !      Initialize indices igq and number of plane waves for q point
    !   * (q + G)_i = q + G_igq;
    !   * i = 1, ngq
    !   * igq = igq(i)
    !
    ! ------------------------------------------------------------------
    
    USE constants,         ONLY : eps8
    USE klist,             ONLY : igk_k, nks, xk, ngk
    USE qpoint,            ONLY : xq
    USE gvect,             ONLY : g
    USE wvfct,             ONLY : npwx
    
    !
    implicit none
    !
    !    variables
    !
    
    integer, allocatable        :: igq (:)
    !    igq array local
    integer                     :: err
    integer                     :: ig
    integer                     :: ik       ! k pts iter
    
    !
    ALLOCATE ( igq (1:npwx) )
    igq (:) = 0
    ngq = 0
    
    !
    !    this is called for every q
    !
    
    !    error index
    
    err = 1
    
    !
    !    iterate over k points
    !
    
    do ik= 1, nks
       if ( abs (xk (1,ik) - xq (1)) < eps8 ) then
          if ( abs (xk (2,ik) - xq (2)) < eps8 ) then
             if ( abs (xk (3,ik) - xq (3)) < eps8 ) then
                !
                err = -1
                igq (:) = igk_k (:,ik)
                ngq = ngk (ik)
                !
             end if
          end if
       end if
    end do
    !
    
    IF ( err == 1 ) call errore ('set_igq','xq not found in xk',err)
    
    !
    !    set igq_q array
    !
    
    ALLOCATE ( igq_q (1:ngq) )
    igq_q (:) = 0
    !
    do ig= 1, ngq
       igq_q (ig) = igq (ig)
    end do
    
    !
    IF ( ALLOCATED ( gq2 ) ) DEALLOCATE ( gq2 )
    ALLOCATE ( gq2 (1:ngq) )
    gq2 (:) = 0._dp
    
    !
    !   local |q+G|^2 ordered array
    !
    do ig= 1, ngq
       gq2 (ig) = SUM ( ( xq (:) + g (:, igq_q (ig)) ) ** 2 )
    end do
    !
    RETURN
    !
  END SUBROUTINE set_igq
  !
  !
  SUBROUTINE set_gg_plus_q ( )
    ! ------------------------------------------------------------------
    !
    !      This subroutine executes the same work done
    !      in global_G_grid module but for G+q vector
    !
    ! ------------------------------------------------------------------
    
    USE constants,          ONLY : eps8
    USE mp_images,          ONLY : intra_image_comm, nproc_image, root_image
    USE qpoint,             ONLY : xq
    USE gvect,              ONLY : g
    USE mp_world,           ONLY : world_comm
    USE mp,                 ONLY : mp_bcast, mp_gather, mp_barrier, mp_set_displs
    
    !
    implicit none
    
    !
    !   internal variables
    !
    
    real(DP), allocatable        :: gq_t (:,:)
    real(DP)                     :: gqsq
    !
    integer                      :: i, ig, igg
    integer                      :: ierr
    integer                      :: ntot
    
    !
    !   compute global ngq -> ngq_g
    !
    
    allocate ( ngqp (1:nproc_image), STAT=ierr )
    if (ierr/=0) call errore ('set_gg_plus_q','ngqp not allocatable', ABS(ierr))
    ngqp (:) = 0
    !
    call mp_gather (ngq, ngqp, root_image, intra_image_comm)
    call mp_bcast (ngqp, root_image, world_comm)
    !
    allocate (displs (1:nproc_image), STAT=ierr )
    if (ierr/=0) call errore ('set_gg_plus_q','displs not allocatable', ABS(ierr))
    displs (:) = 0
    call mp_set_displs (ngqp, displs, ntot, nproc_image)
    !
    ngq_g = sum (ngqp)
    
    !
    !   allocate gq_g array
    !
    
    IF ( ALLOCATED ( gq2_g ) ) DEALLOCATE ( gq2_g )
    !
    ALLOCATE ( gq2_g (1:ngq_g), STAT=ierr )
    IF (ierr/=0) call errore ('set_gg_plus_q','gq2_g not allocatable', ABS(ierr))
    !
    gq2_g (:) = 0._dp
    !
    call mp_gather ( gq2 (:), gq2_g (:), ngqp, displs, root_image, intra_image_comm )
    call mp_bcast ( gq2_g, root_image, world_comm )
    call mp_barrier ( intra_image_comm )
    
    !
    !   allocate igq_g global index (G+q)_i -> G_igq_g(i)
    !
    
    IF ( ALLOCATED ( igq_g ) ) DEALLOCATE ( igq_g )
    ALLOCATE ( igq_g (1:ngq_g), STAT=ierr )
    if (ierr/=0) call errore ('set_gg_plus_q','igq_g not allocatable', ABS(ierr))
    !
    igq_g (:) = 0
    
    !
    !   sort array
    !
    
    call hpsort_eps (ngq_g, gq2_g, igq_g, eps8)
    
    !
    !   set G+q global
    !
    
    allocate ( gq (1:3,1:ngq) )
    gq (:,:) = 0._dp
    
    !
    !   local q+G ordered array
    !
    do ig= 1, ngq
       gq (:,ig) = xq (:) + g (:,igq_q (ig))
    end do
    
    !
    IF ( ALLOCATED ( gq_g ) ) DEALLOCATE ( gq_g )
    ALLOCATE ( gq_g (1:3,1:ngq_g) )
    ALLOCATE ( gq_t (1:3,1:ngq_g) )
    gq_t (:,:) = 0._dp
    gq_g (:,:) = 0._dp
    !
    do i= 1, 3
       call mp_gather ( gq (i,:), gq_t (i,:), ngqp, displs, root_image, intra_image_comm )
    end do
    !
    call mp_bcast ( gq_t, root_image, world_comm )
    call mp_barrier ( intra_image_comm )
    
    !
    !   order gq in global array
    !
    
    do igg= 1, ngq_g
       gq_g (:,igg) = gq_t (:,igq_g (igg))
    end do
    
    !
    !    check the global array order
    !
    
    do igg= 1, ngq_g
       !
       gqsq = gq_g (1,igg) ** 2 + gq_g (2,igg) ** 2 + gq_g (3,igg) ** 2
       IF ( ABS (gqsq - gq2_g (igg)) > eps8 ) THEN
          write(6,*) igg, gqsq, gq2_g (igg)
          call errore ('set_gg_plus_q','gqsq != gq2_g(igg)', 1)
       END IF
       !
    end do
    ! 
    
    RETURN
    !
  END SUBROUTINE set_gg_plus_q
  !
  !
  ! --------------------------------------------------------------------
  SUBROUTINE set_ks_screening ( )
    ! ------------------------------------------------------------------
    !
    !    This subroutine reads from file
    !    the pre-computed susceptibility \Chi_{KS}(G+q)
    !    -> non local field effects not here
    !
    ! ------------------------------------------------------------------

    USE constants,    ONLY : eps8
    USE dynvars,      ONLY : dyncalc, en_max, freq_broadening, nfreq
    USE mp_images,    ONLY : intra_image_comm, me_image, nproc_image, root_image
    USE mp_world,     ONLY : world_comm
    USE mp,           ONLY : mp_gather, mp_bcast
    USE gvect,        ONLY : g
    
    !
    implicit none
    
    !
    !    internal variables
    !
    
    character (LEN=256)   :: buffer      ! buffer file
    integer               :: ig, igg     ! G vec. index
    integer               :: ierr
    integer               :: ios
    integer               :: line        ! line file
    integer               :: iw
    LOGICAL               :: exist       ! file stat
    real(DP)              :: r1          ! aux. value
    real(DP)              :: r2          ! aux. value
    real(DP)              :: dw          ! freq. interval
    real(DP)              :: gv (1:3)    ! G vector
    integer               :: ntot
    
    !
    !    allocate susceptibility
    !
    
    chi%ng = 0
    chi%ng = ngq_g
    !
    IF ( dyncalc ) THEN
       !
       chi%eta= freq_broadening
       chi%nw = nfreq
       !
       allocate ( chi%w (1:nfreq), STAT=ierr )
       if (ierr/=0) call errore ('set_ks_screening','allocating w',ABS(ierr))
       chi%w (:) = 0._dp
       !
       dw = en_max / chi%nw           ! eV
       !
       do iw= 1, chi%nw
          chi%w (iw) = (iw - 1) * dw
       end do
       !
    END IF
    !
    IF ( .not. dyncalc ) THEN
       !
       allocate ( chi%chi_stat (1:chi%ng), STAT=ierr )
       if (ierr/=0) call errore ('set_ks_screening','allocating chi_stat',ABS(ierr))
       chi%chi_stat (:) = cmplx (0._dp, 0._dp, kind=DP)
       !
    ELSE
       !
       allocate ( chi%chi_dyn (1:chi%ng,1:chi%nw), STAT=ierr )
       if (ierr/=0) call errore ('set_ks_screening','allocating chi_dyn',ABS(ierr))
       chi%chi_dyn (:,:) = cmplx (0.d0, 0.d0, kind=DP)
       !
    END IF
    
    !
    !    reading data from file
    !
    
    inquire (file = fil_chi, exist = exist)
    !
    IF ( .not. exist ) THEN
       !
       call errore ('set_ks_screening','susceptibility not found...',1)
       !
    ELSE
       !
       open (100, iostat=ios, file=fil_chi, action="read")
       line = 1
       !
       do while (ios == 0)
          !
          read (100, '(A)', iostat=ios) buffer
          !
          IF ( ios == 0 ) THEN
             !
             IF (line == 1) THEN
                !
                READ (buffer,*) chi%qp (1), chi%qp (2), chi%qp (3)   ! qp (tpiba)
                !
             ELSE
                !
                IF ( .not. dyncalc ) THEN
                   !
                   READ (buffer,*) ig, r1, r2                        ! ig= 1, ..., ngq_g
                   chi%chi_stat (ig) = cmplx (r1, r2, kind=DP)
                   !
                ELSE
                   !
                   READ (buffer,*) iw, ig, r1, r2
                   chi%chi_dyn (ig,iw) = cmplx (r1, r2, kind=DP)
                   !
                END IF
                !
             END IF
             !
          END IF
          !
          line = line + 1
          !
       end do
       !
       close (100, iostat=ios)
       !
    END IF
    
    !
    !   define map from local to global G vec.
    !
    
    allocate ( igl_to_igg (1:ngq), STAT=ierr )
    if (ierr/=0) call errore ('set_ks_screening','allocating igl_to_igg',ABS(ierr))
    igl_to_igg (:) = 0
    !
    do ig= 1, ngq
       !
       gv (:) = 0._dp
       gv (:) = chi%qp (:) + g (:,ig)
       !
       do igg= 1, ngq_g
          IF ( ABS(gv (1) - gq_g (1,igg)) < eps8 ) THEN
             IF ( ABS(gv (2) - gq_g (2,igg)) < eps8 ) THEN
                IF ( ABS(gv (3) - gq_g (3,igg)) < eps8 ) THEN
                   igl_to_igg (ig) = igg
                   EXIT
                END IF
             END IF
          END IF
       end do
       !
    end do
    
    !
    RETURN
    !
  END SUBROUTINE set_ks_screening
  !
  !
END MODULE ks_screening
!
!
MODULE local_field_corrections
  !
  !    This module sets up the local field corrections
  !    G (q) to the screening function
  !
  USE kinds,                     ONLY : DP
  USE lindhard_screening,        ONLY : kf_sp
  !
  SAVE
  
  !
  !    global variables
  !
  
  real(DP), allocatable               :: Gss (:,:)
  !
  !    G_ss (q)
  !
CONTAINS
  !
  !
  SUBROUTINE set_local_field_factors ( )
    ! -------------------------------------------------------------------
    !
    !       This subroutine computes G_up,up
    !       G_dn,dn using Hubbard model (see Vignale, electron liquid)
    !       Gup,dn ~ 0
    !       G_ss(q) = q^2 / (q^2 + kf(s)^2)
    !
    ! -------------------------------------------------------------------
    USE lindhard_screening,         ONLY : kf_sp
    USE cell_base,                  ONLY : tpiba2
    USE gvect,                      ONLY : ngl, gl
    USE noncollin_module,           ONLY : npol
    
    !
    implicit none
    !
    !      internal variables
    !
    integer                    :: igl
    integer                    :: ierr
    real(DP)                   :: g2       ! |G|^2
    
    !
    !      allocate local field factors
    !
    allocate ( Gss (1:ngl, 1:npol), STAT=ierr )
    if (ierr/=0) call errore ('set_local_field_factors','allocating Gss',ABS(ierr))
    Gss (:,:) = 0._dp
    !
    !      iterate over |G|
    !
    do igl= 1, ngl
       !
       g2 = gl (igl) * tpiba2    ! |G|^2
       Gss (igl,:) = g2 / (g2 + kf_sp (:) ** 2)
       !
    end do
    !
    return
    !
  END SUBROUTINE set_local_field_factors
  !
  !
  SUBROUTINE set_localq_field_factors ( )
    ! -------------------------------------------------------------------
    !
    !       This subroutine computes G_up,up
    !       G_dn,dn using Hubbard model (see Vignale, electron liquid)
    !       Gup,dn ~ 0
    !       G_ss(q) = q^2 / (q^2 + kf(s)^2)
    !
    ! -------------------------------------------------------------------
    
    USE gvect,              ONLY : ngm, g
    USE cell_base,          ONLY : tpiba2
    USE noncollin_module,   ONLY : npol
    USE qpoint,             ONLY : xq
    
    !
    implicit none
    !
    !   internal variables
    !
    
    integer              :: ierr
    integer              :: ig
    real(DP)             :: g2a         ! |G+q|^2
    
    !
    !   allocate Gss factors
    !
    
    allocate ( Gss (1:ngm, 1:npol), STAT=ierr )
    if (ierr/=0) call errore ('set_localq_field_factors','allocating Gss',ABS(ierr))
    Gss (:,:) = 0._dp
    !
    !   iterate over G+q
    !
    do ig= 1, ngm
       !
       g2a = (xq (1) + g (1,ig)) ** 2 + (xq (2) + g (2,ig)) ** 2 +             &
            (xq (3) + g (3,ig)) ** 2
       !
       g2a = g2a * tpiba2       ! bohr^-2
       !
       Gss (ig,:) = g2a / (g2a + kf_sp (:) ** 2)
       !
    end do
    !
    return
    !
  END SUBROUTINE set_localq_field_factors
  !
  !
END MODULE local_field_corrections
!
!
! ---------------------------------------------
SUBROUTINE set_screening_vars ( )
  ! -------------------------------------------
  !
  USE constants,       ONLY : BOHR_RADIUS_SI
  USE screening_vars,  ONLY : n_elec
  USE klist,           ONLY : nelec
  USE cell_base,       ONLY : omega
  
  implicit none
  !
  !   set elec. density
  !
  n_elec = nelec / omega                          ! elec. density (bohr^-3)
  !
  n_elec = n_elec / (BOHR_RADIUS_SI * 1.E2) ** 3  ! cm^-3
  !
  return
  !
END SUBROUTINE set_screening_vars
!
MODULE vscr_mod
  USE tf_screening
  USE lindhard_screening
  USE ks_screening
  USE local_field_corrections
END MODULE vscr_mod
!
