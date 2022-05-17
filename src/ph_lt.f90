!
!  author:  Jacopo Simoni
!  module file containing the variables to be set for 
!  phonon lifetimes calculations
!
!  ----------------------------------------------------------
!
MODULE phlt_output
  !
  !   ... output files for phonon lifetime calc.
  !
  SAVE
  !
  CHARACTER (LEN=256)      :: filgt
  !    output file G(t)
  CHARACTER (LEN=256)      :: filgammat
  !    output file \gamma_{xy}(t)
  CHARACTER (LEN=256)      :: filgammast
  !    output file \int_0^t dt' \gamma_{xy}(t')
  !
  CHARACTER (LEN=256)      :: fil_ener
  !    file with phonon energies
  CHARACTER (LEN=256)      :: fil_mode
  !    file with phonon modes
  CHARACTER (LEN=256)      :: fil_elph
  !    file e-ph data
  !
  CHARACTER (LEN=256)      :: relaxq_fil
  !    output file relaxation data
  !
  CHARACTER (LEN=256)      :: fil_exclt
  !    file exciton life times
END MODULE phlt_output
!
!
MODULE ph_lt_vars
  !
  USE kinds,    ONLY : DP
  !
  SAVE
  !
  !     set of defined global variables
  !
  real(DP)                  :: ni
  !     ionic density (m^-3)
  real(DP)                  :: t_max
  !     max. time correl. function
  integer                   :: ntd
  !     numb. of time steps
  real(DP), allocatable     :: t (:)
  !     time array
  real(DP)                  :: M
  !     average system's mass (eV fs^2 / bohr^2)
  complex(DP), allocatable  :: fatx_of_q (:)
  !     f.eq(q)
  complex(DP), allocatable  :: fph_ofq_nm (:,:)
  !     <n,k|f(q)|m,k+q>
  complex(DP), allocatable  :: fk_mat (:,:,:,:)
  !     matr. elements
  !     for each perturb.
  real(DP), allocatable     :: focc (:,:)
  real(DP), allocatable     :: focc_e (:,:)
  real(DP), allocatable     :: focc_h (:,:)
  !     occupation probability
  integer                   :: nmodes
  !     numb. modes
  complex(DP), allocatable  :: eq (:,:)
  !     dyn. matrix
  !     used in this module
CONTAINS
  !
  ! --------------------------------------------------------------------------------------
  SUBROUTINE set_fph_ofq ( im )
    ! ------------------------------------------------------------------------------------
    !
    !     CALCULATION POTENTIAL GRADIENT:
    ! f^x_eq(q) = -i\sum_G e^{-iGr}/Omega \sum_R0 eq(q).(G + q) e^{iGR0} e^{iqR0} v_scr(q+G)
    !
    ! ------------------------------------------------------------------------------------
    
    USE constants,           ONLY : tpi
    USE funct,               ONLY : dft_is_gradient, dft_is_nonlocc
    USE ions_base,           ONLY : nat, tau, ityp
    USE cell_base,           ONLY : bg, tpiba
    USE apawpot,             ONLY : vsqloc_ofG
    USE fft_base,            ONLY : dffts, dfftp
    USE fft_interfaces,      ONLY : fwfft, invfft
    USE gvecs,               ONLY : ngms
    USE gvect,               ONLY : ngm, mill, g
    USE qpoint,              ONLY : xq
    USE uspp_param,          ONLY : upf
    USE uspp,                ONLY : nlcc_any
    USE nlcc_ph,             ONLY : drc
    USE eqv,                 ONLY : dmuxc
    USE scf,                 ONLY : rho, rho_core
    USE gc_lr,               ONLY : dvxc_rr, dvxc_s, dvxc_sr, dvxc_ss, grho
    USE lsda_mod,            ONLY : lsda, nspin
    USE noncollin_module,    ONLY : nspin_gga
    USE modes,               ONLY : u
    USE Coul_cut_2D,         ONLY : do_cutoff_2D
    USE Coul_cut_2D_ph,      ONLY : cutoff_localq
    
    !
    IMPLICIT NONE
    !
    !     internal variables
    !

    integer, intent(in)           :: im
    !     input mode
    integer                       :: ia, ir, mu
    !     atom's index / r / mode
    integer                       :: id, ig
    !     direction's index / g
    integer                       :: n1, n2, n3, nr1, nr2, nr3, nt
    integer                       :: ierr
    !     error index
    complex(DP), allocatable      :: aux1 (:), drhoc (:,:)
    complex(DP), allocatable      :: aux (:,:)
    complex(DP), pointer          :: auxs (:)
    !     aux. vars.
    complex(DP), allocatable      :: str_fac1 (:,:), str_fac2 (:,:), str_fac3 (:,:)
    complex(DP), allocatable      :: str_facq (:)
    !     str. factors
    complex(DP)                   :: u1, u2, u3
    !     u (i) comp.
    real(DP)                      :: bgtau (3), arg
    !     G R0
    complex(DP)                   :: fac, fact, gtau
    !     factors
    complex(DP)                   :: gu, gu0
    !     u.G
    
    ! -----------------------------------------------------------------
    !    if PAW calculation consider
    !    only valence electrons : this part is not implemented yet
    ! -----------------------------------------------------------------
    
    !
    !    start calculation fx_of_q
    !    fx(q) = -i \sum_G e^{-iGr}/Omega \sum_R0 eq(q).(G + q) e^{iGR0} e^{iqR0} v_scr(q+G)
    !
    
    IF ( .NOT. ALLOCATED ( fatx_of_q ) ) THEN
       allocate ( fatx_of_q ( 1:dffts%nnr ), STAT=ierr )
       if (ierr/=0) call errore('set_fph_ofq','allocating fatx_of_q',ABS(ierr))
    END IF
    !
    fatx_of_q (:) = cmplx( 0.d0, 0.d0, kind=DP )
    
    !
    !      allocate aux1
    !
    
    allocate ( aux1 (1:dffts%nnr), STAT=ierr )
    if (ierr/=0) call errore ('set_fph_ofq','allocating aux1',ABS(ierr))
    
    !
    !      the screened potential vie_loc(q+G)
    !      was already computed
    !
    
    fac = cmplx (0.d0, -1.d0, kind=DP) * tpiba
    
    !
    !      define structure factor e^{iqR0}
    !
    
    allocate ( str_facq (1:nat) )
    !
    do ia= 1, nat
       !
       arg = ( xq (1) * tau (1,ia) +        &
            xq (2) * tau (2,ia) +           &
            xq (3) * tau (3,ia) ) * tpi
       !
       str_facq (ia) = cmplx ( cos (arg), sin (arg), kind=DP )
       !
    end do
    
    !
    !      define structure factor e^{iG R0}
    !
    
    nr1 = dfftp%nr1
    nr2 = dfftp%nr2
    nr3 = dfftp%nr3
    !
    allocate ( str_fac1 (-nr1:nr1, 1:nat), str_fac2 (-nr2:nr2, 1:nat),          &
         str_fac3 (-nr3:nr3, 1:nat) )
    !
    str_fac1 (:,:) = cmplx (0.d0, 0.d0, kind=DP)
    str_fac2 (:,:) = cmplx (0.d0, 0.d0, kind=DP)
    str_fac3 (:,:) = cmplx (0.d0, 0.d0, kind=DP)
    !
    do ia= 1, nat
       !
       do id= 1, 3
          bgtau (id) = bg (1,id) * tau (1, ia) +        &
               bg (2,id) * tau (2, ia) +                &
               bg (3,id) * tau (3, ia)
       end do
       !
       do n1= -nr1, nr1
          arg = tpi * n1 * bgtau (1)
          str_fac1 (n1,ia) = cmplx ( cos(arg), sin(arg), kind=DP )
       end do
       !
       do n2= -nr2, nr2
          arg = tpi * n2 * bgtau (2)
          str_fac2 (n2,ia) = cmplx ( cos(arg), sin(arg), kind=DP )
       end do
       !
       do n3= -nr3, nr3
          arg = tpi * n3 * bgtau (3)
          str_fac3 (n3,ia) = cmplx ( cos(arg), sin(arg), kind=DP )
       end do
       !
    end do
    
    !
    !   check non local core corrections
    !
    
    write(6,*) "nlcc= ", nlcc_any
    !
    IF (nlcc_any) THEN
       !
       allocate ( drhoc (1:dfftp%nnr, 1:nspin), STAT=ierr )
       if (ierr/=0) call errore ('set_fph_ofq','allocating drhoc',ABS(ierr))
       !
    END IF
    !
    aux1 (:) = cmplx (0.d0, 0.d0, kind=DP)
    !
    
    !
    !      compute local contribution for each atom
    !
    
    DO ia= 1, nat
       !
       fact = fac * str_facq (ia)
       mu = 3 * (ia - 1)
       !
       IF ( abs (u (mu+1, im)) + abs (u (mu+2, im)) + abs (u (mu+3, im)) .GT. 1.0d-12 ) THEN
          !
          nt = ityp (ia)
          u1 = u (mu+1, im)
          u2 = u (mu+2, im)
          u3 = u (mu+3, im)
          !
          gu0 = xq (1) * u1 + xq (2) * u2 + xq (3) * u3
          !
          !    run over PWs
          !
          do ig= 1, ngms
             !
             gtau = str_fac1 (mill(1,ig), ia) * str_fac2 (mill(2,ig), ia) *          &
                  str_fac3 (mill(3,ig), ia)
             !
             gu = gu0 + g (1,ig) * u1 + g (2,ig) * u2 + g (3,ig) * u3
             aux1 ( dffts%nl (ig) ) = aux1 ( dffts%nl (ig) ) +                       &
                  vsqloc_ofG (ig, nt) * gu * fact * gtau
             !
          end do
          !
          IF (do_cutoff_2D) THEN
             call cutoff_localq (aux1, fact, u1, u2, u3, gu0, nt, ia)
          END IF
          !
       END IF
       !
    END DO
    
    !
    !    add NLCC if present
    !
    
    IF (nlcc_any) THEN
       !
       drhoc (:,:) = cmplx(0.d0, 0.d0, kind=DP)
       !
       allocate ( aux (1:dfftp%nnr, 1:nspin), STAT=ierr )
       if (ierr/=0) call errore ('set_fph_ofq','allocating aux',ABS(ierr))
       aux (:,:) = cmplx (0.d0, 0.d0, kind=DP)
       !
       DO ia= 1, nat
          !
          fact = fac * str_facq (ia)
          mu = 3 * (ia - 1)
          !
          IF ( abs (u (mu+1, im)) + abs (u (mu+2, im)) + abs (u (mu+3, im)) .GT. 1.0d-12 ) THEN
             !
             nt = ityp (ia)
             u1 = u (mu+1, im)
             u2 = u (mu+2, im)
             u3 = u (mu+3, im)
             !
             gu0 = xq (1) * u1 + xq (2) * u2 + xq (3) * u3
             !
             IF ( upf(ityp(ia))%nlcc ) THEN
                !
                do ig= 1, ngm
                   !
                   gtau = str_fac1(mill(1,ig), ia) * str_fac2(mill(2,ig), ia) *            &
                        str_fac3(mill(3,ig), ia)
                   !
                   gu = gu0 + g (1,ig) * u1 + g (2,ig) * u2 + g (3,ig) * u3
                   !
                   drhoc(dfftp%nl (ig),1) = drhoc(dfftp%nl (ig),1) + drc(ig,nt) *          &
                        fact * gtau * gu
                   !
                end do
                !
             END IF
          END IF
          !
       END DO
       !
       call invfft ('Rho', drhoc(:,1), dfftp)
       !
       IF ( .not. lsda ) THEN
          do ir= 1, dfftp%nnr
             aux (ir,1) = drhoc (ir,1) * dmuxc (ir,1,1)
          end do
       ELSE
          !
          !   spin unpolarized calc not implemented here
          !
          call errore ('set_fph_ofq','spin polarized not implemented', 1)
          !
       END IF
       !
       rho%of_r (:,1) = rho%of_r (:,1) + rho_core (:)
       !
       IF ( dft_is_gradient ( ) ) call dgradcorr (dfftp, rho%of_r, grho, dvxc_rr,      &
            dvxc_sr, dvxc_ss, dvxc_s, xq, drhoc, nspin, nspin_gga, g, aux)
       !
       IF (dft_is_nonlocc()) call dnonloccorr (rho%of_r, drhoc, xq, aux)
       deallocate (drhoc)
       !
       rho%of_r (:,1) = rho%of_r (:,1) - rho_core (:)
       !
       call fwfft ('Rho', aux (:,1), dfftp)
       ! 
       !  This is needed also when the smooth and the thick grids coincide to
       !  cut the potential at the cut-off
       !
       allocate ( auxs (1:dffts%nnr), STAT=ierr )
       if (ierr/=0) call errore ('set_fph_ofq','allocating auxs',ABS(ierr))
       auxs (:) = cmplx(0.d0, 0.d0, kind=DP)
       !
       do ig= 1, ngms
          auxs (dffts%nl (ig)) = aux (dfftp%nl (ig), 1)
       end do
       !
       aux1 (:) = aux1 (:) + auxs (:)
       !
       deallocate ( aux, auxs )
       !
    END IF
    
    !
    !    aux1 -> fatx_of_q : in real space
    !
    
    call invfft ('Rho', aux1, dffts)
    !
    fatx_of_q (:) = aux1 (:)                       ! Ry / bohr
    !
    
    return
    !
  END SUBROUTINE set_fph_ofq
  !
  !
  ! ---------------------------------------------------------------------------------
  SUBROUTINE set_eq_matr ( )
    ! -------------------------------------------------------------------------------
    !
    !     This subroutine set eq matrix equivalent to dyn
    !     but with our units : bohr / fs * eV^{-1/2}
    !
    ! -------------------------------------------------------------------------------
    
    USE constants,          ONLY : amu_ry
    USE consts,             ONLY : proton_mass
    USE ions_base,          ONLY : nat, ityp, amass
    USE dynmat,             ONLY : dyn
    
    !
    implicit none
    
    !
    !     internal variables
    !
    
    integer                      :: ia
    !     atom's index
    integer                      :: mu, nu
    !     mode's index
    integer                      :: ierr
    
    !
    !     allocate eq matrix 
    !
    
    IF ( .not. allocated ( eq ) ) THEN
       allocate ( eq (1:nmodes, 1:nmodes), STAT=ierr )
       if (ierr/=0) call errore ('set_eq_matr','allocating eq',ABS(ierr))
    END IF
    !
    eq (:,:) = cmplx (0.d0, 0.d0, kind=DP)
    
    !
    !   get dyn. matrix from dyn
    !
    
    do nu= 1, 3 * nat
       do mu= 1, 3 * nat
          ia = (mu - 1) / 3 + 1
          eq (mu, nu) = dyn (mu, nu) * SQRT ( amu_ry * amass (ityp (ia)) )
          !
          !   modes are now normalized
          !
       end do
    end do
    
    !
    !   & multiply 1/sqrt(M)
    !
    
    do nu= 1, 3 * nat
       do mu= 1, 3 * nat
          ia = (mu - 1) / 3 + 1
          eq (mu, nu) = eq (mu, nu) / SQRT ( amass (ityp(ia)) * proton_mass )
          !
          !   displacement (bohr / fs eV^-1/2)
          !
       end do
    end do
    
    !
    RETURN
    !
  END SUBROUTINE set_eq_matr
  !
  !
  ! ---------------------------------------------------------------------------------
  SUBROUTINE set_fph_mtxel ( ik )
    ! -------------------------------------------------------------------------------
    !
    !  This subroutine computes the f_eq(q) matrix elements
    !  valid with local pseudo potential & first contribution with PAW potentials
    !
    !  f_{n,k; m,k+q} = \delta(k',k+q) 1/Omega \int_\Omega dr u_{n,k}(r)*
    !            f_eq(q) u_{m,k+q}(r)
    !
    ! -------------------------------------------------------------------------------
    
    USE constants,                   ONLY : tpi
    USE buffers,                     ONLY : get_buffer
    USE fft_base,                    ONLY : dffts
    USE klist,                       ONLY : ngk, igk_k, xk
    USE noncollin_module,            ONLY : noncolin, npol
    USE cell_base,                   ONLY : tpiba, omega
    USE lsda_mod,                    ONLY : current_spin, lsda, isk
    USE qe_grid,                     ONLY : dv
    USE units_lr,                    ONLY : iuwfc, lrwfc
    USE control_lr,                  ONLY : lgamma
    USE wvfct,                       ONLY : nbnd, npwx
    USE qpoint,                      ONLY : nksq, ikks, ikqs
    USE uspp,                        ONLY : vkb
    USE mp,                          ONLY : mp_sum
    USE mp_world,                    ONLY : world_comm
    USE wavefunctions,               ONLY : evc
    USE eqv,                         ONLY : evq
    USE gvect,                       ONLY : g, ngm
    USE fft_interfaces,              ONLY : invfft
    
    !
    implicit none
    !
    !    internal variables
    !
    integer, intent(in)                   :: ik
    !    k-point counter
    integer                               :: ikk
    !    index of k-point vector k
    integer                               :: npw
    !    numb. PW for vector k
    integer                               :: npwq
    !    numb. PW for vector k+q
    integer                               :: ikq
    !    index of k-point vector k+q
    integer                               :: ib1, ib2
    !    band indexes
    integer                               :: ierr
    !    error index
    integer                               :: ir
    !    r index
    integer                               :: ig
    !    plane wave index
    complex(DP), allocatable              :: evc_r (:,:), evq_r (:,:)
    !    real space wave functions
    complex(DP)                           :: faux_r
    !    auxiliary matr. el. array
    
    !
    !    allocate wave functions
    !
    
    ALLOCATE ( evc_r (1:dffts%nnr, 1:npol), STAT=ierr )
    if (ierr/=0) call errore ('set_fph_mtxel','allocating evc_r',ABS(ierr))
    !
    ALLOCATE ( evq_r (1:dffts%nnr, 1:npol), STAT=ierr )
    if (ierr/=0) call errore ('set_fph_mtxel','allocating evq_r',ABS(ierr))
    
    !
    !    allocate matrix elements
    !
    
    IF ( .NOT. ALLOCATED ( fph_ofq_nm ) ) THEN
       ALLOCATE ( fph_ofq_nm (1:nbnd,1:nbnd), STAT=ierr )
       if (ierr/=0) call errore ('set_fph_mtxel','allocating fph_ofq_nm',ABS(ierr))
    END IF
    !
    fph_ofq_nm (:,:) = cmplx ( 0.d0, 0.d0, kind=DP )
    
    !
    !    start k-point calculation
    !
    
    !   ik = counter of k-points with vector k
    !   ikk= index of k-point with vector k
    !   ikq= index of k-point with vector k+q
    !        k and k+q are alternated if q!=0, the same if q=0
    
    ikk = ikks (ik)
    ikq = ikqs (ik)
    
    !
    IF ( lsda ) current_spin = isk (ikk)
    !
    npw = ngk (ikk)
    npwq= ngk (ikq)
    !
    call init_us_2 (npwq, igk_k(1,ikq), xk (1, ikq), vkb)
    !
    !   read unperturbed wave functions psi(k) and psi(k+q)
    !
    IF ( nksq .GT. 1 ) THEN
       IF ( lgamma ) THEN
          call get_buffer (evc, lrwfc, iuwfc, ikk)
       ELSE
          call get_buffer (evc, lrwfc, iuwfc, ikk)
          call get_buffer (evq, lrwfc, iuwfc, ikq)
       END IF
    END IF
    
    !
    !   run over bands index
    !
    
    do ib1= 1, nbnd
       !
       !   perform Fourier transform of the wave functions
       !
       evc_r (:,:) = cmplx ( 0.d0, 0.d0, kind=DP )
       !
       do ig= 1, npw
          evc_r (dffts%nl (igk_k(ig,ikk)), 1) = evc (ig,ib1)
       end do
       !
       call invfft ('Wave',evc_r (:,1), dffts)
       !
       IF (noncolin) THEN
          do ig= 1, npw
             evc_r (dffts%nl (igk_k(ig,ikk)), 2) = evc (ig+npwx, ib1)
          end do
          !
          call invfft ('Wave', evc_r (:,2), dffts)
          !
       END IF
       !
       do ib2= 1, nbnd
          !
          !    perform Fourier transform psi(q)
          !
          evq_r (:,:) = cmplx ( 0.d0, 0.d0, kind=DP )
          !
          do ig= 1, npwq
             evq_r (dffts%nl (igk_k(ig,ikq)), 1) = evq (ig,ib2)
          end do
          !
          call invfft ('Wave', evq_r (:,1), dffts)
          !
          IF (noncolin) THEN
             do ig= 1, npwq
                evq_r (dffts%nl (igk_k(ig,ikq)), 2) = evq (ig+npwx, ib2)
             end do
             !
             call invfft ('Wave', evq_r (:,2), dffts)
             !
          END IF
          !
          !    compute the matrix elements: f_{n,k; m,k+q}
          !
          faux_r = cmplx ( 0.d0, 0.d0, kind=DP )
          !
          do ir= 1, dffts%nnr
             !
             faux_r = conjg( evc_r (ir,1) ) * fatx_of_q (ir) * evq_r (ir,1) +       &
                  faux_r
             !
             !   with non collinear calculations
             !
             IF (noncolin) THEN
                faux_r = conjg( evc_r (ir,2) ) * fatx_of_q (ir) * evq_r (ir,2) +    &
                     faux_r
             END IF
             !
          end do
          !
          call mp_sum ( faux_r, world_comm )
          !
          faux_r = faux_r * dv / omega
          !
          !   define matrix elements
          !
          fph_ofq_nm (ib1,ib2) = faux_r
          !
          !   eV / bohr
          !
       end do
       !
    end do
    !
    return
    !
  END SUBROUTINE set_fph_mtxel
  !
  !
  ! ---------------------------------------------------------------------------------
  SUBROUTINE set_fph_mtxel_rec_space ( ik, im )
    ! -------------------------------------------------------------------------------
    !
    !  This subroutine computes the f_eq(q) matrix elements
    !  valid with local pseudo potential & first contribution with PAW potentials
    !
    !  f_{j,k+q; i,k} = \delta(k',k+q) <psi_{k',j}|dv_q*psi_{k,i}>
    !
    ! -------------------------------------------------------------------------------
    
    USE constants,                   ONLY : RYTOEV
    USE buffers,                     ONLY : get_buffer
    USE fft_base,                    ONLY : dffts
    USE klist,                       ONLY : ngk, igk_k, xk
    USE noncollin_module,            ONLY : noncolin, npol
    USE lsda_mod,                    ONLY : current_spin, lsda, isk
    USE units_lr,                    ONLY : iuwfc, lrwfc
    USE control_lr,                  ONLY : lgamma
    USE wvfct,                       ONLY : nbnd, npwx
    USE qpoint,                      ONLY : nksq, ikks, ikqs
    USE uspp,                        ONLY : vkb
    USE mp,                          ONLY : mp_sum
    USE mp_bands,                    ONLY : intra_bgrp_comm
    USE wavefunctions,               ONLY : evc
    USE eqv,                         ONLY : evq, dvpsi
    USE modes,                       ONLY : u
    USE gvect,                       ONLY : g, ngm
    USE fft_interfaces,              ONLY : invfft, fwfft
    
    !
    implicit none
    !
    !    internal variables
    !
    integer, intent(in)                   :: ik
    !    k-point counter
    integer, intent(in)                   :: im
    !    mode index
    integer                               :: ikk
    !    index of k-point vector k
    integer                               :: npw
    !    numb. PW for vector k
    integer                               :: npwq
    !    numb. PW for vector k+q
    integer                               :: ikq
    !    index of k-point vector k+q
    integer                               :: ib1, ib2
    !    band indexes
    integer                               :: ierr
    !    error index
    integer                               :: ir
    !    r index
    integer                               :: ig
    !    plane wave index
    complex(DP), allocatable              :: evc_r (:,:)
    !   wfc_k
    complex(DP), allocatable              :: aux_r (:,:)
    complex(DP)                           :: faux
    !   aux. array
    complex(DP), external                 :: zdotc
    
    !
    !    allocate wave functions
    !
    
    ALLOCATE ( evc_r (1:dffts%nnr, 1:npol), STAT=ierr )
    if (ierr/=0) call errore ('set_fph_mtxel','allocating evc_r',ABS(ierr))
    !
    ALLOCATE ( aux_r (1:dffts%nnr, 1:npol), STAT=ierr )
    if (ierr/=0) call errore ('set_fph_mtxel','allocating aux_r',ABS(ierr))
    !
    
    !
    !    allocate matrix elements
    !
    
    IF ( .NOT. ALLOCATED ( fph_ofq_nm ) ) THEN
       ALLOCATE ( fph_ofq_nm (1:nbnd,1:nbnd), STAT=ierr )
       if (ierr/=0) call errore ('set_fph_mtxel','allocating fph_ofq_nm',ABS(ierr))
    END IF
    !
    fph_ofq_nm (:,:) = cmplx ( 0.d0, 0.d0, kind=DP )
    
    !
    !    start k-point calculation
    !
    
    !   ik = counter of k-points with vector k
    !   ikk= index of k-point with vector k
    !   ikq= index of k-point with vector k+q
    !        k and k+q are alternated if q!=0, the same if q=0
    
    ikk = ikks (ik)
    ikq = ikqs (ik)
    
    !
    IF ( lsda ) current_spin = isk (ikk)
    !
    npw = ngk (ikk)
    npwq= ngk (ikq)
    !
    call init_us_2 (npwq, igk_k(1,ikq), xk (1, ikq), vkb)
    !
    !   read unperturbed wave functions psi(k) and psi(k+q)
    !
    IF ( nksq .GT. 1 ) THEN
       IF ( lgamma ) THEN
          call get_buffer (evc, lrwfc, iuwfc, ikk)
       ELSE
          call get_buffer (evc, lrwfc, iuwfc, ikk)
          call get_buffer (evq, lrwfc, iuwfc, ikq)
       END IF
    END IF
    
    !
    dvpsi (:,:) = cmplx (0.d0, 0.d0, kind=DP)
    !
    
    !
    !   run over bands index
    !
    
    do ib1= 1, nbnd
       !
       !   perform Fourier transform of the wave functions
       !
       evc_r (:,:) = cmplx (0.d0, 0.d0, kind=DP)
       !
       do ig= 1, npw
          evc_r (dffts%nl (igk_k (ig,ikk)), 1) = evc (ig,ib1)
       end do
       !
       call invfft ('Wave',evc_r (:,1), dffts)
       !
       IF (noncolin) THEN
          !
          do ig= 1, npw
             evc_r (dffts%nl (igk_k(ig,ikk)), 2) = evc (ig+npwx, ib1)
          end do
          !
          call invfft ('Wave', evc_r (:,2), dffts)
          !
       END IF
       !
       !  compute dVloc/dtau * psi in real space
       !
       aux_r (:,:) = cmplx (0.d0, 0.d0, kind=DP)
       !
       do ir= 1, dffts%nnr
          aux_r (ir,:) = fatx_of_q (ir) * evc_r (ir,:)
       end do
       !
       !  now transform dVloc/dtau * psi in rec. space
       !
       call fwfft ('Wave', aux_r (:,1), dffts)
       !
       do ig= 1, npwq
          dvpsi (ig,ib1) = aux_r (dffts%nl (igk_k (ig,ikq)), 1 ) +           &
               dvpsi (ig,ib1)
       end do
       !
       IF (noncolin) THEN
          !
          call fwfft ('Wave', aux_r (:,2), dffts)
          !
          do ig= 1, npwq
             dvpsi (ig+npwx,ib1) = aux_r ( dffts%nl (igk_k (ig,ikq)), 2 ) +  &
                  dvpsi (ig+npwx,ib1)
          end do
          !
       END IF
       !
    end do   ! end band index 1
    
    !  we add here the contribution of the nonlocal potential
    !  in the US form
    !
    
    call dvqpsi_us_only (ik, u (1,im))
    
    !
    !  run over band index 1
    !
    
    do ib1= 1, nbnd
       
       !
       !  run over band 2
       !
       
       do ib2= 1, nbnd
          !
          !    compute the matrix elements: f_{j,k+q; i,k}
          !
          faux = cmplx ( 0.d0, 0.d0, kind=DP )
          faux = zdotc (npwq, evq (1,ib2), 1, dvpsi (1,ib1), 1)
          !
          !   with non collinear calculations
          !
          IF (noncolin) THEN
             faux = zdotc (npwq, evq(npwx+1,ib2), 1, dvpsi(npwx+1,ib1), 1) +    &
                  faux
          END IF
          !
          call mp_sum ( faux, intra_bgrp_comm )
          !
          !   define matrix elements
          !
          fph_ofq_nm (ib1,ib2) = CONJG (faux) * RYTOEV
          !
          !   eV / bohr
          !
       end do
       !
    end do
    
    !
    return
    !
  END SUBROUTINE set_fph_mtxel_rec_space
  !
  !
  ! ---------------------------------------------------------------------------
  SUBROUTINE set_fk_mat ( )
    ! -------------------------------------------------------------------------
    !
    !    This subroutine does something equivalent to elphel
    !    it computes fk_mat for the different perturbation modes
    !    for a given q
    ! f_{n,k; m,k+q}^p(k) = \delta(k',k+q) 1/Omega \int_\Omega dr u_{n,k}(r)*
    !            f_eq^p(q) u_{m,k+q}(r)
    !
    ! -------------------------------------------------------------------------
    
    USE qpoint,                    ONLY : nksq, ikks, ikqs
    USE wvfct,                     ONLY : nbnd
    USE mp_pools,                  ONLY : npool
    USE ions_base,                 ONLY : nat
    
    !
    implicit none
    
    !
    !    internal variables
    !
    
    integer                            :: nksqtot
    !    tot. k points
    integer                            :: ik, ikk, ikq
    !    k-pt. counters
    integer                            :: ipert
    !    pert. count
    integer                            :: ierr
    
    !
    !    allocate collected variables
    !
    
    !
    IF (npool==1) THEN
       !
       !   np pool, just copy old variables on the new ones
       !
       nksqtot = nksq
       !
       IF ( ALLOCATED ( fk_mat ) ) DEALLOCATE ( fk_mat )
       allocate ( fk_mat (1:nbnd, 1:nbnd, 1:nksqtot, 1:3*nat), STAT=ierr )
       if (ierr/=0) call errore ('set_fk_mat','allocating fk_mat', ABS(ierr))
       fk_mat (:,:,:,:) = cmplx (0.d0, 0.d0, kind=DP)
       !
    ELSE
       !
       call errore ('set_fk_mat','npool > 1 not implemented', 1)
       !
    END IF
    !
    
    !
    !   run over pert. modes
    !
    
    DO ipert= 1, 3*nat
       
       !
       !   compute directional derivative
       !   f^q_e(q)(r)
       !
       
       call set_fph_ofq ( ipert )
       
       !
       !   run over k points loop
       !
       
       DO ik= 1, nksqtot
          !
          !   ik = counter of k-points with vector k
          !   ikk= index of k-point with vector k
          !   ikq= index of k-point with vector k+q
          !         k and k+q are alternated if q != 0, are the same if q=0
          
          ikk = ikks (ik)
          ikq = ikqs (ik)
          
          !
          !   compute f^x_{nk;mk+q} matrix elements
          !
          
          call set_fph_mtxel_rec_space ( ik, ipert )
          
          !
          !   set matrix array
          !
          
          fk_mat (:,:,ik,ipert) = fph_ofq_nm (:,:)
          
          !
          !   eV / bohr  -- crystal basis
          !
       END DO
       
       !
    END DO
    
    !
    return
    !
  END SUBROUTINE set_fk_mat
  !
  !
  !
END MODULE ph_lt_vars
!
!
MODULE ph_lt
  USE phlt_output
  USE ph_lt_vars
END MODULE ph_lt
