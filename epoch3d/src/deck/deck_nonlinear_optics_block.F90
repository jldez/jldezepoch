MODULE deck_nonlinear_optics_block
#ifdef NONLINEAR_OPTICS

USE strings_advanced

IMPLICIT NONE
  SAVE

  PRIVATE ! Everything is private unless explicitely specified

  PUBLIC :: nonlinear_optics_deck_initialise
  PUBLIC :: nonlinear_optics_block_handle_element
  PUBLIC :: nonlinear_optics_block_check

  INTEGER i ! Loop counter
  INTEGER, PARAMETER :: n = 26 ! Number of parameters to read + 1
  CHARACTER error_msg(n)*100 ! Array of n strings 1000-character long

CONTAINS

  !----------------------------------------------------------------------------
  ! Initialise
  !----------------------------------------------------------------------------
  SUBROUTINE nonlinear_optics_deck_initialise
    
    error_msg(1) = '*** input deck ERROR in nonlinear optics block ***'
    path_medium_mask = ""
    path_rho_incubation = ""
    resonance = -1.0_num 
    gamma     = 0.0_num
    chi1      = -1.0_num
    chi2      = 0.0_num 
    chi3      = 0.0_num
    gamma_D = -1.0_num
    omega_p_static = -1.0_num
    bandGap_drude = -1.0_num
    avalanche_factor = 0.0_num
    atom_density = -1.0_num
    atom_cross_section = 10.0_num**(-20.0_num)
    effective_mass_electron = 1.0_num
    effective_mass_hole = 1.0_num
    recombination_rate = 0.0_num
    nlo_x_min = x_min
    nlo_x_max = x_max
    nlo_y_min = y_min
    nlo_y_max = y_max
    nlo_z_min = z_min
    nlo_z_max = z_max
    rug_thickness = 0.0_num
    omega_laser = 0.0_num
    amp_laser = 0.0_num

    use_Drude = .FALSE.
    Drude_static = .FALSE.

  END SUBROUTINE nonlinear_optics_deck_initialise


  !---------------------------------------------------------------------------
  ! Handle element function
  !----------------------------------------------------------------------------
  FUNCTION nonlinear_optics_block_handle_element(element, value) RESULT(errcode)
 
    CHARACTER(*), INTENT(IN) :: element, value
    INTEGER :: errcode
    
    errcode = c_err_none
    IF (element == blank .OR. value == blank) RETURN

    IF (str_cmp(element, 'path_medium_mask')) THEN
      path_medium_mask = TRIM(ADJUSTL(value))
      RETURN
    ENDIF

    IF (str_cmp(element, 'path_rho_incubation')) THEN
      path_rho_incubation = TRIM(ADJUSTL(value))
      RETURN
    ENDIF

    IF (str_cmp(element, 'resonance')) THEN
      resonance = as_real_print(value, element, errcode)
      RETURN
    ENDIF

    IF (str_cmp(element, 'gamma')) THEN
      gamma = as_real_print(value, element, errcode)
      RETURN
    ENDIF

    IF (str_cmp(element, 'chi1')) THEN
      chi1 = as_real_print(value, element, errcode)
      RETURN
    ENDIF

    IF (str_cmp(element, 'chi2')) THEN
      chi2 = as_real_print(value, element, errcode)
      RETURN
    ENDIF

    IF (str_cmp(element, 'chi3')) THEN
      chi3 = as_real_print(value, element, errcode)
      RETURN
    ENDIF

    IF (str_cmp(element, 'gamma_D')) THEN
      gamma_D = as_real_print(value, element, errcode)
      RETURN
    ENDIF

    IF (str_cmp(element, 'omega_p_static')) THEN
      omega_p_static = as_real_print(value, element, errcode)
      RETURN
    ENDIF

    IF (str_cmp(element, 'bandGap_drude')) THEN
      bandGap_drude = as_real_print(value, element, errcode)
      RETURN
    ENDIF

    IF (str_cmp(element, 'avalanche_factor')) THEN
      avalanche_factor = as_real_print(value, element, errcode)
      RETURN
    ENDIF

    IF (str_cmp(element, 'atom_density')) THEN
      atom_density = as_real_print(value, element, errcode)
      RETURN
    ENDIF

    IF (str_cmp(element, 'atom_cross_section')) THEN
      atom_cross_section = as_real_print(value, element, errcode)
      RETURN
    ENDIF

    IF (str_cmp(element, 'effective_mass_electron')) THEN
      effective_mass_electron = as_real_print(value, element, errcode)
      RETURN
    ENDIF

    IF (str_cmp(element, 'effective_mass_hole')) THEN
      effective_mass_hole = as_real_print(value, element, errcode)
      RETURN
    ENDIF

    IF (str_cmp(element, 'recombination_rate')) THEN
      recombination_rate = as_real_print(value, element, errcode)
      RETURN
    ENDIF

    IF (str_cmp(element, 'nlo_x_min')) THEN
      nlo_x_min = as_real_print(value, element, errcode)
      RETURN
    ENDIF

    IF (str_cmp(element, 'nlo_x_max')) THEN
      nlo_x_max = as_real_print(value, element, errcode)
      RETURN
    ENDIF

    IF (str_cmp(element, 'nlo_y_min')) THEN
      nlo_y_min = as_real_print(value, element, errcode)
      RETURN
    ENDIF

    IF (str_cmp(element, 'nlo_y_max')) THEN
      nlo_y_max = as_real_print(value, element, errcode)
      RETURN
    ENDIF

    IF (str_cmp(element, 'nlo_z_min')) THEN
      nlo_z_min = as_real_print(value, element, errcode)
      RETURN
    ENDIF

    IF (str_cmp(element, 'nlo_z_max')) THEN
      nlo_z_max = as_real_print(value, element, errcode)
      RETURN
    ENDIF

    IF (str_cmp(element, 'rug_thickness')) THEN
      rug_thickness = as_real_print(value, element, errcode)
      RETURN
    ENDIF

    IF (str_cmp(element, 'omega_laser')) THEN
      omega_laser = as_real_print(value, element, errcode)
      RETURN
    ENDIF

    IF (str_cmp(element, 'amp_laser')) THEN
      amp_laser = as_real_print(value, element, errcode)
      RETURN
    ENDIF

    errcode = c_err_missing_elements

  END FUNCTION nonlinear_optics_block_handle_element

  !---------------------------------------------------------------------------
  ! Block check function
  !----------------------------------------------------------------------------
  FUNCTION nonlinear_optics_block_check() RESULT(errcode)

    LOGICAL :: error
    INTEGER :: j
    INTEGER :: io
    INTEGER :: errcode
    errcode = c_err_none

    !-------------------------------------------------------------------------- 
    ! Error messages are initialized to an empty string
    !--------------------------------------------------------------------------
    DO i = 2, n
      error_msg(i) = "" 
    ENDDO
    
    !-------------------------------------------------------------------------- 
    ! Check input and write error messages if needed
    !--------------------------------------------------------------------------
    error = .FALSE.
    j = 2
    IF (resonance .LT. 0.0) THEN
      error_msg(j) = "resonance is negative or missing"
      error = .TRUE.
      j = j + 1
    ENDIF

    IF (gamma .LT. 0.0) THEN
      error_msg(j) = "gamma should be positive."
      error = .TRUE.
      j = j + 1
    ENDIF

    IF (chi1 .LT. 0.0) THEN
      error_msg(j) = "chi1 is negative or missing"
      error = .TRUE.
      j = j + 1
    ENDIF

    IF (chi2 .LT. 0.0) THEN
      error_msg(j) = "chi2 should be positive"
      error = .TRUE.
      j = j + 1
    ENDIF

    IF (chi3 .LT. 0.0) THEN
      error_msg(j) = "chi3 should be positive"
      error = .TRUE.
      j = j + 1
    ENDIF

    IF (nlo_x_min .GE. nlo_x_max) THEN
      error_msg(j) = "Invalid medium definition (nlo_x_min > nlo_x_max)"
      error = .TRUE.
      j = j + 1
    ENDIF

    IF (nlo_y_min .GE. nlo_y_max) THEN
      error_msg(j) = "Invalid medium definition (nlo_y_min > nlo_y_max)"
      error = .TRUE.
      j = j + 1
    ENDIF

    IF (nlo_z_min .GE. nlo_z_max) THEN
      error_msg(j) = "Invalid medium definition (nlo_z_min > nlo_z_max)"
      error = .TRUE.
      j = j + 1
    ENDIF

    IF (gamma_D .GE. 0.0_num) THEN
      use_Drude = .TRUE.
      IF (omega_p_static .GE. 0.0_num) THEN
        Drude_static = .TRUE.
      ELSE
        IF(bandGap_drude .LT. 0.0) THEN
          error_msg(j) = &
            "You need to specify bandGap_drude to use dynamic Drude"
          error = .TRUE.
          j = j + 1
        ENDIF
        IF(atom_density .LT. 0.0) THEN
          error_msg(j) = &
            "You need to specify atom_density to use dynamic Drude"
          error = .TRUE.
          j = j + 1
        ENDIF
        IF(atom_cross_section .LT. 0.0) THEN
          error_msg(j) = &
            "atom_cross_section has to be positive"
          error = .TRUE.
          j = j + 1
        ENDIF
      ENDIF
    ENDIF

    IF (recombination_rate .LT. 0.0) THEN
      error_msg(j) = "Recombination rate has to be positive"
      error = .TRUE.
      j = j + 1
    ENDIF

    IF (avalanche_factor .LT. 0.0) THEN
      error_msg(j) = "avalanche_factor has to be positive"
      error = .TRUE.
      j = j + 1
    ENDIF

    IF (rug_thickness .LT. 0.0) THEN
      error_msg(j) = "Surface roughness has to be positive"
      error = .TRUE.
      j = j + 1
    ENDIF

    !-------------------------------------------------------------------------- 
    ! If error, print and abort
    !--------------------------------------------------------------------------
    IF (error .EQV. .TRUE.) THEN
      IF (rank .EQ. 0) THEN
        WRITE(io,*) ""
        DO i = 1, n
          IF (error_msg(i) .NE. "") THEN
            WRITE(io,*) error_msg(i)
          ENDIF
        ENDDO
        WRITE(io,*) ""
        errcode = c_err_terminate
      ENDIF
    ENDIF

  END FUNCTION nonlinear_optics_block_check

#endif
END MODULE deck_nonlinear_optics_block
