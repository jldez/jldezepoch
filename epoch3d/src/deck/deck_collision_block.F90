! Copyright (C) 2010-2015 Keith Bennett <K.Bennett@warwick.ac.uk>
! Copyright (C) 2009      Chris Brady <C.S.Brady@warwick.ac.uk>
!
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.

MODULE deck_collision_block

  USE strings_advanced
  USE collisions

  IMPLICIT NONE
  SAVE

  PRIVATE
  PUBLIC :: collision_deck_initialise, collision_deck_finalise
  PUBLIC :: collision_block_start, collision_block_end
  PUBLIC :: collision_block_handle_element, collision_block_check

  LOGICAL, ALLOCATABLE, DIMENSION(:,:) :: coll_pairs_touched

CONTAINS

  SUBROUTINE collision_deck_initialise

    IF (deck_state == c_ds_first) THEN
      use_collisions = .FALSE.
      use_collisional_ionisation = .FALSE.
    ELSE
      ALLOCATE(coll_pairs_touched(1:n_species, 1:n_species))
      coll_pairs_touched = .FALSE.
      CALL setup_collisions
    ENDIF

  END SUBROUTINE collision_deck_initialise



  SUBROUTINE collision_deck_finalise

    INTEGER :: i, j

    IF (deck_state == c_ds_first) RETURN
    DEALLOCATE(coll_pairs_touched)

    IF (use_collisions) THEN
      use_collisions = .FALSE.
      DO j = 1, n_species
        DO i = 1, n_species
          IF (coll_pairs(i,j) > 0) THEN
            use_collisions = .TRUE.
            EXIT
          ENDIF
        ENDDO
      ENDDO
      use_particle_lists = use_particle_lists .OR. use_collisions
      need_random_state = .TRUE.
    ENDIF

  END SUBROUTINE collision_deck_finalise



  SUBROUTINE collision_block_start

  END SUBROUTINE collision_block_start



  SUBROUTINE collision_block_end

  END SUBROUTINE collision_block_end



  FUNCTION collision_block_handle_element(element, value) RESULT(errcode)

    CHARACTER(*), INTENT(IN) :: element, value
    INTEGER :: errcode

    errcode = c_err_none
    IF (element == blank .OR. value == blank) RETURN

    ! Performed on second parse to ensure that species are set up first.
    IF (str_cmp(element, 'collide')) THEN
      IF (deck_state /= c_ds_first) THEN
        CALL set_collision_matrix(TRIM(ADJUSTL(value)), errcode)
      ENDIF
      RETURN
    ENDIF

    IF (str_cmp(element, 'use_collisions')) THEN
      use_collisions = as_logical_print(value, element, errcode)
      RETURN
    ENDIF

    IF (str_cmp(element, 'coulomb_log')) THEN
      IF (str_cmp(value, 'auto')) THEN
        coulomb_log_auto = .TRUE.
      ELSE
        coulomb_log_auto = .FALSE.
        coulomb_log = as_real_print(value, element, errcode)
      ENDIF
      RETURN
    ENDIF

    IF (str_cmp(element, 'collisional_ionisation')) THEN
      use_collisional_ionisation = as_logical_print(value, element, errcode)
      IF (use_collisional_ionisation) THEN
#ifndef PER_SPECIES_WEIGHT
        use_collisions = .TRUE.
#else
        IF (as_logical_print(value, element, errcode)) THEN
          errcode = c_err_pp_options_wrong
          extended_error_string = '-DPER_SPECIES_WEIGHT'
        ENDIF
        use_collisional_ionisation = .FALSE.
#endif
      ENDIF
      RETURN
    ENDIF

    errcode = c_err_unknown_element

  END FUNCTION collision_block_handle_element



  FUNCTION collision_block_check() RESULT(errcode)

    INTEGER :: errcode

    errcode = c_err_none

  END FUNCTION collision_block_check



! The following code is all about reading the coll_pairs from the input deck

  SUBROUTINE get_token(str_in, str_out, token_out, err)

    CHARACTER(*), INTENT(IN) :: str_in
    CHARACTER(*), INTENT(OUT) :: str_out
    CHARACTER(*), INTENT(OUT) :: token_out
    INTEGER, INTENT(INOUT) :: err
    INTEGER :: str_len, char, pos
    CHARACTER(1) :: c

    str_len = LEN(str_in)
    pos = str_len

    DO char = 1, str_len
      c = str_in(char:char)
      IF (c == ' ')  THEN
        pos = char
        EXIT
      ENDIF
    ENDDO

    IF (pos < str_len) THEN
      str_out = TRIM(ADJUSTL(str_in(pos+1:str_len)))
    ELSE
      str_out = ''
    ENDIF

    token_out = TRIM(str_in(1:pos))

  END SUBROUTINE get_token



  SUBROUTINE set_collision_matrix(str_in, errcode)

    CHARACTER(*), INTENT(IN) :: str_in
    INTEGER, INTENT(INOUT) :: errcode
    CHARACTER(LEN=string_length) :: tstr1, tstr2
    CHARACTER(LEN=string_length) :: species1, species2
    REAL(num) :: collstate
    INTEGER :: io, iu, sp1, sp2

    IF (deck_state /= c_ds_last) RETURN

    IF (str_cmp(TRIM(str_in), 'all')) THEN
      coll_pairs = 1.0_num
      RETURN
    ENDIF

    IF (str_cmp(TRIM(str_in), 'none')) THEN
      coll_pairs = -1.0_num
      RETURN
    ENDIF

    CALL get_token(str_in, tstr1, species1, errcode)
    IF (errcode /= 0) RETURN

    sp1 = as_integer(species1, errcode)
    IF (errcode /= 0) RETURN

    CALL get_token(tstr1, tstr2, species2, errcode)
    IF (errcode /= 0) RETURN

    sp2 = as_integer(species2, errcode)
    IF (errcode /= 0) RETURN

    collstate = 1.0_num
    IF (str_cmp(TRIM(tstr2), 'on') .OR. str_cmp(TRIM(tstr2), '')) THEN
      collstate = 1.0_num
    ELSEIF (str_cmp(TRIM(tstr2), 'off')) THEN
      collstate = -1.0_num
    ELSE
      collstate = as_real(tstr2, errcode)
      IF (errcode /= 0) RETURN
    ENDIF

    IF (coll_pairs_touched(sp1, sp2) .AND. rank == 0) THEN
      DO iu = 1, nio_units ! Print to stdout and to file
        io = io_units(iu)
        WRITE(io,*)
        WRITE(io,*) '*** WARNING ***'
        WRITE(io,*) 'The collide parameter for ' // TRIM(species1) // ' <-> ' &
            // TRIM(species2)
        WRITE(io,*) 'has been set multiple times!'
        WRITE(io,*) 'Collisions will only be carried out once per species pair.'
        WRITE(io,*) 'Later specifications will always override earlier ones.'
        WRITE(io,*)
      ENDDO
    ENDIF

    coll_pairs(sp1, sp2) = collstate
    coll_pairs(sp2, sp1) = collstate
    coll_pairs_touched(sp1, sp2) = .TRUE.
    coll_pairs_touched(sp2, sp1) = .TRUE.

  END SUBROUTINE set_collision_matrix

END MODULE deck_collision_block
