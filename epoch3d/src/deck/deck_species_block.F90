! Copyright (C) 2010-2016 Keith Bennett <K.Bennett@warwick.ac.uk>
! Copyright (C) 2009-2012 Chris Brady <C.S.Brady@warwick.ac.uk>
! Copyright (C) 2012      Martin Ramsay <M.G.Ramsay@warwick.ac.uk>
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

MODULE deck_species_block

  USE strings_advanced
  USE setup
  USE simple_io
  USE utilities
  USE partlist

  IMPLICIT NONE
  SAVE

  PRIVATE
  PUBLIC :: species_deck_initialise, species_deck_finalise
  PUBLIC :: species_block_start, species_block_end
  PUBLIC :: species_block_handle_element, species_block_check

  INTEGER :: species_id, current_block
  INTEGER(KIND=MPI_OFFSET_KIND) :: offset = 0
  CHARACTER(LEN=string_length), DIMENSION(:), POINTER :: species_names
  INTEGER, DIMENSION(:), POINTER :: species_blocks
  LOGICAL :: got_name
  INTEGER :: check_block = c_err_none
  LOGICAL, DIMENSION(:), ALLOCATABLE :: species_charge_set
  INTEGER :: n_secondary_species_in_block
  CHARACTER(LEN=string_length) :: release_species_list
  CHARACTER(LEN=string_length), DIMENSION(:), POINTER :: release_species
  REAL(num), DIMENSION(:), POINTER :: species_ionisation_energies
  REAL(num), DIMENSION(:), POINTER :: ionisation_energies, ionise_to_species
  REAL(num), DIMENSION(:), POINTER :: mass, charge, angular, part_count
  REAL(num), DIMENSION(:), POINTER :: principle, dumpmask_array
  REAL(num) :: species_mass, species_charge, species_dumpmask

CONTAINS

  SUBROUTINE species_deck_initialise

    current_block = 0
    IF (deck_state == c_ds_first) THEN
      n_species = 0
      ALLOCATE(species_names(4))
      ALLOCATE(species_blocks(4))
      ! All the following information is required during c_ds_first so that the
      ! derived ionisation species can be correctly set up
      ALLOCATE(ionise_to_species(4))
      ALLOCATE(release_species(4))
      ALLOCATE(ionisation_energies(4))
      ALLOCATE(mass(4))
      ALLOCATE(charge(4))
      ALLOCATE(principle(4))
      ALLOCATE(angular(4))
      ALLOCATE(part_count(4))
      ALLOCATE(dumpmask_array(4))
      release_species = ''
    ENDIF

  END SUBROUTINE species_deck_initialise



  SUBROUTINE species_deck_finalise

    INTEGER :: i, j, io, iu, nlevels, nrelease
    CHARACTER(LEN=8) :: string
    INTEGER :: errcode
    TYPE(primitive_stack) :: stack

    IF (deck_state == c_ds_first) THEN
      CALL setup_species
      ALLOCATE(species_charge_set(n_species))
      species_charge_set = .FALSE.

      DO i = 1, n_species
        species_list(i)%name = species_names(i)
        IF (rank == 0) THEN
          CALL integer_as_string(i, string)
          PRINT*, 'Name of species ', TRIM(ADJUSTL(string)), ' is ', &
              TRIM(species_names(i))
        ENDIF
        ! This would usually be set after c_ds_first but all of this is required
        ! during setup of derived ionisation species
        species_list(i)%ionise_to_species = INT(ionise_to_species(i))
        species_list(i)%ionisation_energy = ionisation_energies(i)
        species_list(i)%n = INT(principle(i))
        species_list(i)%l = INT(angular(i))
        species_list(i)%mass = mass(i)
        species_list(i)%charge = charge(i)
        species_list(i)%count = INT(part_count(i),i8)
        species_list(i)%dumpmask = INT(dumpmask_array(i))
        IF (species_list(i)%ionise_to_species > 0) &
            species_list(i)%ionise = .TRUE.
      ENDDO

      DEALLOCATE(dumpmask_array)
      DEALLOCATE(part_count)
      DEALLOCATE(principle)
      DEALLOCATE(angular)
      DEALLOCATE(charge)
      DEALLOCATE(mass)
      DEALLOCATE(ionisation_energies)

      DO i = 1, n_species
        IF (TRIM(release_species(i)) == '') CYCLE

        CALL initialise_stack(stack)
        CALL tokenize(release_species(i), stack, errcode)
        nlevels = 0
        j = i
        ! Count number of ionisation levels of species i
        DO WHILE(species_list(j)%ionise)
          nlevels = nlevels + 1
          j = species_list(j)%ionise_to_species
        ENDDO

        ! Count number of release species listed for species i; we need to do
        ! this because sometimes extra values are returned on the stack
        nrelease = 0
        DO j = 1, SIZE(stack%entries)
          IF (stack%entries(j)%value > 0 &
              .AND. stack%entries(j)%value <= n_species) &
                  nrelease = nrelease + 1
        ENDDO

        ! If there's only one release species use it for all ionisation levels
        IF (SIZE(stack%entries) == 1) THEN
          j = i
          species_list(stack%entries(1)%value)%electron = .TRUE.
          DO WHILE(species_list(j)%ionise)
            species_list(j)%release_species = stack%entries(1)%value
            j = species_list(j)%ionise_to_species
          ENDDO
        ! If there's a list of release species use it
        ELSEIF (nlevels == nrelease) THEN
          nlevels = 1
          j = i
          DO WHILE(species_list(j)%ionise)
            species_list(j)%release_species = stack%entries(nlevels)%value
            species_list(stack%entries(nlevels)%value)%electron = .TRUE.
            nlevels = nlevels + 1
            j = species_list(j)%ionise_to_species
          ENDDO
        ! If there's too many or not enough release species specified use the
        ! first one only and throw an error
        ELSE
          j = i
          species_list(stack%entries(1)%value)%electron = .TRUE.
          DO WHILE(species_list(j)%ionise)
            species_list(j)%release_species = stack%entries(1)%value
            j = species_list(j)%ionise_to_species
          ENDDO
          IF (rank == 0) THEN
            DO iu = 1, nio_units ! Print to stdout and to file
              io = io_units(iu)
              WRITE(io,*) '*** WARNING ***'
              WRITE(io,*) 'Incorrect number of release species specified ', &
                  'for ', TRIM(species_names(i)), '. Using only first ', &
                  'specified.'
            ENDDO
          ENDIF
        ENDIF

        CALL deallocate_stack(stack)
      ENDDO
      DEALLOCATE(release_species)
      DEALLOCATE(ionise_to_species)
      DEALLOCATE(species_names)
    ELSE
      DEALLOCATE(species_charge_set)
      DEALLOCATE(species_blocks)

      ! Sanity check
      DO i = 1, n_species
        IF (species_list(i)%species_type == c_species_id_photon) CYCLE
        IF (species_list(i)%mass > c_tiny) CYCLE
        IF (rank == 0) THEN
          DO iu = 1, nio_units ! Print to stdout and to file
            io = io_units(iu)
            WRITE(io,*) '*** ERROR ***'
            WRITE(io,*) 'The species named "' // TRIM(species_list(i)%name) &
                // '" must have a positive mass.'
          ENDDO
        ENDIF
        CALL abort_code(c_err_bad_value)
      ENDDO

      IF (track_ejected_particles) THEN
        ALLOCATE(ejected_list(n_species))
        ejected_list = species_list
        DO i = 1, n_species
          ejected_list(i)%count = 0
          ejected_list(i)%name = &
              'ejected_' // TRIM(ADJUSTL(species_list(i)%name))
          CALL create_empty_partlist(ejected_list(i)%attached_list)
        ENDDO
      ENDIF
    ENDIF

    IF (use_field_ionisation) need_random_state = .TRUE.

  END SUBROUTINE species_deck_finalise



  SUBROUTINE species_block_start

    n_secondary_species_in_block = 0
    current_block = current_block + 1
    got_name = .FALSE.
    species_dumpmask = c_io_always
    IF (deck_state == c_ds_first) RETURN
    species_id = species_blocks(current_block)
    offset = 0

  END SUBROUTINE species_block_start



  SUBROUTINE species_block_end

    CHARACTER(LEN=8) :: id_string
    CHARACTER(LEN=string_length) :: name
    INTEGER :: i, io, iu, block_species_id

    IF (.NOT.got_name) THEN
      IF (rank == 0) THEN
        CALL integer_as_string(current_block, id_string)
        DO iu = 1, nio_units ! Print to stdout and to file
          io = io_units(iu)
          WRITE(io,*) '*** ERROR ***'
          WRITE(io,*) 'Species block number ', TRIM(id_string), &
              ' has no "name" element.'
        ENDDO
      ENDIF

      check_block = c_err_missing_elements
    ENDIF

    IF (deck_state == c_ds_first) THEN
      block_species_id = n_species
      charge(n_species) = species_charge
      mass(n_species) = species_mass
      IF (n_secondary_species_in_block > 0) THEN
        ! Create an empty species for each ionisation energy listed in species
        ! block
        release_species(n_species) = release_species_list
        DO i = 1, n_secondary_species_in_block
          CALL integer_as_string(i, id_string)
          name = TRIM(TRIM(species_names(block_species_id))//id_string)
          CALL create_ionisation_species_from_name(name, &
              species_ionisation_energies(i), &
              n_secondary_species_in_block + 1 - i)
        ENDDO
        DEALLOCATE(species_ionisation_energies)
      ENDIF
    ENDIF

  END SUBROUTINE species_block_end



  FUNCTION species_block_handle_element(element, value) RESULT(errcode)

    CHARACTER(*), INTENT(IN) :: element, value
    INTEGER :: errcode
    TYPE(primitive_stack) :: stack
    REAL(num) :: dmin, mult
    CHARACTER(LEN=string_length) :: filename, mult_string
    LOGICAL :: got_file, dump
    INTEGER :: i, j, io, iu, n

    errcode = c_err_none
    IF (value == blank .OR. element == blank) RETURN

    IF (str_cmp(element, 'name')) THEN
      IF (got_name) THEN
        errcode = c_err_preset_element
        RETURN
      ENDIF
      got_name = .TRUE.
      IF (deck_state /= c_ds_first) RETURN
      CALL grow_array(species_blocks, current_block)
      species_blocks(current_block) = create_species_number_from_name(value)
      RETURN
    ENDIF

    ! Collect ionisation energies for the species
    IF (str_cmp(element, 'ionisation_energies')) THEN
      IF (deck_state == c_ds_first) THEN
        NULLIFY(species_ionisation_energies)
        CALL initialise_stack(stack)
        CALL tokenize(value, stack, errcode)
        CALL evaluate_and_return_all(stack, 0, 0, 0, &
            n_secondary_species_in_block, species_ionisation_energies, errcode)
        CALL deallocate_stack(stack)
      ENDIF
      RETURN
    ENDIF

    IF (str_cmp(element, 'ionisation_electron_species') &
        .OR. str_cmp(element, 'electron_species') &
        .OR. str_cmp(element, 'electron')) THEN
      release_species_list = value
      RETURN
    ENDIF

    IF (str_cmp(element, 'mass')) THEN
      species_mass = as_real_print(value, element, errcode) * m0
    ENDIF

    IF (str_cmp(element, 'charge')) THEN
      species_charge = as_real_print(value, element, errcode) * q0
    ENDIF

    IF (str_cmp(element, 'dump')) THEN
      dump = as_logical_print(value, element, errcode)
      IF (dump) THEN
        species_dumpmask = c_io_always
      ELSE
        species_dumpmask = c_io_never
      ENDIF
    ENDIF

    IF (str_cmp(element, 'dumpmask')) THEN
      species_dumpmask = as_integer_print(value, element, errcode)
    ENDIF

    IF (deck_state == c_ds_first) RETURN

    ! *************************************************************
    ! This section identifies a species. Generic
    ! but currently only used in photon production
    ! *************************************************************
    IF (str_cmp(element, 'identify')) THEN
      CALL identify_species(value, errcode)
      RETURN
    ENDIF

    IF (str_cmp(element, 'mass')) THEN
      species_list(species_id)%mass = species_mass
      ! Find the release species for each ionising species and subtract the
      ! release mass from the ionising species and each child species. Doing it
      ! like this ensures the right number of electron masses is removed for
      ! each ion.
      DO i = 1, n_species
        IF (species_id == species_list(i)%release_species) THEN
          j = species_list(i)%ionise_to_species
          DO WHILE(j > 0)
            species_list(j)%mass = species_list(j)%mass &
                - species_list(species_id)%mass
            j = species_list(j)%ionise_to_species
          ENDDO
        ENDIF
      ENDDO
      IF (species_list(species_id)%mass < 0) THEN
        IF (rank == 0) THEN
          DO iu = 1, nio_units ! Print to stdout and to file
            io = io_units(iu)
            WRITE(io,*) '*** ERROR ***'
            WRITE(io,*) 'Particle species cannot have negative mass.'
          ENDDO
        ENDIF
        errcode = c_err_bad_value
      ENDIF
      RETURN
    ENDIF

    IF (str_cmp(element, 'charge')) THEN
      species_list(species_id)%charge = species_charge
      species_charge_set(species_id) = .TRUE.
      ! Find the release species for each ionising species and subtract the
      ! release charge from the ionising species and each child species. Doing
      ! it like this ensures the right number of electron charges is removed for
      ! each ion. The species charge is considered set for the derived ionised
      ! species if it is touched in this routine.
      DO i = 1, n_species
        IF (species_id == species_list(i)%release_species) THEN
          j = species_list(i)%ionise_to_species
          DO WHILE(j > 0)
            species_list(j)%charge = species_list(j)%charge &
                - species_list(species_id)%charge
            species_charge_set(j) = .TRUE.
            j = species_list(j)%ionise_to_species
          ENDDO
        ENDIF
      ENDDO
      RETURN
    ENDIF

    IF (str_cmp(element, 'frac') .OR. str_cmp(element, 'fraction')) THEN
      IF (npart_global >= 0) THEN
        species_list(species_id)%count = &
            INT(as_real_print(value, element, errcode) * npart_global, i8)
      ELSE
        species_list(species_id)%count = 0
      ENDIF
      RETURN
    ENDIF

    IF (str_cmp(element, 'npart')) THEN
      species_list(species_id)%count = &
          as_long_integer_print(value, element, errcode)
      RETURN
    ENDIF

    IF (str_cmp(element, 'npart_per_cell')) THEN
      species_list(species_id)%npart_per_cell = &
          as_real_print(value, element, errcode)
      RETURN
    ENDIF

    IF (str_cmp(element, 'dump') .OR. str_cmp(element, 'dumpmask')) THEN
      species_list(species_id)%dumpmask = INT(species_dumpmask)
      RETURN
    ENDIF

    IF (str_cmp(element, 'immobile')) THEN
      species_list(species_id)%immobile = &
          as_logical_print(value, element, errcode)
      RETURN
    ENDIF

    ! *************************************************************
    ! This section sets properties for tracer particles
    ! *************************************************************
    IF (str_cmp(element, 'tracer')) THEN
#ifndef NO_TRACER_PARTICLES
      species_list(species_id)%tracer = &
          as_logical_print(value, element, errcode)
#else
      IF (as_logical_print(value, element, errcode)) THEN
        errcode = c_err_pp_options_wrong
        extended_error_string = '-DNO_TRACER_PARTICLES'
      ENDIF
#endif
      RETURN
    ENDIF

    ! *************************************************************
    ! This section sets properties for particle splitting
    ! *************************************************************
    IF (str_cmp(element, 'split')) THEN
      species_list(species_id)%split = as_logical_print(value, element, errcode)
      RETURN
    ENDIF

    IF (str_cmp(element, 'npart_max')) THEN
      species_list(species_id)%npart_max = &
          as_long_integer_print(value, element, errcode)
      RETURN
    ENDIF

    ! *************************************************************
    ! This section sets properties for migration
    ! *************************************************************

    IF (str_cmp(element, 'migrate')) THEN
      species_list(species_id)%migrate%this_species = &
          as_logical_print(value, element, errcode)
      RETURN
    ENDIF

    IF (str_cmp(element, 'promote_to')) THEN
      species_list(species_id)%migrate%promote_to_species = &
          as_integer_print(value, element, errcode)
      RETURN
    ENDIF

    IF (str_cmp(element, 'demote_to')) THEN
      species_list(species_id)%migrate%demote_to_species = &
          as_integer_print(value, element, errcode)
      RETURN
    ENDIF

    IF (str_cmp(element, 'promote_multiplier')) THEN
      species_list(species_id)%migrate%promotion_energy_factor = &
          as_real_print(value, element, errcode)
      RETURN
    ENDIF

    IF (str_cmp(element, 'demote_multiplier')) THEN
      species_list(species_id)%migrate%demotion_energy_factor = &
          as_real_print(value, element, errcode)
      RETURN
    ENDIF

    IF (str_cmp(element, 'promote_density')) THEN
      species_list(species_id)%migrate%promotion_density = &
          as_real_print(value, element, errcode)
      RETURN
    ENDIF

    IF (str_cmp(element, 'demote_density')) THEN
      species_list(species_id)%migrate%demotion_density = &
          as_real_print(value, element, errcode)
      RETURN
    ENDIF

    ! Initial conditions

    IF (str_cmp(element, 'offset')) THEN
      offset = as_long_integer_print(value, element, errcode)
      RETURN
    ENDIF

    IF (str_cmp(element, 'density_min') .OR. str_cmp(element, 'minrho')) THEN
      dmin = as_real_print(value, element, errcode)
      IF (dmin <= 0.0_num) dmin = EPSILON(1.0_num)
      initial_conditions(species_id)%density_min = dmin
      RETURN
    ENDIF

    IF (str_cmp(element, 'density_max') .OR. str_cmp(element, 'maxrho')) THEN
      initial_conditions(species_id)%density_max = &
          as_real_print(value, element, errcode)
      RETURN
    ENDIF

    CALL get_filename(value, filename, got_file, errcode)

    mult = 1.0_num

    IF (str_cmp(element, 'density') .OR. str_cmp(element, 'rho') &
        .OR. str_cmp(element, 'mass_density')) THEN

      IF (str_cmp(element, 'mass_density')) THEN
        mult = 1.0_num / species_list(species_id)%mass
        mult_string = '/ species_list(species_id)%mass'
      ENDIF

      CALL fill_array(species_list(species_id)%density_function, &
          initial_conditions(species_id)%density, &
          mult, mult_string, element, value, filename, got_file)
      RETURN
    ENDIF

    IF (str_cmp(element, 'drift_x')) THEN
      n = 1
      CALL fill_array(species_list(species_id)%drift_function(n), &
          initial_conditions(species_id)%drift(:,:,:,n), &
          mult, mult_string, element, value, filename, got_file)
      RETURN
    ENDIF

    IF (str_cmp(element, 'drift_y')) THEN
      n = 2
      CALL fill_array(species_list(species_id)%drift_function(n), &
          initial_conditions(species_id)%drift(:,:,:,n), &
          mult, mult_string, element, value, filename, got_file)
      RETURN
    ENDIF

    IF (str_cmp(element, 'drift_z')) THEN
      n = 3
      CALL fill_array(species_list(species_id)%drift_function(n), &
          initial_conditions(species_id)%drift(:,:,:,n), &
          mult, mult_string, element, value, filename, got_file)
      RETURN
    ENDIF

    mult_string = '* ev / kb'

    IF (str_cmp(element, 'temp') .OR. str_cmp(element, 'temp_k') &
        .OR. str_cmp(element, 'temp_ev')) THEN
      IF (str_cmp(element, 'temp_ev')) mult = ev / kb

      n = 1
      CALL fill_array(species_list(species_id)%temperature_function(n), &
          initial_conditions(species_id)%temp(:,:,:,n), &
          mult, mult_string, element, value, filename, got_file)
      n = 2
      CALL fill_array(species_list(species_id)%temperature_function(n), &
          initial_conditions(species_id)%temp(:,:,:,n), &
          mult, mult_string, element, value, filename, got_file)
      n = 3
      CALL fill_array(species_list(species_id)%temperature_function(n), &
          initial_conditions(species_id)%temp(:,:,:,n), &
          mult, mult_string, element, value, filename, got_file)

      debug_mode = .FALSE.
      initial_conditions(species_id)%temp(:,:,:,2) = &
          initial_conditions(species_id)%temp(:,:,:,n)
      initial_conditions(species_id)%temp(:,:,:,3) = &
          initial_conditions(species_id)%temp(:,:,:,n)
      RETURN
    ENDIF

    IF (str_cmp(element, 'temp_x') .OR. str_cmp(element, 'temp_x_k') &
        .OR. str_cmp(element, 'temp_x_ev')) THEN
      IF (str_cmp(element, 'temp_x_ev')) mult = ev / kb

      n = 1
      CALL fill_array(species_list(species_id)%temperature_function(n), &
          initial_conditions(species_id)%temp(:,:,:,n), &
          mult, mult_string, element, value, filename, got_file)
      RETURN
    ENDIF

    IF (str_cmp(element, 'temp_y') .OR. str_cmp(element, 'temp_y_k') &
        .OR. str_cmp(element, 'temp_y_ev')) THEN
      IF (str_cmp(element, 'temp_y_ev')) mult = ev / kb

      n = 2
      CALL fill_array(species_list(species_id)%temperature_function(n), &
          initial_conditions(species_id)%temp(:,:,:,n), &
          mult, mult_string, element, value, filename, got_file)
      RETURN
    ENDIF

    IF (str_cmp(element, 'temp_z') .OR. str_cmp(element, 'temp_z_k') &
        .OR. str_cmp(element, 'temp_z_ev')) THEN
      IF (str_cmp(element, 'temp_z_ev')) mult = ev / kb

      n = 3
      CALL fill_array(species_list(species_id)%temperature_function(n), &
          initial_conditions(species_id)%temp(:,:,:,n), &
          mult, mult_string, element, value, filename, got_file)
      RETURN
    ENDIF

    errcode = c_err_unknown_element

  END FUNCTION species_block_handle_element



  FUNCTION species_block_check() RESULT(errcode)

    INTEGER :: errcode
    INTEGER :: i, io, iu

    errcode = check_block

    IF (deck_state == c_ds_first) RETURN

    DO i = 1, n_species
      IF (species_list(i)%mass < 0) THEN
        IF (rank == 0) THEN
          DO iu = 1, nio_units ! Print to stdout and to file
            io = io_units(iu)
            WRITE(io,*) '*** ERROR ***'
            WRITE(io,*) 'No mass specified for particle species "', &
                TRIM(species_list(i)%name), '"'
          ENDDO
        ENDIF
        errcode = c_err_missing_elements
      ENDIF
      IF (.NOT. species_charge_set(i)) THEN
        IF (rank == 0) THEN
          DO iu = 1, nio_units ! Print to stdout and to file
            io = io_units(iu)
            WRITE(io,*) '*** ERROR ***'
            WRITE(io,*) 'No charge specified for particle species "', &
                TRIM(species_list(i)%name), '"'
          ENDDO
        ENDIF
        errcode = c_err_missing_elements
      ENDIF
      IF (species_list(i)%npart_per_cell >= 0) THEN
        IF (species_list(i)%count >= 0 .AND. rank == 0) THEN
          DO iu = 1, nio_units ! Print to stdout and to file
            io = io_units(iu)
            WRITE(io,*) '*** WARNING ***'
            WRITE(io,*) 'Two forms of npart used for particle species "', &
                TRIM(species_list(i)%name), '"'
            WRITE(io,*) 'Just using "npart_per_cell".'
          ENDDO
        ENDIF
        species_list(i)%count = INT(species_list(i)%npart_per_cell, i8)
      ENDIF
    ENDDO

  END FUNCTION species_block_check



  FUNCTION create_species_number_from_name(name)

    CHARACTER(*), INTENT(IN) :: name
    INTEGER :: create_species_number_from_name
    INTEGER :: i, io, iu
    TYPE(stack_element) :: block

    DO i = 1, n_species
      IF (str_cmp(name, species_names(i))) THEN
        create_species_number_from_name = i
        RETURN
      ENDIF
    ENDDO

    ! If we're here then then named species doesn't yet exist

    ! First issue a warning message if the name overrides a built-in one
    CALL load_block(name, block)
    IF (block%ptype /= c_pt_bad .AND. block%ptype /= c_pt_null) THEN
      IF (rank == 0) THEN
        DO iu = 1, nio_units ! Print to stdout and to file
          io = io_units(iu)
          WRITE(io,*) '*** ERROR ***'
          WRITE(io,*) 'The species name "' // TRIM(name) // '" is not valid.'
          WRITE(io,*) 'Please choose a different name and try again.'
        ENDDO
      ENDIF
      CALL abort_code(c_err_bad_value)
    ENDIF

    n_species = n_species + 1
    create_species_number_from_name = n_species

    CALL grow_array(species_names, n_species)
    CALL grow_array(ionise_to_species, n_species)
    CALL grow_array(release_species, n_species)
    CALL grow_array(mass, n_species)
    CALL grow_array(charge, n_species)
    CALL grow_array(ionisation_energies, n_species)
    CALL grow_array(principle, n_species)
    CALL grow_array(angular, n_species)
    CALL grow_array(part_count, n_species)
    CALL grow_array(dumpmask_array, n_species)

    species_names(n_species) = TRIM(name)
    ionise_to_species(n_species) = -1
    release_species(n_species) = ''
    mass(n_species) = -1.0_num
    charge(n_species) = 0.0_num
    ionisation_energies(n_species) = HUGE(0.0_num)
    principle(n_species) = -1
    angular(n_species) = -1
    part_count(n_species) = -1
    dumpmask_array(n_species) = species_dumpmask

    RETURN

  END FUNCTION create_species_number_from_name



  SUBROUTINE create_ionisation_species_from_name(name, ionisation_energy, &
      n_electrons)

    CHARACTER(*), INTENT(IN) :: name
    REAL(num), INTENT(IN) :: ionisation_energy
    INTEGER, INTENT(IN) :: n_electrons
    INTEGER :: i, n, l

    DO i = 1, n_species
      IF (str_cmp(name, species_names(i))) RETURN
    ENDDO
    ! This calculates the principle and angular quantum number based on the
    ! assumption that shells are filled as they would be in the ground state
    ! e.g. 1s, 2s, 2p, 3s, 3p, 4s, 3d, 4p, 5s, 4d, 5p, 6s, 4f, 5d, 6p, 7s, etc
    n = 0
    l = 0
    i = 0
    DO WHILE(n_electrons > i)
      n = n + 1
      DO l = (n - 1) / 2, 0, -1
        i = i + 4 * l + 2
        IF (n_electrons <= i) THEN
          n = n - l
          EXIT
        ENDIF
      ENDDO
    ENDDO
    principle(n_species) = n
    angular(n_species) = l
    ionisation_energies(n_species) = ionisation_energy
    ionise_to_species(n_species) = n_species + 1
    n_species = n_species + 1
    CALL grow_array(species_names, n_species)
    species_names(n_species) = TRIM(name)
    CALL grow_array(ionise_to_species, n_species)
    ionise_to_species(n_species) = -1
    CALL grow_array(release_species, n_species)
    release_species(n_species) = ''
    CALL grow_array(mass, n_species)
    mass(n_species) = species_mass
    CALL grow_array(charge, n_species)
    charge(n_species) = species_charge
    CALL grow_array(ionisation_energies, n_species)
    ionisation_energies(n_species) = HUGE(0.0_num)
    CALL grow_array(principle, n_species)
    principle(n_species) = -1
    CALL grow_array(angular, n_species)
    angular(n_species) = -1
    CALL grow_array(part_count, n_species)
    part_count(n_species) = 0
    CALL grow_array(dumpmask_array, n_species)
    dumpmask_array(n_species) = species_dumpmask
    RETURN

  END SUBROUTINE create_ionisation_species_from_name



  FUNCTION species_number_from_name(name)

    CHARACTER(*), INTENT(IN) :: name
    INTEGER :: species_number_from_name
    INTEGER :: i

    DO i = 1, n_species
      IF (str_cmp(name, species_list(i)%name)) THEN
        species_number_from_name = i
        RETURN
      ENDIF
    ENDDO
    species_number_from_name = -1
    RETURN

  END FUNCTION species_number_from_name



  SUBROUTINE fill_array(output, array, mult, mult_string, element, value, &
      filename, got_file)

    TYPE(primitive_stack), INTENT(INOUT) :: output
    REAL(num), DIMENSION(-2:,-2:,-2:), INTENT(INOUT) :: array
    REAL(num), INTENT(IN) :: mult
    CHARACTER(LEN=*), INTENT(IN) :: mult_string, element, value, filename
    LOGICAL, INTENT(IN) :: got_file
    TYPE(stack_element) :: block
    TYPE(primitive_stack) :: stack
    INTEGER :: io, iu, ix, iy, iz

    CALL initialise_stack(stack)
    IF (got_file) THEN
      IF (move_window) THEN
        IF (rank == 0) THEN
          DO iu = 1, nio_units ! Print to stdout and to file
            io = io_units(iu)
            WRITE(io,*) '*** ERROR ***'
            WRITE(io,*) 'Cannot load from file whilst using a moving window.'
          ENDDO
        ENDIF
        errcode = c_err_bad_value
        RETURN
      ENDIF

      CALL load_single_array_from_file(filename, array, offset, errcode)

      CALL load_block(species_list(species_id)%name, block)
      CALL push_to_stack(stack, block)
      CALL load_block(element, block)
      CALL push_to_stack(stack, block)
      IF (ABS(mult - 1.0_num) > c_tiny) array = mult * array
    ELSE
      CALL tokenize(value, stack, errcode)
      IF (ABS(mult - 1.0_num) > c_tiny) &
          CALL tokenize(mult_string, stack, errcode)

      ! Sanity check
      array(1,1,1) = evaluate_at_point(stack, 1, 1, 1, errcode)
      IF (errcode /= c_err_none) THEN
        IF (rank == 0) THEN
          DO iu = 1, nio_units ! Print to stdout and to file
            io = io_units(iu)
            WRITE(io,*) '*** ERROR ***'
            WRITE(io,*) 'Unable to parse input deck.'
          ENDDO
        ENDIF
        CALL abort_code(errcode)
      ENDIF

      DO iz = -2, nz+3
      DO iy = -2, ny+3
      DO ix = -2, nx+3
        array(ix,iy,iz) = evaluate_at_point(stack, ix, iy, iz, errcode)
      ENDDO
      ENDDO
      ENDDO
    ENDIF

    CALL deallocate_stack(output)
    output = stack

  END SUBROUTINE fill_array



  SUBROUTINE identify_species(value, errcode)

    CHARACTER(*), INTENT(IN) :: value
    INTEGER, INTENT(INOUT) :: errcode

    ! Just a plain old electron
    IF (str_cmp(value, 'electron')) THEN
      species_list(species_id)%charge = -q0
      species_list(species_id)%mass = m0
      species_charge_set(species_id) = .TRUE.
      species_list(species_id)%species_type = c_species_id_electron
      RETURN
    ENDIF

    IF (str_cmp(value, 'proton')) THEN
      species_list(species_id)%charge = q0
      species_list(species_id)%mass = m0 * 1836.2_num
      species_charge_set(species_id) = .TRUE.
      species_list(species_id)%species_type = c_species_id_proton
      RETURN
    ENDIF

    IF (str_cmp(value, 'positron')) THEN
      species_list(species_id)%charge = q0
      species_list(species_id)%mass = m0
      species_charge_set(species_id) = .TRUE.
      species_list(species_id)%species_type = c_species_id_positron
      RETURN
    ENDIF

#ifndef PHOTONS
    extended_error_string = 'Cannot identify species "' &
        // TRIM(species_list(species_id)%name) // '" as "' // TRIM(value) &
        // '" because compiler option -DPHOTONS has not been set.'
#else
#ifndef TRIDENT_PHOTONS
    extended_error_string = 'Cannot identify species "' &
        // TRIM(species_list(species_id)%name) // '" as "' // TRIM(value) &
        // '" because compiler option -DTRIDENT_PHOTONS has not been set.'
#endif
#endif

    ! trident process electron
    IF (str_cmp(value, 'trident_electron')) THEN
      species_list(species_id)%charge = -q0
      species_list(species_id)%mass = m0
      species_list(species_id)%species_type = c_species_id_electron
      species_charge_set(species_id) = .TRUE.
#if defined(PHOTONS) && defined(TRIDENT_PHOTONS)
      trident_electron_species = species_id
#else
      errcode = c_err_generic_warning
#endif
      RETURN
    ENDIF

    ! Breit Wheeler process electron
    IF (str_cmp(value, 'bw_electron')) THEN
      species_list(species_id)%charge = -q0
      species_list(species_id)%mass = m0
      species_list(species_id)%species_type = c_species_id_electron
      species_charge_set(species_id) = .TRUE.
#ifdef PHOTONS
      breit_wheeler_electron_species = species_id
#else
      errcode = c_err_generic_warning
#endif
      RETURN
    ENDIF

    IF (str_cmp(value, 'bw_positron')) THEN
      species_list(species_id)%charge = q0
      species_list(species_id)%mass = m0
      species_charge_set(species_id) = .TRUE.
      species_list(species_id)%species_type = c_species_id_positron
#ifdef PHOTONS
      breit_wheeler_positron_species = species_id
#endif
      RETURN
    ENDIF

    ! trident process positron
    IF (str_cmp(value, 'trident_positron')) THEN
      species_list(species_id)%charge = q0
      species_list(species_id)%mass = m0
      species_list(species_id)%species_type = c_species_id_positron
      species_charge_set(species_id) = .TRUE.
#if defined(PHOTONS) && defined(TRIDENT_PHOTONS)
      trident_electron_species = species_id
#else
      errcode = c_err_generic_warning
#endif
      RETURN
    ENDIF

    IF (str_cmp(value, 'photon')) THEN
      species_list(species_id)%charge = 0.0_num
      species_list(species_id)%mass = 0.0_num
      species_list(species_id)%species_type = c_species_id_photon
      species_charge_set(species_id) = .TRUE.
#ifdef PHOTONS
      IF (photon_species == -1) photon_species = species_id
#else
      errcode = c_err_generic_warning
#endif
      RETURN
    ENDIF

    errcode = IOR(errcode, c_err_bad_value)

  END SUBROUTINE identify_species

END MODULE deck_species_block
