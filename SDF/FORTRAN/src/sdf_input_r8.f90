!
! SDF (Self-Describing Format) Fortran Library
! Copyright (c) 2011-2016, SDF Development Team
!
! Distributed under the terms of the BSD 3-clause License.
! See the LICENSE file for details.
!

MODULE sdf_input_r8

  USE sdf_input_ru

  IMPLICIT NONE

CONTAINS

  SUBROUTINE read_header_r8(h, step, time, code_name, code_io_version, &
      string_length, restart_flag, other_domains)

    TYPE(sdf_file_handle) :: h
    INTEGER, INTENT(OUT) :: step
    REAL(r8), INTENT(OUT) :: time
    CHARACTER(LEN=*), INTENT(OUT), OPTIONAL :: code_name
    INTEGER, INTENT(OUT), OPTIONAL :: code_io_version, string_length
    LOGICAL, INTENT(OUT), OPTIONAL :: restart_flag, other_domains

    CALL read_header_ru(h, step)

    time = REAL(h%time,r8)
    IF (PRESENT(code_io_version)) code_io_version = h%code_io_version
    IF (PRESENT(restart_flag)) restart_flag = h%restart_flag
    IF (PRESENT(other_domains)) other_domains = h%other_domains
    IF (PRESENT(code_name)) CALL sdf_safe_copy_string(h%code_name, code_name)
    IF (PRESENT(string_length)) string_length = h%string_length

  END SUBROUTINE read_header_r8



  SUBROUTINE read_constant_real_r8(h, value)

    TYPE(sdf_file_handle) :: h
    REAL(r8), INTENT(OUT) :: value
    REAL(r4) :: real4
    REAL(r8) :: real8
    TYPE(sdf_block_type), POINTER :: b

    CALL sdf_read_constant(h)

    b => h%current_block

    IF (b%datatype == c_datatype_real4) THEN
      real4 = TRANSFER(b%const_value, real4)
      value = REAL(real4,r8)
    ELSE IF (b%datatype == c_datatype_real8) THEN
      real8 = TRANSFER(b%const_value, real8)
      value = REAL(real8,r8)
    ENDIF

  END SUBROUTINE read_constant_real_r8



  SUBROUTINE read_1d_array_real_r8(h, values)

    TYPE(sdf_file_handle) :: h
    REAL(r8), DIMENSION(:), INTENT(OUT) :: values
    INTEGER, DIMENSION(c_maxdims) :: dims
    INTEGER :: errcode, n1
    TYPE(sdf_block_type), POINTER :: b

    IF (sdf_check_block_header(h)) RETURN

    b => h%current_block
    IF (.NOT. b%done_info) CALL sdf_read_array_info(h, dims)

    n1 = b%dims(1)

    h%current_location = b%data_location

    IF (h%rank == h%rank_master) THEN
      CALL MPI_FILE_READ_AT(h%filehandle, h%current_location, values, n1, &
          b%mpitype, MPI_STATUS_IGNORE, errcode)
    ENDIF

    CALL MPI_BCAST(values, n1, b%mpitype, h%rank_master, h%comm, errcode)

    h%current_location = b%next_block_location
    b%done_data = .TRUE.

  END SUBROUTINE read_1d_array_real_r8



  SUBROUTINE read_2d_array_real_r8(h, values)

    TYPE(sdf_file_handle) :: h
    REAL(r8), DIMENSION(:,:), INTENT(OUT) :: values
    INTEGER, DIMENSION(c_maxdims) :: dims
    INTEGER :: errcode, n1
    TYPE(sdf_block_type), POINTER :: b

    IF (sdf_check_block_header(h)) RETURN

    b => h%current_block
    IF (.NOT. b%done_info) CALL sdf_read_array_info(h, dims)

    n1 = b%dims(1) * b%dims(2)

    h%current_location = b%data_location

    IF (h%rank == h%rank_master) THEN
      CALL MPI_FILE_READ_AT(h%filehandle, h%current_location, values, n1, &
          b%mpitype, MPI_STATUS_IGNORE, errcode)
    ENDIF

    CALL MPI_BCAST(values, n1, b%mpitype, h%rank_master, h%comm, errcode)

    h%current_location = b%next_block_location
    b%done_data = .TRUE.

  END SUBROUTINE read_2d_array_real_r8



  SUBROUTINE read_3d_array_real_r8(h, values)

    TYPE(sdf_file_handle) :: h
    REAL(r8), DIMENSION(:,:,:), INTENT(OUT) :: values
    INTEGER, DIMENSION(c_maxdims) :: dims
    INTEGER :: errcode, n1
    TYPE(sdf_block_type), POINTER :: b

    IF (sdf_check_block_header(h)) RETURN

    b => h%current_block
    IF (.NOT. b%done_info) CALL sdf_read_array_info(h, dims)

    n1 = b%dims(1) * b%dims(2) * b%dims(3)

    h%current_location = b%data_location

    IF (h%rank == h%rank_master) THEN
      CALL MPI_FILE_READ_AT(h%filehandle, h%current_location, values, n1, &
          b%mpitype, MPI_STATUS_IGNORE, errcode)
    ENDIF

    CALL MPI_BCAST(values, n1, b%mpitype, h%rank_master, h%comm, errcode)

    h%current_location = b%next_block_location
    b%done_data = .TRUE.

  END SUBROUTINE read_3d_array_real_r8

END MODULE sdf_input_r8
