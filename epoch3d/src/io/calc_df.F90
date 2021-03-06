! Copyright (C) 2010-2015 Keith Bennett <K.Bennett@warwick.ac.uk>
! Copyright (C) 2009-2012 Chris Brady <C.S.Brady@warwick.ac.uk>
! Copyright (C) 2011-2012 Martin Ramsay <M.G.Ramsay@warwick.ac.uk>
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

MODULE calc_df

  USE boundary

  IMPLICIT NONE

CONTAINS

  SUBROUTINE calc_boundary(data_array)

    REAL(num), DIMENSION(-2:,-2:,-2:), INTENT(OUT) :: data_array
    INTEGER :: i, j, k

    CALL processor_summation_bcs(data_array, ng)

    IF (x_min_boundary .AND. bc_particle(c_bd_x_min) /= c_bc_periodic &
        .AND. bc_particle(c_bd_x_min) /= c_bc_reflect) THEN
      DO k = -2, nz+3
      DO j = -2, ny+3
        data_array(1,j,k) = data_array(1,j,k) + data_array( 0,j,k)
        data_array(2,j,k) = data_array(2,j,k) + data_array(-1,j,k)
        data_array(3,j,k) = data_array(3,j,k) + data_array(-2,j,k)
      ENDDO
      ENDDO
    ENDIF
    IF (x_max_boundary .AND. bc_particle(c_bd_x_max) /= c_bc_periodic &
        .AND. bc_particle(c_bd_x_max) /= c_bc_reflect) THEN
      DO k = -2, nz+3
      DO j = -2, ny+3
        data_array(nx-2,j,k) = data_array(nx-2,j,k) + data_array(nx+3,j,k)
        data_array(nx-1,j,k) = data_array(nx-1,j,k) + data_array(nx+2,j,k)
        data_array(nx  ,j,k) = data_array(nx  ,j,k) + data_array(nx+1,j,k)
      ENDDO
      ENDDO
    ENDIF

    IF (y_min_boundary .AND. bc_particle(c_bd_y_min) /= c_bc_periodic &
        .AND. bc_particle(c_bd_y_min) /= c_bc_reflect) THEN
      DO k = -2, nz+3
      DO i = -2, nx+3
        data_array(i,1,k) = data_array(i,1,k) + data_array(i, 0,k)
        data_array(i,2,k) = data_array(i,2,k) + data_array(i,-1,k)
        data_array(i,3,k) = data_array(i,3,k) + data_array(i,-2,k)
      ENDDO
      ENDDO
    ENDIF
    IF (y_max_boundary .AND. bc_particle(c_bd_y_max) /= c_bc_periodic &
        .AND. bc_particle(c_bd_y_max) /= c_bc_reflect) THEN
      DO k = -2, nz+3
      DO i = -2, nx+3
        data_array(i,ny-2,k) = data_array(i,ny-2,k) + data_array(i,ny+3,k)
        data_array(i,ny-1,k) = data_array(i,ny-1,k) + data_array(i,ny+2,k)
        data_array(i,ny  ,k) = data_array(i,ny  ,k) + data_array(i,ny+1,k)
      ENDDO
      ENDDO
    ENDIF

    IF (z_min_boundary .AND. bc_particle(c_bd_z_min) /= c_bc_periodic &
        .AND. bc_particle(c_bd_z_min) /= c_bc_reflect) THEN
      DO j = -2, ny+3
      DO i = -2, nx+3
        data_array(i,j,1) = data_array(i,j,1) + data_array(i,j, 0)
        data_array(i,j,2) = data_array(i,j,2) + data_array(i,j,-1)
        data_array(i,j,3) = data_array(i,j,3) + data_array(i,j,-2)
      ENDDO
      ENDDO
    ENDIF
    IF (z_max_boundary .AND. bc_particle(c_bd_z_max) /= c_bc_periodic &
        .AND. bc_particle(c_bd_z_max) /= c_bc_reflect) THEN
      DO j = -2, ny+3
      DO i = -2, nx+3
        data_array(i,j,nz-2) = data_array(i,j,nz-2) + data_array(i,j,nz+3)
        data_array(i,j,nz-1) = data_array(i,j,nz-1) + data_array(i,j,nz+2)
        data_array(i,j,nz  ) = data_array(i,j,nz  ) + data_array(i,j,nz+1)
      ENDDO
      ENDDO
    ENDIF

  END SUBROUTINE calc_boundary



  SUBROUTINE calc_mass_density(data_array, current_species)

    REAL(num), DIMENSION(-2:,-2:,-2:), INTENT(OUT) :: data_array
    INTEGER, INTENT(IN) :: current_species
    ! Properties of the current particle. Copy out of particle arrays for speed
    REAL(num) :: part_m
    ! The data to be weighted onto the grid
    REAL(num) :: wdata
    REAL(num) :: fac, idx
    INTEGER :: ispecies, ix, iy, iz, spec_start, spec_end
    TYPE(particle), POINTER :: current
    LOGICAL :: spec_sum
#include "particle_head.inc"

    data_array = 0.0_num
    part_m = 0.0_num
    fac = 0.0_num

    idx = 1.0_num / dx / dy / dz

    spec_start = current_species
    spec_end = current_species
    spec_sum = .FALSE.

    IF (current_species <= 0) THEN
      spec_start = 1
      spec_end = n_species
      spec_sum = .TRUE.
    ENDIF

    DO ispecies = spec_start, spec_end
      IF (io_list(ispecies)%species_type == c_species_id_photon) CYCLE
#ifndef NO_TRACER_PARTICLES
      IF (spec_sum .AND. io_list(ispecies)%tracer) CYCLE
#endif
      current => io_list(ispecies)%attached_list%head
      part_m  = io_list(ispecies)%mass
      fac = io_list(ispecies)%weight
      wdata = part_m * fac

      DO WHILE (ASSOCIATED(current))
        ! Copy the particle properties out for speed
#ifdef PER_PARTICLE_CHARGE_MASS
        part_m  = current%mass
#ifndef PER_SPECIES_WEIGHT
        fac = current%weight
#endif
        wdata = part_m * fac
#else
#ifndef PER_SPECIES_WEIGHT
        fac = current%weight
        wdata = part_m * fac
#endif
#endif

#include "particle_to_grid.inc"

        DO iz = sf_min, sf_max
        DO iy = sf_min, sf_max
        DO ix = sf_min, sf_max
          data_array(cell_x+ix, cell_y+iy, cell_z+iz) = &
              data_array(cell_x+ix, cell_y+iy, cell_z+iz) &
              + gx(ix) * gy(iy) * gz(iz) * wdata
        ENDDO
        ENDDO
        ENDDO

        current => current%next
      ENDDO
    ENDDO

    data_array = data_array * idx

    CALL calc_boundary(data_array)
    DO ix = 1, 2*c_ndims
      CALL field_zero_gradient(data_array, c_stagger_centre, ix)
    ENDDO

  END SUBROUTINE calc_mass_density



  SUBROUTINE calc_ekbar(data_array, current_species)

    REAL(num), DIMENSION(-2:,-2:,-2:), INTENT(OUT) :: data_array
    INTEGER, INTENT(IN) :: current_species
    ! Properties of the current particle. Copy out of particle arrays for speed
    REAL(num) :: part_ux, part_uy, part_uz, part_mc
    ! The weight of a particle
    REAL(num) :: part_w
    ! The data to be weighted onto the grid
    REAL(num) :: wdata
    REAL(num) :: fac, gamma
    REAL(num), DIMENSION(:,:,:), ALLOCATABLE :: wt
    INTEGER :: ispecies, ix, iy, iz, spec_start, spec_end
    TYPE(particle), POINTER :: current
    LOGICAL :: spec_sum
#include "particle_head.inc"

    ALLOCATE(wt(-2:nx+3,-2:ny+3,-2:nz+3))
    data_array = 0.0_num
    wt = 0.0_num
    part_mc  = 1.0_num
    part_w = 1.0_num

    spec_start = current_species
    spec_end = current_species
    spec_sum = .FALSE.

    IF (current_species <= 0) THEN
      spec_start = 1
      spec_end = n_species
      spec_sum = .TRUE.
    ENDIF

    DO ispecies = spec_start, spec_end
#ifndef NO_TRACER_PARTICLES
      IF (spec_sum .AND. io_list(ispecies)%tracer) CYCLE
#endif
      current => io_list(ispecies)%attached_list%head
      part_mc = c * io_list(ispecies)%mass
      part_w = io_list(ispecies)%weight
      fac = part_mc * part_w * c

      DO WHILE (ASSOCIATED(current))
        ! Copy the particle properties out for speed
#ifdef PER_PARTICLE_CHARGE_MASS
        part_mc = c * current%mass
#ifndef PER_SPECIES_WEIGHT
        part_w = current%weight
#endif
        fac = part_mc * part_w * c
#else
#ifndef PER_SPECIES_WEIGHT
        part_w = current%weight
        fac = part_mc * part_w * c
#endif
#endif

        IF (io_list(ispecies)%species_type /= c_species_id_photon) THEN
          part_ux = current%part_p(1) / part_mc
          part_uy = current%part_p(2) / part_mc
          part_uz = current%part_p(3) / part_mc
          gamma = SQRT(part_ux**2 + part_uy**2 + part_uz**2 + 1.0_num)
          wdata = (gamma - 1.0_num) * fac
        ELSE
#ifdef PHOTONS
          wdata = current%particle_energy * part_w
#else
          wdata = 0.0_num
#endif
        ENDIF

#include "particle_to_grid.inc"

        DO iz = sf_min, sf_max
        DO iy = sf_min, sf_max
        DO ix = sf_min, sf_max
          data_array(cell_x+ix, cell_y+iy, cell_z+iz) = &
              data_array(cell_x+ix, cell_y+iy, cell_z+iz) &
              + gx(ix) * gy(iy) * gz(iz) * wdata
          wt(cell_x+ix, cell_y+iy, cell_z+iz) = &
              wt(cell_x+ix, cell_y+iy, cell_z+iz) &
              + gx(ix) * gy(iy) * gz(iz) * part_w
        ENDDO
        ENDDO
        ENDDO

        current => current%next
      ENDDO
    ENDDO

    CALL calc_boundary(data_array)
    CALL calc_boundary(wt)

    data_array = data_array / MAX(wt, c_tiny)
    DO ix = 1, 2*c_ndims
      CALL field_zero_gradient(data_array, c_stagger_centre, ix)
    ENDDO

    DEALLOCATE(wt)

  END SUBROUTINE calc_ekbar



  SUBROUTINE calc_ekflux(data_array, current_species, direction)

    REAL(num), DIMENSION(-2:,-2:,-2:), INTENT(OUT) :: data_array
    INTEGER, INTENT(IN) :: current_species, direction
    ! Properties of the current particle. Copy out of particle arrays for speed
    REAL(num) :: part_ux, part_uy, part_uz, part_mc
    ! The weight of a particle
    REAL(num) :: part_w
    ! The data to be weighted onto the grid
    REAL(num) :: wdata = 0.0_num
    REAL(num) :: fac, gamma, ek, part_flux, xfac, yfac, zfac
    REAL(num), DIMENSION(:,:,:), ALLOCATABLE :: wt
    INTEGER :: ispecies, ix, iy, iz, spec_start, spec_end
    TYPE(particle), POINTER :: current
    LOGICAL :: spec_sum
#include "particle_head.inc"

    ALLOCATE(wt(-2:nx+3,-2:ny+3,-2:nz+3))
    data_array = 0.0_num
    wt = 0.0_num
    part_mc  = 1.0_num
    part_w = 1.0_num

    xfac = c * dy * dz
    yfac = c * dx * dz
    zfac = c * dx * dy

    spec_start = current_species
    spec_end = current_species
    spec_sum = .FALSE.

    IF (current_species <= 0) THEN
      spec_start = 1
      spec_end = n_species
      spec_sum = .TRUE.
    ENDIF

    DO ispecies = spec_start, spec_end
#ifndef NO_TRACER_PARTICLES
      IF (spec_sum .AND. io_list(ispecies)%tracer) CYCLE
#endif
      current => io_list(ispecies)%attached_list%head
      part_mc = c * io_list(ispecies)%mass
      part_w = io_list(ispecies)%weight
      fac = part_mc * part_w * c

      DO WHILE (ASSOCIATED(current))
        ! Copy the particle properties out for speed
#ifdef PER_PARTICLE_CHARGE_MASS
        part_mc = c * current%mass
#ifndef PER_SPECIES_WEIGHT
        part_w = current%weight
#endif
        fac = part_mc * part_w * c
#else
#ifndef PER_SPECIES_WEIGHT
        part_w = current%weight
        fac = part_mc * part_w * c
#endif
#endif

        IF (io_list(ispecies)%species_type /= c_species_id_photon) THEN
          part_ux = current%part_p(1) / part_mc
          part_uy = current%part_p(2) / part_mc
          part_uz = current%part_p(3) / part_mc
          gamma = SQRT(part_ux**2 + part_uy**2 + part_uz**2 + 1.0_num)
          ek = (gamma - 1.0_num) * fac
        ELSE
#ifdef PHOTONS
          fac = c / current%particle_energy
          part_ux = current%part_p(1) * fac
          part_uy = current%part_p(2) * fac
          part_uz = current%part_p(3) * fac
          ek = current%particle_energy * part_w
          gamma = 1.0_num
#else
          ek = 0.0_num
          gamma = 1.0_num
#endif
        ENDIF

        SELECT CASE(direction)
        CASE(-c_dir_x)
          ! negative flux in x
          part_flux = xfac * part_ux / gamma
          wdata = -ek * MIN(part_flux, 0.0_num)
        CASE( c_dir_x)
          ! positive flux in x
          part_flux = xfac * part_ux / gamma
          wdata =  ek * MAX(part_flux, 0.0_num)
        CASE(-c_dir_y)
          ! negative flux in y
          part_flux = yfac * part_uy / gamma
          wdata = -ek * MIN(part_flux, 0.0_num)
        CASE( c_dir_y)
          ! positive flux in y
          part_flux = yfac * part_uy / gamma
          wdata =  ek * MAX(part_flux, 0.0_num)
        CASE(-c_dir_z)
          ! negative flux in z
          part_flux = zfac * part_uz / gamma
          wdata = -ek * MIN(part_flux, 0.0_num)
        CASE( c_dir_z)
          ! positive flux in z
          part_flux = zfac * part_uz / gamma
          wdata =  ek * MAX(part_flux, 0.0_num)
        END SELECT

#include "particle_to_grid.inc"

        DO iz = sf_min, sf_max
        DO iy = sf_min, sf_max
        DO ix = sf_min, sf_max
          data_array(cell_x+ix, cell_y+iy, cell_z+iz) = &
              data_array(cell_x+ix, cell_y+iy, cell_z+iz) &
              + gx(ix) * gy(iy) * gz(iz) * wdata
          wt(cell_x+ix, cell_y+iy, cell_z+iz) = &
              wt(cell_x+ix, cell_y+iy, cell_z+iz) &
              + gx(ix) * gy(iy) * gz(iz) * part_w
        ENDDO
        ENDDO
        ENDDO

        current => current%next
      ENDDO
    ENDDO

    CALL calc_boundary(data_array)
    CALL calc_boundary(wt)

    data_array = data_array / MAX(wt, c_tiny)
    DO ix = 1, 2*c_ndims
      CALL field_zero_gradient(data_array, c_stagger_centre, ix)
    ENDDO

    DEALLOCATE(wt)

  END SUBROUTINE calc_ekflux



  SUBROUTINE calc_poynt_flux(data_array, direction)

    REAL(num), DIMENSION(-2:,-2:,-2:), INTENT(OUT) :: data_array
    INTEGER, INTENT(IN) :: direction
    INTEGER :: ix, iy, iz
    REAL(num) :: ex_cc, ey_cc, ez_cc, bx_cc, by_cc, bz_cc

    SELECT CASE(direction)
    CASE(c_dir_x)
      DO iz = 1, nz
      DO iy = 1, ny
      DO ix = 1, nx
        ey_cc = 0.5_num  * (ey(ix  , iy-1, iz  ) + ey(ix, iy  , iz  ))
        ez_cc = 0.5_num  * (ez(ix  , iy  , iz-1) + ez(ix, iy  , iz  ))
        by_cc = 0.25_num * (by(ix-1, iy  , iz-1) + by(ix, iy  , iz-1) &
                         +  by(ix-1, iy  , iz  ) + by(ix, iy  , iz  ))
        bz_cc = 0.25_num * (bz(ix-1, iy-1, iz  ) + bz(ix, iy-1, iz  ) &
                         +  bz(ix-1, iy  , iz  ) + bz(ix, iy  , iz  ))
        data_array(ix,iy,iz) = (ey_cc * bz_cc - ez_cc * by_cc) / mu0
      ENDDO
      ENDDO
      ENDDO
    CASE(c_dir_y)
      DO iz = 1, nz
      DO iy = 1, ny
      DO ix = 1, nx
        ex_cc = 0.5_num  * (ex(ix-1, iy  , iz  ) + ex(ix, iy  , iz  ))
        ez_cc = 0.5_num  * (ez(ix  , iy  , iz-1) + ez(ix, iy  , iz  ))
        bx_cc = 0.25_num * (bx(ix  , iy-1, iz-1) + bx(ix, iy  , iz-1) &
                         +  bx(ix  , iy-1, iz  ) + bx(ix, iy  , iz  ))
        bz_cc = 0.25_num * (bz(ix-1, iy-1, iz  ) + bz(ix, iy-1, iz  ) &
                         +  bz(ix-1, iy  , iz  ) + bz(ix, iy  , iz  ))
        data_array(ix,iy,iz) = (ez_cc * bx_cc - ex_cc * bz_cc) / mu0
      ENDDO
      ENDDO
      ENDDO
    CASE(c_dir_z)
      DO iz = 1, nz
      DO iy = 1, ny
      DO ix = 1, nx
        ex_cc = 0.5_num  * (ex(ix-1, iy  , iz  ) + ex(ix, iy, iz  ))
        ey_cc = 0.5_num  * (ey(ix  , iy-1, iz  ) + ey(ix, iy, iz  ))
        bx_cc = 0.25_num * (bx(ix  , iy-1, iz-1) + bx(ix, iy, iz-1) &
                         +  bx(ix  , iy-1, iz  ) + bx(ix, iy, iz  ))
        by_cc = 0.25_num * (by(ix-1, iy  , iz-1) + by(ix, iy, iz-1) &
                         +  by(ix-1, iy  , iz  ) + by(ix, iy, iz  ))
        data_array(ix,iy,iz) = (ex_cc * by_cc - ey_cc * bx_cc) / mu0
      ENDDO
      ENDDO
      ENDDO
    END SELECT

  END SUBROUTINE calc_poynt_flux



  SUBROUTINE calc_charge_density(data_array, current_species)

    REAL(num), DIMENSION(-2:,-2:,-2:), INTENT(OUT) :: data_array
    INTEGER, INTENT(IN) :: current_species
    ! Properties of the current particle. Copy out of particle arrays for speed
    REAL(num) :: part_q
    ! The data to be weighted onto the grid
    REAL(num) :: wdata
    REAL(num) :: fac, idx
    INTEGER :: ispecies, ix, iy, iz, spec_start, spec_end
    TYPE(particle), POINTER :: current
    LOGICAL :: spec_sum
#include "particle_head.inc"

    data_array = 0.0_num
    part_q = 0.0_num
    fac = 0.0_num

    idx = 1.0_num / dx / dy / dz

    spec_start = current_species
    spec_end = current_species
    spec_sum = .FALSE.

    IF (current_species <= 0) THEN
      spec_start = 1
      spec_end = n_species
      spec_sum = .TRUE.
    ENDIF

    DO ispecies = spec_start, spec_end
      IF (io_list(ispecies)%species_type == c_species_id_photon) CYCLE
#ifndef NO_TRACER_PARTICLES
      IF (spec_sum .AND. io_list(ispecies)%tracer) CYCLE
#endif
      current => io_list(ispecies)%attached_list%head
      part_q  = io_list(ispecies)%charge
      fac = io_list(ispecies)%weight
      wdata = part_q * fac

      DO WHILE (ASSOCIATED(current))
        ! Copy the particle properties out for speed
#ifdef PER_PARTICLE_CHARGE_MASS
        part_q  = current%charge
#ifndef PER_SPECIES_WEIGHT
        fac = current%weight
#endif
        wdata = part_q * fac
#else
#ifndef PER_SPECIES_WEIGHT
        fac = current%weight
        wdata = part_q * fac
#endif
#endif

#include "particle_to_grid.inc"

        DO iz = sf_min, sf_max
        DO iy = sf_min, sf_max
        DO ix = sf_min, sf_max
          data_array(cell_x+ix, cell_y+iy, cell_z+iz) = &
              data_array(cell_x+ix, cell_y+iy, cell_z+iz) &
              + gx(ix) * gy(iy) * gz(iz) * wdata
        ENDDO
        ENDDO
        ENDDO

        current => current%next
      ENDDO
    ENDDO

    data_array = data_array * idx

    CALL calc_boundary(data_array)
    DO ix = 1, 2*c_ndims
      CALL field_zero_gradient(data_array, c_stagger_centre, ix)
    ENDDO

  END SUBROUTINE calc_charge_density



  SUBROUTINE calc_number_density(data_array, current_species)

    REAL(num), DIMENSION(-2:,-2:,-2:), INTENT(OUT) :: data_array
    INTEGER, INTENT(IN) :: current_species
    ! The data to be weighted onto the grid
    REAL(num) :: wdata
    REAL(num) :: idx
    INTEGER :: ispecies, ix, iy, iz, spec_start, spec_end
    TYPE(particle), POINTER :: current
    LOGICAL :: spec_sum
#include "particle_head.inc"

    data_array = 0.0_num

    idx = 1.0_num / dx / dy / dz

    spec_start = current_species
    spec_end = current_species
    spec_sum = .FALSE.

    IF (current_species <= 0) THEN
      spec_start = 1
      spec_end = n_species
      spec_sum = .TRUE.
    ENDIF

    DO ispecies = spec_start, spec_end
#ifndef NO_TRACER_PARTICLES
      IF (spec_sum .AND. io_list(ispecies)%tracer) CYCLE
#endif
      current => io_list(ispecies)%attached_list%head
      wdata = io_list(ispecies)%weight

      DO WHILE (ASSOCIATED(current))
#ifndef PER_SPECIES_WEIGHT
        wdata = current%weight
#endif

#include "particle_to_grid.inc"

        DO iz = sf_min, sf_max
        DO iy = sf_min, sf_max
        DO ix = sf_min, sf_max
          data_array(cell_x+ix, cell_y+iy, cell_z+iz) = &
              data_array(cell_x+ix, cell_y+iy, cell_z+iz) &
              + gx(ix) * gy(iy) * gz(iz) * wdata
        ENDDO
        ENDDO
        ENDDO

        current => current%next
      ENDDO
    ENDDO

    data_array = data_array * idx

    CALL calc_boundary(data_array)
    DO ix = 1, 2*c_ndims
      CALL field_zero_gradient(data_array, c_stagger_centre, ix)
    ENDDO

  END SUBROUTINE calc_number_density



  SUBROUTINE calc_temperature(sigma, current_species)

    REAL(num), DIMENSION(-2:,-2:,-2:), INTENT(OUT) :: sigma
    INTEGER, INTENT(IN) :: current_species
    ! Properties of the current particle. Copy out of particle arrays for speed
    REAL(num) :: part_pmx, part_pmy, part_pmz, sqrt_part_m
    ! The weight of a particle
    REAL(num) :: part_w
    REAL(num) :: gf
    REAL(num), DIMENSION(:,:,:), ALLOCATABLE :: part_count, meanx, meany, meanz
    INTEGER :: ispecies, ix, iy, iz, spec_start, spec_end
    TYPE(particle), POINTER :: current
    LOGICAL :: spec_sum
#include "particle_head.inc"

    spec_start = current_species
    spec_end = current_species
    spec_sum = .FALSE.

    IF (current_species <= 0) THEN
      spec_start = 1
      spec_end = n_species
      spec_sum = .TRUE.
    ENDIF

    ALLOCATE(meanx(-2:nx+3,-2:ny+3,-2:nz+3))
    ALLOCATE(meany(-2:nx+3,-2:ny+3,-2:nz+3))
    ALLOCATE(meanz(-2:nx+3,-2:ny+3,-2:nz+3))
    ALLOCATE(part_count(-2:nx+3,-2:ny+3,-2:nz+3))
    meanx = 0.0_num
    meany = 0.0_num
    meanz = 0.0_num
    part_count = 0.0_num
    sigma = 0.0_num

    DO ispecies = spec_start, spec_end
#ifndef NO_TRACER_PARTICLES
      IF (spec_sum .AND. io_list(ispecies)%tracer) CYCLE
#endif
      current => io_list(ispecies)%attached_list%head
      sqrt_part_m  = SQRT(io_list(ispecies)%mass)
      part_w = io_list(ispecies)%weight

      DO WHILE(ASSOCIATED(current))
#ifdef PER_PARTICLE_CHARGE_MASS
        sqrt_part_m  = SQRT(current%mass)
#endif
#ifndef PER_SPECIES_WEIGHT
        part_w = current%weight
#endif
        ! Copy the particle properties out for speed
        part_pmx = current%part_p(1) / sqrt_part_m
        part_pmy = current%part_p(2) / sqrt_part_m
        part_pmz = current%part_p(3) / sqrt_part_m

#include "particle_to_grid.inc"

        DO iz = sf_min, sf_max
        DO iy = sf_min, sf_max
        DO ix = sf_min, sf_max
          gf = gx(ix) * gy(iy) * gz(iz) * part_w
          meanx(cell_x+ix, cell_y+iy, cell_z+iz) = &
              meanx(cell_x+ix, cell_y+iy, cell_z+iz) + gf * part_pmx
          meany(cell_x+ix, cell_y+iy, cell_z+iz) = &
              meany(cell_x+ix, cell_y+iy, cell_z+iz) + gf * part_pmy
          meanz(cell_x+ix, cell_y+iy, cell_z+iz) = &
              meanz(cell_x+ix, cell_y+iy, cell_z+iz) + gf * part_pmz
          part_count(cell_x+ix, cell_y+iy, cell_z+iz) = &
              part_count(cell_x+ix, cell_y+iy, cell_z+iz) + gf
        ENDDO
        ENDDO
        ENDDO
        current => current%next
      ENDDO
    ENDDO

    CALL calc_boundary(meanx)
    CALL calc_boundary(meany)
    CALL calc_boundary(meanz)
    CALL calc_boundary(part_count)

    part_count = MAX(part_count, 1.e-6_num)

    meanx = meanx / part_count
    meany = meany / part_count
    meanz = meanz / part_count

    part_count = 0.0_num
    DO ispecies = spec_start, spec_end
#ifndef NO_TRACER_PARTICLES
      IF (spec_sum .AND. io_list(ispecies)%tracer) CYCLE
#endif
      current => io_list(ispecies)%attached_list%head
      sqrt_part_m  = SQRT(io_list(ispecies)%mass)

      DO WHILE(ASSOCIATED(current))
#ifdef PER_PARTICLE_CHARGE_MASS
        sqrt_part_m  = SQRT(current%mass)
#endif
        ! Copy the particle properties out for speed
        part_pmx = current%part_p(1) / sqrt_part_m
        part_pmy = current%part_p(2) / sqrt_part_m
        part_pmz = current%part_p(3) / sqrt_part_m

#include "particle_to_grid.inc"

        DO iz = sf_min, sf_max
        DO iy = sf_min, sf_max
        DO ix = sf_min, sf_max
          gf = gx(ix) * gy(iy) * gz(iz)
          sigma(cell_x+ix, cell_y+iy, cell_z+iz) = &
              sigma(cell_x+ix, cell_y+iy, cell_z+iz) + gf &
              * ((part_pmx - meanx(cell_x+ix, cell_y+iy, cell_z+iz))**2 &
              + (part_pmy - meany(cell_x+ix, cell_y+iy, cell_z+iz))**2 &
              + (part_pmz - meanz(cell_x+ix, cell_y+iy, cell_z+iz))**2)
          part_count(cell_x+ix, cell_y+iy, cell_z+iz) = &
              part_count(cell_x+ix, cell_y+iy, cell_z+iz) + gf
        ENDDO
        ENDDO
        ENDDO
        current => current%next
      ENDDO
    ENDDO

    CALL calc_boundary(sigma)
    CALL calc_boundary(part_count)

    ! 3/2 kT = <p^2>/(2m)
    sigma = sigma / MAX(part_count, 1.e-6_num) / kb / 3.0_num

    DEALLOCATE(part_count, meanx, meany, meanz)

  END SUBROUTINE calc_temperature



  SUBROUTINE calc_on_grid_with_evaluator(data_array, current_species, evaluator)

    REAL(num), DIMENSION(-2:,-2:,-2:), INTENT(OUT) :: data_array
    INTEGER, INTENT(IN) :: current_species
    ! The data to be weighted onto the grid
    REAL(num) :: wdata
    INTEGER :: ispecies, ix, iy, iz, spec_start, spec_end
    TYPE(particle), POINTER :: current
    LOGICAL :: spec_sum
#include "particle_head.inc"

    INTERFACE
      FUNCTION evaluator(a_particle, species_eval)
        USE shared_data
        TYPE(particle), POINTER :: a_particle
        INTEGER, INTENT(IN) :: species_eval
        REAL(num) :: evaluator
      END FUNCTION evaluator
    END INTERFACE

    data_array = 0.0_num

    spec_start = current_species
    spec_end = current_species
    spec_sum = .FALSE.

    IF (current_species <= 0) THEN
      spec_start = 1
      spec_end = n_species
      spec_sum = .TRUE.
    ENDIF

    DO ispecies = spec_start, spec_end
#ifndef NO_TRACER_PARTICLES
      IF (spec_sum .AND. io_list(ispecies)%tracer) CYCLE
#endif
      current => io_list(ispecies)%attached_list%head
      DO WHILE (ASSOCIATED(current))
#include "particle_to_grid.inc"

        wdata = evaluator(current, ispecies)
        DO iz = sf_min, sf_max
        DO iy = sf_min, sf_max
        DO ix = sf_min, sf_max
          data_array(cell_x+ix, cell_y+iy, cell_z+iz) = &
              data_array(cell_x+ix, cell_y+iy, cell_z+iz) &
              + gx(ix) * gy(iy) * gz(iz) * wdata
        ENDDO
        ENDDO
        ENDDO

        current => current%next
      ENDDO
    ENDDO

    CALL calc_boundary(data_array)
    DO ix = 1, 2*c_ndims
      CALL field_zero_gradient(data_array, c_stagger_centre, ix)
    ENDDO

  END SUBROUTINE calc_on_grid_with_evaluator



  SUBROUTINE calc_per_species_current(data_array, current_species, direction)

    REAL(num), DIMENSION(-2:,-2:,-2:), INTENT(INOUT) :: data_array
    INTEGER, INTENT(IN) :: current_species, direction

    REAL(num) :: part_px, part_py, part_pz
    REAL(num) :: part_q, part_mc, part_w
    REAL(num) :: part_j = 0.0_num
    REAL(num) :: idx, root, fac
    INTEGER :: ispecies, spec_start, spec_end, ix, iy, iz
    INTEGER(i8) :: ipart
    LOGICAL :: spec_sum

    TYPE (particle), POINTER :: current, next
#include "particle_head.inc"

    data_array = 0.0_num

    idx = 1.0_num / dx / dy / dz

    spec_start = current_species
    spec_end = current_species
    spec_sum = .FALSE.

    IF (current_species <= 0) THEN
      spec_start = 1
      spec_end = n_species
      spec_sum = .TRUE.
    ENDIF

    DO ispecies = spec_start, spec_end
      IF (io_list(ispecies)%species_type == c_species_id_photon) CYCLE
#ifndef NO_TRACER_PARTICLES
      IF (spec_sum .AND. io_list(ispecies)%tracer) CYCLE
#endif
      current => io_list(ispecies)%attached_list%head

      part_w = io_list(ispecies)%weight
      part_q  = io_list(ispecies)%charge
      part_mc = c * io_list(ispecies)%mass
      fac = part_q * part_w

      DO ipart = 1, io_list(ispecies)%attached_list%count
        next => current%next
#ifndef PER_SPECIES_WEIGHT
        part_w = current%weight
#ifdef PER_PARTICLE_CHARGE_MASS
        part_q  = current%charge
        part_mc = c * current%mass
#endif
        fac = part_q * part_w
#else
#ifdef PER_PARTICLE_CHARGE_MASS
        part_q  = current%charge
        part_mc = c * current%mass
        fac = part_q * part_w
#endif
#endif

        ! Copy the particle properties out for speed
        part_px = current%part_p(1)
        part_py = current%part_p(2)
        part_pz = current%part_p(3)
        root = 1.0_num / SQRT(part_mc**2 + part_px**2 + part_py**2 + part_pz**2)
        SELECT CASE (direction)
          CASE(c_dir_x)
            part_j = fac * part_px * root
          CASE(c_dir_y)
            part_j = fac * part_py * root
          CASE(c_dir_z)
            part_j = fac * part_pz * root
        END SELECT

#include "particle_to_grid.inc"

        DO iz = sf_min, sf_max
        DO iy = sf_min, sf_max
        DO ix = sf_min, sf_max
          data_array(cell_x+ix, cell_y+iy, cell_z+iz) = &
              data_array(cell_x+ix, cell_y+iy, cell_z+iz) &
              + gx(ix) * gy(iy) * gz(iz) * part_j
        ENDDO
        ENDDO
        ENDDO
        current => next
      ENDDO
    ENDDO

    fac = c * idx
    data_array = data_array * fac
    CALL processor_summation_bcs(data_array, ng, direction)

  END SUBROUTINE calc_per_species_current



  SUBROUTINE calc_per_species_jx(data_array, current_species)

    REAL(num), DIMENSION(-2:,-2:,-2:), INTENT(OUT) :: data_array
    INTEGER, INTENT(IN) :: current_species

    CALL calc_per_species_current(data_array, current_species, c_dir_x)

  END SUBROUTINE calc_per_species_jx



  SUBROUTINE calc_per_species_jy(data_array, current_species)

    REAL(num), DIMENSION(-2:,-2:,-2:), INTENT(OUT) :: data_array
    INTEGER, INTENT(IN) :: current_species

    CALL calc_per_species_current(data_array, current_species, c_dir_y)

  END SUBROUTINE calc_per_species_jy



  SUBROUTINE calc_per_species_jz(data_array, current_species)

    REAL(num), DIMENSION(-2:,-2:,-2:), INTENT(OUT) :: data_array
    INTEGER, INTENT(IN) :: current_species

    CALL calc_per_species_current(data_array, current_species, c_dir_z)

  END SUBROUTINE calc_per_species_jz



  SUBROUTINE calc_total_energy_sum

    REAL(num) :: particle_energy, field_energy
    REAL(num) :: part_ux, part_uy, part_uz
    REAL(num) :: part_mc, part_w, fac, gamma
    REAL(num) :: sum_out(2), sum_in(2)
    REAL(num), PARAMETER :: c2 = c**2
    INTEGER :: ispecies, i, j, k
    TYPE(particle), POINTER :: current

    particle_energy = 0.0_num

    ! Sum over all particles to calculate total kinetic energy
    DO ispecies = 1, n_species
#ifndef NO_TRACER_PARTICLES
      IF (species_list(ispecies)%tracer) CYCLE
#endif
      current => species_list(ispecies)%attached_list%head
      part_mc = c * species_list(ispecies)%mass
      part_w = species_list(ispecies)%weight
      fac = part_mc * part_w * c

      DO WHILE (ASSOCIATED(current))
        ! Copy the particle properties out for speed
#ifdef PER_PARTICLE_CHARGE_MASS
        part_mc = c * current%mass
#ifndef PER_SPECIES_WEIGHT
        part_w = current%weight
#endif
        fac = part_mc * part_w * c
#else
#ifndef PER_SPECIES_WEIGHT
        part_w = current%weight
        fac = part_mc * part_w * c
#endif
#endif

        IF (species_list(ispecies)%species_type /= c_species_id_photon) THEN
          part_ux = current%part_p(1) / part_mc
          part_uy = current%part_p(2) / part_mc
          part_uz = current%part_p(3) / part_mc
          gamma = SQRT(part_ux**2 + part_uy**2 + part_uz**2 + 1.0_num)
          particle_energy = particle_energy + (gamma - 1.0_num) * fac
#ifdef PHOTONS
        ELSE
          particle_energy = particle_energy + current%particle_energy * part_w
#endif
        ENDIF

        current => current%next
      ENDDO
    ENDDO

    ! EM field energy
    field_energy = 0.0_num
    DO k = 1, nz
    DO j = 1, ny
    DO i = 1, nx
      field_energy = field_energy + ex(i,j,k)**2 + ey(i,j,k)**2 &
          + ez(i,j,k)**2 + c2 * (bx(i,j,k)**2 + by(i,j,k)**2 + bz(i,j,k)**2)
    ENDDO
    ENDDO
    ENDDO
    field_energy = 0.5_num * epsilon0 * field_energy * dx * dy * dz

    sum_out(1) = particle_energy
    sum_out(2) = field_energy
    CALL MPI_REDUCE(sum_out, sum_in, 2, mpireal, MPI_SUM, 0, comm, errcode)
    total_particle_energy = sum_in(1)
    total_field_energy = sum_in(2)

  END SUBROUTINE calc_total_energy_sum

END MODULE calc_df
