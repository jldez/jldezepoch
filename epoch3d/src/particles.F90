! Copyright (C) 2010-2015 Keith Bennett <K.Bennett@warwick.ac.uk>
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

MODULE particles

  USE current_smooth
  USE boundary
  USE partlist
#ifdef PREFETCH
  USE prefetch
#endif

  IMPLICIT NONE

CONTAINS

  SUBROUTINE push_particles

    ! 2nd order accurate particle pusher using parabolic weighting
    ! on and off the grid. The calculation of J looks rather odd
    ! Since it works by solving d(rho)/dt = div(J) and doing a 1st order
    ! Estimate of rho(t+1.5*dt) rather than calculating J directly
    ! This gives exact charge conservation on the grid

    ! Contains the integer cell position of the particle in x, y, z
    INTEGER :: cell_x1, cell_x2, cell_x3
    INTEGER :: cell_y1, cell_y2, cell_y3
    INTEGER :: cell_z1, cell_z2, cell_z3

    ! Xi (space factor see page 38 in manual)
    ! The code now uses gx and hx instead of xi0 and xi1

    ! J from a given particle, can be spread over up to 3 cells in
    ! Each direction due to parabolic weighting. We allocate 4 or 5
    ! Cells because the position of the particle at t = t+1.5dt is not
    ! known until later. This part of the algorithm could probably be
    ! Improved, but at the moment, this is just a straight copy of
    ! The core of the PSC algorithm
    INTEGER, PARAMETER :: sf0 = sf_min, sf1 = sf_max
    REAL(num) :: jxh
    REAL(num), DIMENSION(sf0-1:sf1+1) :: jyh
    REAL(num), DIMENSION(sf0-1:sf1+1,sf0-1:sf1+1) :: jzh

    ! Properties of the current particle. Copy out of particle arrays for speed
    REAL(num) :: part_x, part_y, part_z
    REAL(num) :: part_ux, part_uy, part_uz
    REAL(num) :: part_q, part_mc, ipart_mc, part_weight

    ! Used for particle probes (to see of probe conditions are satisfied)
#ifndef NO_PARTICLE_PROBES
    REAL(num) :: init_part_x, final_part_x
    REAL(num) :: init_part_y, final_part_y
    REAL(num) :: init_part_z, final_part_z
    TYPE(particle_probe), POINTER :: current_probe
    TYPE(particle), POINTER :: particle_copy
    REAL(num) :: d_init, d_final
    REAL(num) :: probe_energy, part_mc2
#endif

    ! Contains the floating point version of the cell number (never actually
    ! used)
    REAL(num) :: cell_x_r, cell_y_r, cell_z_r

    ! The fraction of a cell between the particle position and the cell boundary
    REAL(num) :: cell_frac_x, cell_frac_y, cell_frac_z

    ! Weighting factors as Eqn 4.77 page 25 of manual
    ! Eqn 4.77 would be written as
    ! F(j-1) * gmx + F(j) * g0x + F(j+1) * gpx
    ! Defined at the particle position
    REAL(num), DIMENSION(sf_min-1:sf_max+1) :: gx, gy, gz

    ! Defined at the particle position - 0.5 grid cell in each direction
    ! This is to deal with the grid stagger
    REAL(num), DIMENSION(sf_min-1:sf_max+1) :: hx, hy, hz

    ! Fields at particle location
    REAL(num) :: ex_part, ey_part, ez_part, bx_part, by_part, bz_part

    ! P+, P- and Tau variables from Boris1970, page27 of manual
    REAL(num) :: uxp, uxm, uyp, uym, uzp, uzm
    REAL(num) :: tau, taux, tauy, tauz, taux2, tauy2, tauz2

    ! charge to mass ratio modified by normalisation
    REAL(num) :: cmratio, ccmratio

    ! Used by J update
    INTEGER :: xmin, xmax, ymin, ymax, zmin, zmax
    REAL(num) :: wx, wy, wz

    ! Temporary variables
    REAL(num) :: idx, idy, idz
    REAL(num) :: idtyz, idtxz, idtxy
    REAL(num) :: idt, dto2, dtco2
    REAL(num) :: fcx, fcy, fcz, fjx, fjy, fjz
    REAL(num) :: root, dtfac, gamma, third
    REAL(num) :: delta_x, delta_y, delta_z
    REAL(num) :: xfac1, xfac2, yfac1, yfac2, zfac1, zfac2
    REAL(num) :: gz_iz, hz_iz, hygz, hyhz, hzyfac1, hzyfac2, yzfac
    INTEGER :: ispecies, ix, iy, iz, dcellx, dcelly, dcellz, cx, cy, cz
    INTEGER(i8) :: ipart
#ifndef NO_PARTICLE_PROBES
    LOGICAL :: probes_for_species
#endif
#ifndef NO_TRACER_PARTICLES
    LOGICAL :: not_tracer_species
#endif
    ! Particle weighting multiplication factor
#ifdef PARTICLE_SHAPE_BSPLINE3
    REAL(num) :: cf2
    REAL(num), PARAMETER :: fac = (1.0_num / 24.0_num)**c_ndims
#elif  PARTICLE_SHAPE_TOPHAT
    REAL(num), PARAMETER :: fac = 1.0_num
#else
    REAL(num) :: cf2
    REAL(num), PARAMETER :: fac = (0.5_num)**c_ndims
#endif

    TYPE(particle), POINTER :: current, next

#ifdef PREFETCH
    CALL prefetch_particle(species_list(1)%attached_list%head)
#endif

    jx = 0.0_num
    jy = 0.0_num
    jz = 0.0_num

    gx = 0.0_num
    gy = 0.0_num
    gz = 0.0_num

    ! Unvarying multiplication factors

    idx = 1.0_num / dx
    idy = 1.0_num / dy
    idz = 1.0_num / dz
    idt = 1.0_num / dt
    dto2 = dt / 2.0_num
    dtco2 = c * dto2
    dtfac = 0.5_num * dt * fac
    third = 1.0_num / 3.0_num

    idtyz = idt * idy * idz * fac
    idtxz = idt * idx * idz * fac
    idtxy = idt * idx * idy * fac

    DO ispecies = 1, n_species
      current => species_list(ispecies)%attached_list%head
      IF (species_list(ispecies)%immobile) CYCLE
      IF (species_list(ispecies)%species_type == c_species_id_photon) THEN
#ifdef PHOTONS
        IF (photon_dynamics) CALL push_photons(ispecies)
#endif
        CYCLE
      ENDIF
#ifndef NO_PARTICLE_PROBES
      current_probe => species_list(ispecies)%attached_probes
      probes_for_species = ASSOCIATED(current_probe)
#endif
#ifndef NO_TRACER_PARTICLES
      not_tracer_species = .NOT. species_list(ispecies)%tracer
#endif

#ifdef PER_SPECIES_WEIGHT
      part_weight = species_list(ispecies)%weight
      fcx = idtyz * part_weight
      fcy = idtxz * part_weight
      fcz = idtxy * part_weight
#endif
#ifndef PER_PARTICLE_CHARGE_MASS
      part_q   = species_list(ispecies)%charge
      part_mc  = c * species_list(ispecies)%mass
      ipart_mc = 1.0_num / part_mc
      cmratio  = part_q * dtfac * ipart_mc
      ccmratio = c * cmratio
#ifndef NO_PARTICLE_PROBES
      part_mc2 = c * part_mc
#endif
#endif
      !DEC$ VECTOR ALWAYS
      DO ipart = 1, species_list(ispecies)%attached_list%count
        next => current%next
#ifdef PREFETCH
        CALL prefetch_particle(next)
#endif
#ifndef PER_SPECIES_WEIGHT
        part_weight = current%weight
        fcx = idtyz * part_weight
        fcy = idtxz * part_weight
        fcz = idtxy * part_weight
#endif
#ifndef NO_PARTICLE_PROBES
        init_part_x = current%part_pos(1)
        init_part_y = current%part_pos(2)
        init_part_z = current%part_pos(3)
#endif
#ifdef PER_PARTICLE_CHARGE_MASS
        part_q   = current%charge
        part_mc  = c * current%mass
        ipart_mc = 1.0_num / part_mc
        cmratio  = part_q * dtfac * ipart_mc
        ccmratio = c * cmratio
#ifndef NO_PARTICLE_PROBES
        part_mc2 = c * part_mc
#endif
#endif
        ! Copy the particle properties out for speed
        part_x  = current%part_pos(1) - x_grid_min_local
        part_y  = current%part_pos(2) - y_grid_min_local
        part_z  = current%part_pos(3) - z_grid_min_local
        part_ux = current%part_p(1) * ipart_mc
        part_uy = current%part_p(2) * ipart_mc
        part_uz = current%part_p(3) * ipart_mc

        ! Calculate v(t) from p(t)
        ! See PSC manual page (25-27)
        root = dtco2 / SQRT(part_ux**2 + part_uy**2 + part_uz**2 + 1.0_num)

        ! Move particles to half timestep position to first order
        part_x = part_x + part_ux * root
        part_y = part_y + part_uy * root
        part_z = part_z + part_uz * root

        ! Grid cell position as a fraction.
#ifdef PARTICLE_SHAPE_TOPHAT
        cell_x_r = part_x * idx - 0.5_num
        cell_y_r = part_y * idy - 0.5_num
        cell_z_r = part_z * idz - 0.5_num
#else
        cell_x_r = part_x * idx
        cell_y_r = part_y * idy
        cell_z_r = part_z * idz
#endif
        ! Round cell position to nearest cell
        cell_x1 = FLOOR(cell_x_r + 0.5_num)
        ! Calculate fraction of cell between nearest cell boundary and particle
        cell_frac_x = REAL(cell_x1, num) - cell_x_r
        cell_x1 = cell_x1 + 1

        cell_y1 = FLOOR(cell_y_r + 0.5_num)
        cell_frac_y = REAL(cell_y1, num) - cell_y_r
        cell_y1 = cell_y1 + 1

        cell_z1 = FLOOR(cell_z_r + 0.5_num)
        cell_frac_z = REAL(cell_z1, num) - cell_z_r
        cell_z1 = cell_z1 + 1

        ! Particle weight factors as described in the manual, page25
        ! These weight grid properties onto particles
        ! Also used to weight particle properties onto grid, used later
        ! to calculate J
        ! NOTE: These weights require an additional multiplication factor!
#ifdef PARTICLE_SHAPE_BSPLINE3
#include "bspline3/gx.inc"
#elif  PARTICLE_SHAPE_TOPHAT
#include "tophat/gx.inc"
#else
#include "triangle/gx.inc"
#endif

        ! Now redo shifted by half a cell due to grid stagger.
        ! Use shifted version for ex in X, ey in Y, ez in Z
        ! And in Y&Z for bx, X&Z for by, X&Y for bz
        cell_x2 = FLOOR(cell_x_r)
        cell_frac_x = REAL(cell_x2, num) - cell_x_r + 0.5_num
        cell_x2 = cell_x2 + 1

        cell_y2 = FLOOR(cell_y_r)
        cell_frac_y = REAL(cell_y2, num) - cell_y_r + 0.5_num
        cell_y2 = cell_y2 + 1

        cell_z2 = FLOOR(cell_z_r)
        cell_frac_z = REAL(cell_z2, num) - cell_z_r + 0.5_num
        cell_z2 = cell_z2 + 1

        dcellx = 0
        dcelly = 0
        dcellz = 0
        ! NOTE: These weights require an additional multiplication factor!
#ifdef PARTICLE_SHAPE_BSPLINE3
#include "bspline3/hx_dcell.inc"
#elif  PARTICLE_SHAPE_TOPHAT
#include "tophat/hx_dcell.inc"
#else
#include "triangle/hx_dcell.inc"
#endif

        ! These are the electric and magnetic fields interpolated to the
        ! particle position. They have been checked and are correct.
        ! Actually checking this is messy.
#ifdef PARTICLE_SHAPE_BSPLINE3
#include "bspline3/e_part.inc"
#include "bspline3/b_part.inc"
#elif  PARTICLE_SHAPE_TOPHAT
#include "tophat/e_part.inc"
#include "tophat/b_part.inc"
#else
#include "triangle/e_part.inc"
#include "triangle/b_part.inc"
#endif

        ! update particle momenta using weighted fields
        uxm = part_ux + cmratio * ex_part
        uym = part_uy + cmratio * ey_part
        uzm = part_uz + cmratio * ez_part

        ! Half timestep, then use Boris1970 rotation, see Birdsall and Langdon
        root = ccmratio / SQRT(uxm**2 + uym**2 + uzm**2 + 1.0_num)

        taux = bx_part * root
        tauy = by_part * root
        tauz = bz_part * root

        taux2 = taux**2
        tauy2 = tauy**2
        tauz2 = tauz**2

        tau = 1.0_num / (1.0_num + taux2 + tauy2 + tauz2)

        uxp = ((1.0_num + taux2 - tauy2 - tauz2) * uxm &
            + 2.0_num * ((taux * tauy + tauz) * uym &
            + (taux * tauz - tauy) * uzm)) * tau
        uyp = ((1.0_num - taux2 + tauy2 - tauz2) * uym &
            + 2.0_num * ((tauy * tauz + taux) * uzm &
            + (tauy * taux - tauz) * uxm)) * tau
        uzp = ((1.0_num - taux2 - tauy2 + tauz2) * uzm &
            + 2.0_num * ((tauz * taux + tauy) * uxm &
            + (tauz * tauy - taux) * uym)) * tau

        ! Rotation over, go to full timestep
        part_ux = uxp + cmratio * ex_part
        part_uy = uyp + cmratio * ey_part
        part_uz = uzp + cmratio * ez_part

        ! Calculate particle velocity from particle momentum
        gamma = SQRT(part_ux**2 + part_uy**2 + part_uz**2 + 1.0_num)
        root = dtco2 / gamma

        delta_x = part_ux * root
        delta_y = part_uy * root
        delta_z = part_uz * root

        ! Move particles to end of time step at 2nd order accuracy
        part_x = part_x + delta_x
        part_y = part_y + delta_y
        part_z = part_z + delta_z

        ! particle has now finished move to end of timestep, so copy back
        ! into particle array
        current%part_pos = (/ part_x + x_grid_min_local, &
            part_y + y_grid_min_local, part_z + z_grid_min_local /)
        current%part_p   = part_mc * (/ part_ux, part_uy, part_uz /)

#ifndef NO_PARTICLE_PROBES
        final_part_x = current%part_pos(1)
        final_part_y = current%part_pos(2)
        final_part_z = current%part_pos(3)
#endif
        ! Original code calculates densities of electrons, ions and neutrals
        ! here. This has been removed to reduce memory footprint

        ! If the code is compiled with tracer particle support then put in an
        ! IF statement so that the current is not calculated for this species
#ifndef NO_TRACER_PARTICLES
        IF (not_tracer_species) THEN
#endif
          ! Now advance to t+1.5dt to calculate current. This is detailed in
          ! the manual between pages 37 and 41. The version coded up looks
          ! completely different to that in the manual, but is equivalent.
          ! Use t+1.5 dt so that can update J to t+dt at 2nd order
          part_x = part_x + delta_x
          part_y = part_y + delta_y
          part_z = part_z + delta_z

#ifdef PARTICLE_SHAPE_TOPHAT
          cell_x_r = part_x * idx - 0.5_num
          cell_y_r = part_y * idy - 0.5_num
          cell_z_r = part_z * idz - 0.5_num
#else
          cell_x_r = part_x * idx
          cell_y_r = part_y * idy
          cell_z_r = part_z * idz
#endif
          cell_x3 = FLOOR(cell_x_r + 0.5_num)
          cell_frac_x = REAL(cell_x3, num) - cell_x_r
          cell_x3 = cell_x3 + 1

          cell_y3 = FLOOR(cell_y_r + 0.5_num)
          cell_frac_y = REAL(cell_y3, num) - cell_y_r
          cell_y3 = cell_y3 + 1

          cell_z3 = FLOOR(cell_z_r + 0.5_num)
          cell_frac_z = REAL(cell_z3, num) - cell_z_r
          cell_z3 = cell_z3 + 1

          hx = 0.0_num
          hy = 0.0_num
          hz = 0.0_num

          dcellx = cell_x3 - cell_x1
          dcelly = cell_y3 - cell_y1
          dcellz = cell_z3 - cell_z1
          ! NOTE: These weights require an additional multiplication factor!
#ifdef PARTICLE_SHAPE_BSPLINE3
#include "bspline3/hx_dcell.inc"
#elif  PARTICLE_SHAPE_TOPHAT
#include "tophat/hx_dcell.inc"
#else
#include "triangle/hx_dcell.inc"
#endif

          ! Now change Xi1* to be Xi1*-Xi0*. This makes the representation of
          ! the current update much simpler
          hx = hx - gx
          hy = hy - gy
          hz = hz - gz

          ! Remember that due to CFL condition particle can never cross more
          ! than one gridcell in one timestep

          xmin = sf_min + (dcellx - 1) / 2
          xmax = sf_max + (dcellx + 1) / 2

          ymin = sf_min + (dcelly - 1) / 2
          ymax = sf_max + (dcelly + 1) / 2

          zmin = sf_min + (dcellz - 1) / 2
          zmax = sf_max + (dcellz + 1) / 2

          fjx = fcx * part_q
          fjy = fcy * part_q
          fjz = fcz * part_q

          jzh = 0.0_num
          DO iz = zmin, zmax
            cz = cell_z1 + iz
            zfac1 =         gz(iz) + 0.5_num * hz(iz)
            zfac2 = third * hz(iz) + 0.5_num * gz(iz)

            gz_iz = gz(iz)
            hz_iz = hz(iz)

            jyh = 0.0_num
            DO iy = ymin, ymax
              cy = cell_y1 + iy
              yfac1 =         gy(iy) + 0.5_num * hy(iy)
              yfac2 = third * hy(iy) + 0.5_num * gy(iy)

              hygz = hy(iy) * gz_iz
              hyhz = hy(iy) * hz_iz
              yzfac = gy(iy) * zfac1 + hy(iy) * zfac2
              hzyfac1 = hz_iz * yfac1
              hzyfac2 = hz_iz * yfac2

              jxh = 0.0_num
              DO ix = xmin, xmax
                cx = cell_x1 + ix
                xfac1 =         gx(ix) + 0.5_num * hx(ix)
                xfac2 = third * hx(ix) + 0.5_num * gx(ix)

                wx = hx(ix) * yzfac
                wy = xfac1 * hygz + xfac2 * hyhz
                wz = gx(ix) * hzyfac1 + hx(ix) * hzyfac2

                ! This is the bit that actually solves d(rho)/dt = -div(J)
                jxh = jxh - fjx * wx
                jyh(ix) = jyh(ix) - fjy * wy
                jzh(ix, iy) = jzh(ix, iy) - fjz * wz

                jx(cx, cy, cz) = jx(cx, cy, cz) + jxh
                jy(cx, cy, cz) = jy(cx, cy, cz) + jyh(ix)
                jz(cx, cy, cz) = jz(cx, cy, cz) + jzh(ix, iy)
              ENDDO
            ENDDO
          ENDDO
#ifndef NO_TRACER_PARTICLES
        ENDIF
#endif
#ifndef NO_PARTICLE_PROBES
        IF (probes_for_species) THEN
          ! Compare the current particle with the parameters of any probes in
          ! the system. These particles are copied into a separate part of the
          ! output file.

          current_probe => species_list(ispecies)%attached_probes

          ! Cycle through probes
          DO WHILE(ASSOCIATED(current_probe))
            ! Note that this is the energy of a single REAL particle in the
            ! pseudoparticle, NOT the energy of the pseudoparticle
            probe_energy = (gamma - 1.0_num) * part_mc2

            ! Unidirectional probe
            IF (probe_energy > current_probe%ek_min) THEN
              IF (probe_energy < current_probe%ek_max) THEN

                d_init  = SUM(current_probe%normal * (current_probe%point &
                    - (/init_part_x, init_part_y, init_part_z/)))
                d_final = SUM(current_probe%normal * (current_probe%point &
                    - (/final_part_x, final_part_y, final_part_z/)))
                IF (d_final < 0.0_num .AND. d_init >= 0.0_num) THEN
                  ! this particle is wanted so copy it to the list associated
                  ! with this probe
                  ALLOCATE(particle_copy)
                  particle_copy = current
                  CALL add_particle_to_partlist(&
                      current_probe%sampled_particles, particle_copy)
                  NULLIFY(particle_copy)
                ENDIF

              ENDIF
            ENDIF
            current_probe => current_probe%next
          ENDDO
        ENDIF
#endif
        current => next
      ENDDO
    ENDDO

    IF (.NOT.use_field_ionisation) THEN
      CALL current_bcs
      CALL particle_bcs

      IF (smooth_currents) CALL smooth_current()
    END IF

  END SUBROUTINE push_particles

#ifdef PHOTONS
  SUBROUTINE push_photons(ispecies)

    ! Very simple photon pusher
    ! Properties of the current particle. Copy out of particle arrays for speed
    REAL(num) :: delta_x, delta_y, delta_z
    INTEGER,INTENT(IN) :: ispecies
    TYPE(particle), POINTER :: current

    ! Used for particle probes (to see of probe conditions are satisfied)
#ifndef NO_PARTICLE_PROBES
    REAL(num) :: init_part_x, final_part_x
    REAL(num) :: init_part_y, final_part_y
    REAL(num) :: init_part_z, final_part_z
    TYPE(particle_probe), POINTER :: current_probe
    TYPE(particle), POINTER :: particle_copy
    REAL(num) :: d_init, d_final
    REAL(num) :: probe_energy, dtfac, fac
    LOGICAL :: probes_for_species
#endif

#ifndef NO_PARTICLE_PROBES
    current_probe => species_list(ispecies)%attached_probes
    probes_for_species = ASSOCIATED(current_probe)
#endif
    dtfac = dt * c**2

    ! set current to point to head of list
    current => species_list(ispecies)%attached_list%head
    ! loop over photons
    DO WHILE(ASSOCIATED(current))
      ! Note that this is the energy of a single REAL particle in the
      ! pseudoparticle, NOT the energy of the pseudoparticle
      probe_energy = current%particle_energy

      fac = dtfac / probe_energy
      delta_x = current%part_p(1) * fac
      delta_y = current%part_p(2) * fac
      delta_z = current%part_p(3) * fac
#ifndef NO_PARTICLE_PROBES
      init_part_x = current%part_pos(1)
      init_part_y = current%part_pos(2)
      init_part_z = current%part_pos(3)
#endif
      current%part_pos = current%part_pos + (/delta_x, delta_y, delta_z/)
#ifndef NO_PARTICLE_PROBES
      final_part_x = current%part_pos(1)
      final_part_y = current%part_pos(2)
      final_part_z = current%part_pos(3)
#endif

#ifndef NO_PARTICLE_PROBES
      IF (probes_for_species) THEN
        ! Compare the current particle with the parameters of any probes in
        ! the system. These particles are copied into a separate part of the
        ! output file.

        current_probe => species_list(ispecies)%attached_probes

        ! Cycle through probes
        DO WHILE(ASSOCIATED(current_probe))
          ! Unidirectional probe
          IF (probe_energy > current_probe%ek_min) THEN
            IF (probe_energy < current_probe%ek_max) THEN

              d_init  = SUM(current_probe%normal * (current_probe%point &
                  - (/init_part_x, init_part_y, init_part_z/)))
              d_final = SUM(current_probe%normal * (current_probe%point &
                  - (/final_part_x, final_part_y, final_part_z/)))
              IF (d_final < 0.0_num .AND. d_init >= 0.0_num) THEN
                ! this particle is wanted so copy it to the list associated
                ! with this probe
                ALLOCATE(particle_copy)
                particle_copy = current
                CALL add_particle_to_partlist(&
                    current_probe%sampled_particles, particle_copy)
                NULLIFY(particle_copy)
              ENDIF

            ENDIF
          ENDIF
          current_probe => current_probe%next
        ENDDO
      ENDIF
#endif

      current => current%next
    ENDDO

  END SUBROUTINE push_photons
#endif

END MODULE particles
