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

MODULE fields

  USE boundary

  IMPLICIT NONE

  INTEGER :: field_order
  REAL(num) :: hdt, fac
  REAL(num) :: hdtx, hdty, hdtz
  REAL(num) :: cnx, cny, cnz

CONTAINS

  SUBROUTINE set_field_order(order)

    INTEGER, INTENT(IN) :: order

    field_order = order
    fng = field_order / 2

    IF (field_order == 2) THEN
      cfl = 1.0_num
    ELSE IF (field_order == 4) THEN
      cfl = 6.0_num / 7.0_num
    ELSE
      cfl = 120.0_num / 149.0_num
    ENDIF

  END SUBROUTINE set_field_order



  SUBROUTINE update_e_field

    INTEGER :: ix, iy, iz
    REAL(num) :: cpml_x, cpml_y, cpml_z
    REAL(num) :: c1, c2, c3
    REAL(num) :: cx1, cx2, cx3
    REAL(num) :: cy1, cy2, cy3
    REAL(num) :: cz1, cz2, cz3

    IF (cpml_boundaries) THEN
      cpml_x = cnx
      cpml_y = cny
      cpml_z = cnz

      IF (field_order == 2) THEN
        DO iz = 1, nz
          cpml_z = cnz / cpml_kappa_ez(iz)
          DO iy = 1, ny
            cpml_y = cny / cpml_kappa_ey(iy)
            DO ix = 1, nx
              cpml_x = cnx / cpml_kappa_ex(ix)

              ex(ix, iy, iz) = ex(ix, iy, iz) &
                  + cpml_y * (bz(ix  , iy  , iz  ) - bz(ix  , iy-1, iz  )) &
                  - cpml_z * (by(ix  , iy  , iz  ) - by(ix  , iy  , iz-1)) &
                  - fac * jx(ix, iy, iz)

              ey(ix, iy, iz) = ey(ix, iy, iz) &
                  + cpml_z * (bx(ix  , iy  , iz  ) - bx(ix  , iy  , iz-1)) &
                  - cpml_x * (bz(ix  , iy  , iz  ) - bz(ix-1, iy  , iz  )) &
                  - fac * jy(ix, iy, iz)

              ez(ix, iy, iz) = ez(ix, iy, iz) &
                  + cpml_x * (by(ix  , iy  , iz  ) - by(ix-1, iy  , iz  )) &
                  - cpml_y * (bx(ix  , iy  , iz  ) - bx(ix  , iy-1, iz  )) &
                  - fac * jz(ix, iy, iz)
            ENDDO
          ENDDO
        ENDDO
      ELSE IF (field_order == 4) THEN
        c1 = 9.0_num / 8.0_num
        c2 = -1.0_num / 24.0_num

        DO iz = 1, nz
          cpml_z = cnz / cpml_kappa_ez(iz)
          cz1 = c1 * cpml_z
          cz2 = c2 * cpml_z
          DO iy = 1, ny
            cpml_y = cny / cpml_kappa_ey(iy)
            cy1 = c1 * cpml_y
            cy2 = c2 * cpml_y
            DO ix = 1, nx
              cpml_x = cnx / cpml_kappa_ex(ix)
              cx1 = c1 * cpml_x
              cx2 = c2 * cpml_x

              ex(ix, iy, iz) = ex(ix, iy, iz) &
                  + cy1 * (bz(ix  , iy  , iz  ) - bz(ix  , iy-1, iz  )) &
                  + cy2 * (bz(ix  , iy+1, iz  ) - bz(ix  , iy-2, iz  )) &
                  - cz1 * (by(ix  , iy  , iz  ) - by(ix  , iy  , iz-1)) &
                  - cz2 * (by(ix  , iy  , iz+1) - by(ix  , iy  , iz-2)) &
                  - fac * jx(ix, iy, iz)

              ey(ix, iy, iz) = ey(ix, iy, iz) &
                  + cz1 * (bx(ix  , iy  , iz  ) - bx(ix  , iy  , iz-1)) &
                  + cz2 * (bx(ix  , iy  , iz+1) - bx(ix  , iy  , iz-2)) &
                  - cx1 * (bz(ix  , iy  , iz  ) - bz(ix-1, iy  , iz  )) &
                  - cx2 * (bz(ix+1, iy  , iz  ) - bz(ix-2, iy  , iz  )) &
                  - fac * jy(ix, iy, iz)

              ez(ix, iy, iz) = ez(ix, iy, iz) &
                  + cx1 * (by(ix  , iy  , iz  ) - by(ix-1, iy  , iz  )) &
                  + cx2 * (by(ix+1, iy  , iz  ) - by(ix-2, iy  , iz  )) &
                  - cy1 * (bx(ix  , iy  , iz  ) - bx(ix  , iy-1, iz  )) &
                  - cy2 * (bx(ix  , iy+1, iz  ) - bx(ix  , iy-2, iz  )) &
                  - fac * jz(ix, iy, iz)
            ENDDO
          ENDDO
        ENDDO
      ELSE
        c1 = 75.0_num / 64.0_num
        c2 = -25.0_num / 384.0_num
        c3 = 3.0_num / 640.0_num

        DO iz = 1, nz
          cpml_z = cnz / cpml_kappa_ez(iz)
          cz1 = c1 * cpml_z
          cz2 = c2 * cpml_z
          cz3 = c3 * cpml_z
          DO iy = 1, ny
            cpml_y = cny / cpml_kappa_ey(iy)
            cy1 = c1 * cpml_y
            cy2 = c2 * cpml_y
            cy3 = c3 * cpml_y
            DO ix = 1, nx
              cpml_x = cnx / cpml_kappa_ex(ix)
              cx1 = c1 * cpml_x
              cx2 = c2 * cpml_x
              cx3 = c3 * cpml_x

              ex(ix, iy, iz) = ex(ix, iy, iz) &
                  + cy1 * (bz(ix  , iy  , iz  ) - bz(ix  , iy-1, iz  )) &
                  + cy2 * (bz(ix  , iy+1, iz  ) - bz(ix  , iy-2, iz  )) &
                  + cy3 * (bz(ix  , iy+2, iz  ) - bz(ix  , iy-3, iz  )) &
                  - cz1 * (by(ix  , iy  , iz  ) - by(ix  , iy  , iz-1)) &
                  - cz2 * (by(ix  , iy  , iz+1) - by(ix  , iy  , iz-2)) &
                  - cz3 * (by(ix  , iy  , iz+2) - by(ix  , iy  , iz-3)) &
                  - fac * jx(ix, iy, iz)

              ey(ix, iy, iz) = ey(ix, iy, iz) &
                  + cz1 * (bx(ix  , iy  , iz  ) - bx(ix  , iy  , iz-1)) &
                  + cz2 * (bx(ix  , iy  , iz+1) - bx(ix  , iy  , iz-2)) &
                  + cz3 * (bx(ix  , iy  , iz+2) - bx(ix  , iy  , iz-3)) &
                  - cx1 * (bz(ix  , iy  , iz  ) - bz(ix-1, iy  , iz  )) &
                  - cx2 * (bz(ix+1, iy  , iz  ) - bz(ix-2, iy  , iz  )) &
                  - cx3 * (bz(ix+2, iy  , iz  ) - bz(ix-3, iy  , iz  )) &
                  - fac * jy(ix, iy, iz)

              ez(ix, iy, iz) = ez(ix, iy, iz) &
                  + cx1 * (by(ix  , iy  , iz  ) - by(ix-1, iy  , iz  )) &
                  + cx2 * (by(ix+1, iy  , iz  ) - by(ix-2, iy  , iz  )) &
                  + cx3 * (by(ix+2, iy  , iz  ) - by(ix-3, iy  , iz  )) &
                  - cy1 * (bx(ix  , iy  , iz  ) - bx(ix  , iy-1, iz  )) &
                  - cy2 * (bx(ix  , iy+1, iz  ) - bx(ix  , iy-2, iz  )) &
                  - cy3 * (bx(ix  , iy+2, iz  ) - bx(ix  , iy-3, iz  )) &
                  - fac * jz(ix, iy, iz)
            ENDDO
          ENDDO
        ENDDO
      ENDIF

      CALL cpml_advance_e_currents(hdt)
    ELSE
      IF (field_order == 2) THEN
        DO iz = 1, nz
          DO iy = 1, ny
            DO ix = 1, nx

              ex(ix, iy, iz) = ex(ix, iy, iz) &
                  + cny * (bz(ix  , iy  , iz  ) - bz(ix  , iy-1, iz  )) &
                  - cnz * (by(ix  , iy  , iz  ) - by(ix  , iy  , iz-1)) &
                  - fac * jx(ix, iy, iz)

              ey(ix, iy, iz) = ey(ix, iy, iz) &
                  + cnz * (bx(ix  , iy  , iz  ) - bx(ix  , iy  , iz-1)) &
                  - cnx * (bz(ix  , iy  , iz  ) - bz(ix-1, iy  , iz  )) &
                  - fac * jy(ix, iy, iz)

              ez(ix, iy, iz) = ez(ix, iy, iz) &
                  + cnx * (by(ix  , iy  , iz  ) - by(ix-1, iy  , iz  )) &
                  - cny * (bx(ix  , iy  , iz  ) - bx(ix  , iy-1, iz  )) &
                  - fac * jz(ix, iy, iz)
            ENDDO
          ENDDO
        ENDDO
      ELSE IF (field_order == 4) THEN
        c1 = 9.0_num / 8.0_num
        c2 = -1.0_num / 24.0_num

        DO iz = 1, nz
          cz1 = c1 * cnz
          cz2 = c2 * cnz
          DO iy = 1, ny
            cy1 = c1 * cny
            cy2 = c2 * cny
            DO ix = 1, nx
              cx1 = c1 * cnx
              cx2 = c2 * cnx

              ex(ix, iy, iz) = ex(ix, iy, iz) &
                  + cy1 * (bz(ix  , iy  , iz  ) - bz(ix  , iy-1, iz  )) &
                  + cy2 * (bz(ix  , iy+1, iz  ) - bz(ix  , iy-2, iz  )) &
                  - cz1 * (by(ix  , iy  , iz  ) - by(ix  , iy  , iz-1)) &
                  - cz2 * (by(ix  , iy  , iz+1) - by(ix  , iy  , iz-2)) &
                  - fac * jx(ix, iy, iz)

              ey(ix, iy, iz) = ey(ix, iy, iz) &
                  + cz1 * (bx(ix  , iy  , iz  ) - bx(ix  , iy  , iz-1)) &
                  + cz2 * (bx(ix  , iy  , iz+1) - bx(ix  , iy  , iz-2)) &
                  - cx1 * (bz(ix  , iy  , iz  ) - bz(ix-1, iy  , iz  )) &
                  - cx2 * (bz(ix+1, iy  , iz  ) - bz(ix-2, iy  , iz  )) &
                  - fac * jy(ix, iy, iz)

              ez(ix, iy, iz) = ez(ix, iy, iz) &
                  + cx1 * (by(ix  , iy  , iz  ) - by(ix-1, iy  , iz  )) &
                  + cx2 * (by(ix+1, iy  , iz  ) - by(ix-2, iy  , iz  )) &
                  - cy1 * (bx(ix  , iy  , iz  ) - bx(ix  , iy-1, iz  )) &
                  - cy2 * (bx(ix  , iy+1, iz  ) - bx(ix  , iy-2, iz  )) &
                  - fac * jz(ix, iy, iz)
            ENDDO
          ENDDO
        ENDDO
      ELSE
        c1 = 75.0_num / 64.0_num
        c2 = -25.0_num / 384.0_num
        c3 = 3.0_num / 640.0_num

        DO iz = 1, nz
          cz1 = c1 * cnz
          cz2 = c2 * cnz
          cz3 = c3 * cnz
          DO iy = 1, ny
            cy1 = c1 * cny
            cy2 = c2 * cny
            cy3 = c3 * cny
            DO ix = 1, nx
              cx1 = c1 * cnx
              cx2 = c2 * cnx
              cx3 = c3 * cnx

              ex(ix, iy, iz) = ex(ix, iy, iz) &
                  + cy1 * (bz(ix  , iy  , iz  ) - bz(ix  , iy-1, iz  )) &
                  + cy2 * (bz(ix  , iy+1, iz  ) - bz(ix  , iy-2, iz  )) &
                  + cy3 * (bz(ix  , iy+2, iz  ) - bz(ix  , iy-3, iz  )) &
                  - cz1 * (by(ix  , iy  , iz  ) - by(ix  , iy  , iz-1)) &
                  - cz2 * (by(ix  , iy  , iz+1) - by(ix  , iy  , iz-2)) &
                  - cz3 * (by(ix  , iy  , iz+2) - by(ix  , iy  , iz-3)) &
                  - fac * jx(ix, iy, iz)

              ey(ix, iy, iz) = ey(ix, iy, iz) &
                  + cz1 * (bx(ix  , iy  , iz  ) - bx(ix  , iy  , iz-1)) &
                  + cz2 * (bx(ix  , iy  , iz+1) - bx(ix  , iy  , iz-2)) &
                  + cz3 * (bx(ix  , iy  , iz+2) - bx(ix  , iy  , iz-3)) &
                  - cx1 * (bz(ix  , iy  , iz  ) - bz(ix-1, iy  , iz  )) &
                  - cx2 * (bz(ix+1, iy  , iz  ) - bz(ix-2, iy  , iz  )) &
                  - cx3 * (bz(ix+2, iy  , iz  ) - bz(ix-3, iy  , iz  )) &
                  - fac * jy(ix, iy, iz)

              ez(ix, iy, iz) = ez(ix, iy, iz) &
                  + cx1 * (by(ix  , iy  , iz  ) - by(ix-1, iy  , iz  )) &
                  + cx2 * (by(ix+1, iy  , iz  ) - by(ix-2, iy  , iz  )) &
                  + cx3 * (by(ix+2, iy  , iz  ) - by(ix-3, iy  , iz  )) &
                  - cy1 * (bx(ix  , iy  , iz  ) - bx(ix  , iy-1, iz  )) &
                  - cy2 * (bx(ix  , iy+1, iz  ) - bx(ix  , iy-2, iz  )) &
                  - cy3 * (bx(ix  , iy+2, iz  ) - bx(ix  , iy-3, iz  )) &
                  - fac * jz(ix, iy, iz)
            ENDDO
          ENDDO
        ENDDO
      ENDIF
    ENDIF

  END SUBROUTINE update_e_field



  SUBROUTINE update_b_field

    INTEGER :: ix, iy, iz
    REAL(num) :: cpml_x, cpml_y, cpml_z
    REAL(num) :: c1, c2, c3
    REAL(num) :: cx1, cx2, cx3
    REAL(num) :: cy1, cy2, cy3
    REAL(num) :: cz1, cz2, cz3

    IF (cpml_boundaries) THEN
      cpml_x = hdtx
      cpml_y = hdty
      cpml_z = hdtz

      IF (field_order == 2) THEN
        DO iz = 1, nz
          cpml_z = hdtz / cpml_kappa_bz(iz)
          DO iy = 1, ny
            cpml_y = hdty / cpml_kappa_by(iy)
            DO ix = 1, nx
              cpml_x = hdtx / cpml_kappa_bx(ix)

              bx(ix, iy, iz) = bx(ix, iy, iz) &
                  - cpml_y * (ez(ix  , iy+1, iz  ) - ez(ix  , iy  , iz  )) &
                  + cpml_z * (ey(ix  , iy  , iz+1) - ey(ix  , iy  , iz  ))

              by(ix, iy, iz) = by(ix, iy, iz) &
                  - cpml_z * (ex(ix  , iy  , iz+1) - ex(ix  , iy  , iz  )) &
                  + cpml_x * (ez(ix+1, iy  , iz  ) - ez(ix  , iy  , iz  ))

              bz(ix, iy, iz) = bz(ix, iy, iz) &
                  - cpml_x * (ey(ix+1, iy  , iz  ) - ey(ix  , iy  , iz  )) &
                  + cpml_y * (ex(ix  , iy+1, iz  ) - ex(ix  , iy  , iz  ))
            ENDDO
          ENDDO
        ENDDO
      ELSE IF (field_order == 4) THEN
        c1 = 9.0_num / 8.0_num
        c2 = -1.0_num / 24.0_num

        DO iz = 1, nz
          cpml_z = hdtz / cpml_kappa_bz(iz)
          cz1 = c1 * cpml_z
          cz2 = c2 * cpml_z
          DO iy = 1, ny
            cpml_y = hdty / cpml_kappa_by(iy)
            cy1 = c1 * cpml_y
            cy2 = c2 * cpml_y
            DO ix = 1, nx
              cpml_x = hdtx / cpml_kappa_bx(ix)
              cx1 = c1 * cpml_x
              cx2 = c2 * cpml_x

              bx(ix, iy, iz) = bx(ix, iy, iz) &
                  - cy1 * (ez(ix  , iy+1, iz  ) - ez(ix  , iy  , iz  )) &
                  - cy2 * (ez(ix  , iy+2, iz  ) - ez(ix  , iy-1, iz  )) &
                  + cz1 * (ey(ix  , iy  , iz+1) - ey(ix  , iy  , iz  )) &
                  + cz2 * (ey(ix  , iy  , iz+2) - ey(ix  , iy  , iz-1))

              by(ix, iy, iz) = by(ix, iy, iz) &
                  - cz1 * (ex(ix  , iy  , iz+1) - ex(ix  , iy  , iz  )) &
                  - cz2 * (ex(ix  , iy  , iz+2) - ex(ix  , iy  , iz-1)) &
                  + cx1 * (ez(ix+1, iy  , iz  ) - ez(ix  , iy  , iz  )) &
                  + cx2 * (ez(ix+2, iy  , iz  ) - ez(ix-1, iy  , iz  ))

              bz(ix, iy, iz) = bz(ix, iy, iz) &
                  - cx1 * (ey(ix+1, iy  , iz  ) - ey(ix  , iy  , iz  )) &
                  - cx2 * (ey(ix+2, iy  , iz  ) - ey(ix-1, iy  , iz  )) &
                  + cy1 * (ex(ix  , iy+1, iz  ) - ex(ix  , iy  , iz  )) &
                  + cy2 * (ex(ix  , iy+2, iz  ) - ex(ix  , iy-1, iz  ))
            ENDDO
          ENDDO
        ENDDO
      ELSE
        c1 = 75.0_num / 64.0_num
        c2 = -25.0_num / 384.0_num
        c3 = 3.0_num / 640.0_num

        DO iz = 1, nz
          cpml_z = hdtz / cpml_kappa_bz(iz)
          cz1 = c1 * cpml_z
          cz2 = c2 * cpml_z
          cz3 = c3 * cpml_z
          DO iy = 1, ny
            cpml_y = hdty / cpml_kappa_by(iy)
            cy1 = c1 * cpml_y
            cy2 = c2 * cpml_y
            cy3 = c3 * cpml_y
            DO ix = 1, nx
              cpml_x = hdtx / cpml_kappa_bx(ix)
              cx1 = c1 * cpml_x
              cx2 = c2 * cpml_x
              cx3 = c3 * cpml_x

              bx(ix, iy, iz) = bx(ix, iy, iz) &
                  - cy1 * (ez(ix  , iy+1, iz  ) - ez(ix  , iy  , iz  )) &
                  - cy2 * (ez(ix  , iy+2, iz  ) - ez(ix  , iy-1, iz  )) &
                  - cy3 * (ez(ix  , iy+3, iz  ) - ez(ix  , iy-2, iz  )) &
                  + cz1 * (ey(ix  , iy  , iz+1) - ey(ix  , iy  , iz  )) &
                  + cz2 * (ey(ix  , iy  , iz+2) - ey(ix  , iy  , iz-1)) &
                  + cz3 * (ey(ix  , iy  , iz+3) - ey(ix  , iy  , iz-2))

              by(ix, iy, iz) = by(ix, iy, iz) &
                  - cz1 * (ex(ix  , iy  , iz+1) - ex(ix  , iy  , iz  )) &
                  - cz2 * (ex(ix  , iy  , iz+2) - ex(ix  , iy  , iz-1)) &
                  - cz3 * (ex(ix  , iy  , iz+3) - ex(ix  , iy  , iz-2)) &
                  + cx1 * (ez(ix+1, iy  , iz  ) - ez(ix  , iy  , iz  )) &
                  + cx2 * (ez(ix+2, iy  , iz  ) - ez(ix-1, iy  , iz  )) &
                  + cx3 * (ez(ix+3, iy  , iz  ) - ez(ix-2, iy  , iz  ))

              bz(ix, iy, iz) = bz(ix, iy, iz) &
                  - cx1 * (ey(ix+1, iy  , iz  ) - ey(ix  , iy  , iz  )) &
                  - cx2 * (ey(ix+2, iy  , iz  ) - ey(ix-1, iy  , iz  )) &
                  - cx3 * (ey(ix+3, iy  , iz  ) - ey(ix-2, iy  , iz  )) &
                  + cy1 * (ex(ix  , iy+1, iz  ) - ex(ix  , iy  , iz  )) &
                  + cy2 * (ex(ix  , iy+2, iz  ) - ex(ix  , iy-1, iz  )) &
                  + cy3 * (ex(ix  , iy+3, iz  ) - ex(ix  , iy-2, iz  ))
            ENDDO
          ENDDO
        ENDDO
      ENDIF

      CALL cpml_advance_b_currents(hdt)
    ELSE
      IF (field_order == 2) THEN
        DO iz = 1, nz
          DO iy = 1, ny
            DO ix = 1, nx
              bx(ix, iy, iz) = bx(ix, iy, iz) &
                  - hdty * (ez(ix  , iy+1, iz  ) - ez(ix  , iy  , iz  )) &
                  + hdtz * (ey(ix  , iy  , iz+1) - ey(ix  , iy  , iz  ))

              by(ix, iy, iz) = by(ix, iy, iz) &
                  - hdtz * (ex(ix  , iy  , iz+1) - ex(ix  , iy  , iz  )) &
                  + hdtx * (ez(ix+1, iy  , iz  ) - ez(ix  , iy  , iz  ))

              bz(ix, iy, iz) = bz(ix, iy, iz) &
                  - hdtx * (ey(ix+1, iy  , iz  ) - ey(ix  , iy  , iz  )) &
                  + hdty * (ex(ix  , iy+1, iz  ) - ex(ix  , iy  , iz  ))
            ENDDO
          ENDDO
        ENDDO
      ELSE IF (field_order == 4) THEN
        c1 = 9.0_num / 8.0_num
        c2 = -1.0_num / 24.0_num

        DO iz = 1, nz
          cz1 = c1 * hdtz
          cz2 = c2 * hdtz
          DO iy = 1, ny
            cy1 = c1 * hdty
            cy2 = c2 * hdty
            DO ix = 1, nx
              cx1 = c1 * hdtx
              cx2 = c2 * hdtx

              bx(ix, iy, iz) = bx(ix, iy, iz) &
                  - cy1 * (ez(ix  , iy+1, iz  ) - ez(ix  , iy  , iz  )) &
                  - cy2 * (ez(ix  , iy+2, iz  ) - ez(ix  , iy-1, iz  )) &
                  + cz1 * (ey(ix  , iy  , iz+1) - ey(ix  , iy  , iz  )) &
                  + cz2 * (ey(ix  , iy  , iz+2) - ey(ix  , iy  , iz-1))

              by(ix, iy, iz) = by(ix, iy, iz) &
                  - cz1 * (ex(ix  , iy  , iz+1) - ex(ix  , iy  , iz  )) &
                  - cz2 * (ex(ix  , iy  , iz+2) - ex(ix  , iy  , iz-1)) &
                  + cx1 * (ez(ix+1, iy  , iz  ) - ez(ix  , iy  , iz  )) &
                  + cx2 * (ez(ix+2, iy  , iz  ) - ez(ix-1, iy  , iz  ))

              bz(ix, iy, iz) = bz(ix, iy, iz) &
                  - cx1 * (ey(ix+1, iy  , iz  ) - ey(ix  , iy  , iz  )) &
                  - cx2 * (ey(ix+2, iy  , iz  ) - ey(ix-1, iy  , iz  )) &
                  + cy1 * (ex(ix  , iy+1, iz  ) - ex(ix  , iy  , iz  )) &
                  + cy2 * (ex(ix  , iy+2, iz  ) - ex(ix  , iy-1, iz  ))
            ENDDO
          ENDDO
        ENDDO
      ELSE
        c1 = 75.0_num / 64.0_num
        c2 = -25.0_num / 384.0_num
        c3 = 3.0_num / 640.0_num

        DO iz = 1, nz
          cz1 = c1 * hdtz
          cz2 = c2 * hdtz
          cz3 = c3 * hdtz
          DO iy = 1, ny
            cy1 = c1 * hdty
            cy2 = c2 * hdty
            cy3 = c3 * hdty
            DO ix = 1, nx
              cx1 = c1 * hdtx
              cx2 = c2 * hdtx
              cx3 = c3 * hdtx

              bx(ix, iy, iz) = bx(ix, iy, iz) &
                  - cy1 * (ez(ix  , iy+1, iz  ) - ez(ix  , iy  , iz  )) &
                  - cy2 * (ez(ix  , iy+2, iz  ) - ez(ix  , iy-1, iz  )) &
                  - cy3 * (ez(ix  , iy+3, iz  ) - ez(ix  , iy-2, iz  )) &
                  + cz1 * (ey(ix  , iy  , iz+1) - ey(ix  , iy  , iz  )) &
                  + cz2 * (ey(ix  , iy  , iz+2) - ey(ix  , iy  , iz-1)) &
                  + cz3 * (ey(ix  , iy  , iz+3) - ey(ix  , iy  , iz-2))

              by(ix, iy, iz) = by(ix, iy, iz) &
                  - cz1 * (ex(ix  , iy  , iz+1) - ex(ix  , iy  , iz  )) &
                  - cz2 * (ex(ix  , iy  , iz+2) - ex(ix  , iy  , iz-1)) &
                  - cz3 * (ex(ix  , iy  , iz+3) - ex(ix  , iy  , iz-2)) &
                  + cx1 * (ez(ix+1, iy  , iz  ) - ez(ix  , iy  , iz  )) &
                  + cx2 * (ez(ix+2, iy  , iz  ) - ez(ix-1, iy  , iz  )) &
                  + cx3 * (ez(ix+3, iy  , iz  ) - ez(ix-2, iy  , iz  ))

              bz(ix, iy, iz) = bz(ix, iy, iz) &
                  - cx1 * (ey(ix+1, iy  , iz  ) - ey(ix  , iy  , iz  )) &
                  - cx2 * (ey(ix+2, iy  , iz  ) - ey(ix-1, iy  , iz  )) &
                  - cx3 * (ey(ix+3, iy  , iz  ) - ey(ix-2, iy  , iz  )) &
                  + cy1 * (ex(ix  , iy+1, iz  ) - ex(ix  , iy  , iz  )) &
                  + cy2 * (ex(ix  , iy+2, iz  ) - ex(ix  , iy-1, iz  )) &
                  + cy3 * (ex(ix  , iy+3, iz  ) - ex(ix  , iy-2, iz  ))
            ENDDO
          ENDDO
        ENDDO
      ENDIF
    ENDIF

  END SUBROUTINE update_b_field



  SUBROUTINE update_eb_fields_half

    hdt  = 0.5_num * dt
    hdtx = hdt / dx
    hdty = hdt / dy
    hdtz = hdt / dz

    cnx = hdtx * c**2
    cny = hdty * c**2
    cnz = hdtz * c**2

    fac = hdt / epsilon0

    ! Update E field to t+dt/2
    CALL update_e_field

    ! Now have E(t+dt/2), do boundary conditions on E
    CALL efield_bcs

    ! Update B field to t+dt/2 using E(t+dt/2)
    CALL update_b_field

    ! Now have B field at t+dt/2. Do boundary conditions on B
    CALL bfield_bcs(.TRUE.)

    ! Now have E&B fields at t = t+dt/2
    ! Move to particle pusher

  END SUBROUTINE update_eb_fields_half



  SUBROUTINE update_eb_fields_final

    hdt  = 0.5_num * dt
    hdtx = hdt / dx
    hdty = hdt / dy
    hdtz = hdt / dz

    cnx = hdtx * c**2
    cny = hdty * c**2
    cnz = hdtz * c**2

    fac = hdt / epsilon0

    CALL update_b_field

    CALL bfield_final_bcs

    CALL update_e_field

    CALL efield_bcs

  END SUBROUTINE update_eb_fields_final

#ifdef NONLINEAR_OPTICS
  SUBROUTINE update_medium_polarisation()
    !----------------------------------------------------------------------------
    ! Nonlinear optics - Written by C. Varin and R. Emms
    !----------------------------------------------------------------------------
    REAL(num) :: omega_0
    REAL(num) :: a, b, g
    INTEGER :: start_rug
    INTEGER :: ix,iy,iz

    omega_0 = 2.0*pi*c/resonance
    a = (1.0 - gamma*dt/2.0)/(1.0 + gamma*dt/2.0)
    b = omega_0**2*dt/(1.0 + gamma*dt/2.0)

    intensity = ex**2 + ey**2 + ez**2
    eps0chi = epsilon0*(chi1 + chi2*sqrt(intensity) + chi3*intensity)*medium_mask
    ! Note, medium_mask is an array of values from 0.0 to 1.0 where
    ! 0.0 = no medium
    ! 1.0 = medium at full density

    !----------------------------------------------------------------------------
    ! Calculate medium current and polarization densities
    !----------------------------------------------------------------------------

    jx_nlo = a*jx_nlo + b*( eps0chi*ex - px_nlo )
    px_nlo = px_nlo + dt*jx_nlo
    
    jy_nlo = a*jy_nlo + b*( eps0chi*ey - py_nlo )
    py_nlo = py_nlo + dt*jy_nlo
    
    jz_nlo = a*jz_nlo + b*( eps0chi*ez - pz_nlo )
    pz_nlo = pz_nlo + dt*jz_nlo

    !----------------------------------------------------------------------------
    ! Drude Model - Written by C. Varin and J.-L. Deziel
    !----------------------------------------------------------------------------
    ! If Drude model is static, omega_p is constant. If Drude model is dynamic,
    ! we calculate omega_p with electron_density_Drude.

    IF (use_Drude) THEN
      IF (Drude_static) THEN
        omega_p = omega_p_static*medium_mask
      ELSE
        ! this subroutine updates the electron density using the
        ! provided ionization rates (keldysh model)
        IF (use_dynamic_gamma_drude) THEN
          CALL update_gamma_drude
        ENDIF
        CALL update_electron_density_drude
      ENDIF
      reduced_mass = 1.0/(1.0/effective_mass_electron + 1.0/effective_mass_hole)
      omega_p = (electron_density_Drude*q0**2.0/(reduced_mass*m0*epsilon0))**0.5 &
        *medium_mask

      DO iz = 1, nz
        DO iy = 1, ny
          DO ix = 1, nx
            g = dynamic_gamma_drude(ix,iy,iz)
            jx_D(ix,iy,iz) = (1.0 - g*dt/2.0)/(1.0 + g*dt/2.0)*jx_D(ix,iy,iz) &
              + epsilon0*dt/(1.0 + g*dt/2.0)*omega_p(ix,iy,iz)**2.0*ex(ix,iy,iz)
            jy_D(ix,iy,iz) = (1.0 - g*dt/2.0)/(1.0 + g*dt/2.0)*jy_D(ix,iy,iz) &
              + epsilon0*dt/(1.0 + g*dt/2.0)*omega_p(ix,iy,iz)**2.0*ey(ix,iy,iz)
            jz_D(ix,iy,iz) = (1.0 - g*dt/2.0)/(1.0 + g*dt/2.0)*jz_D(ix,iy,iz) &
              + epsilon0*dt/(1.0 + g*dt/2.0)*omega_p(ix,iy,iz)**2.0*ez(ix,iy,iz)
          ENDDO
        ENDDO
      ENDDO
    ENDIF

    !----------------------------------------------------------------------------
    ! Add contribution to the global (pic) current density
    !---------------------------------------------------------------------------- 
    jz(1-ng:nx+ng, 1-ng:ny+ng, 1-ng:nz+ng) = &
                 jz(1-ng:nx+ng, 1-ng:ny+ng, 1-ng:nz+ng) &
                  + jz_nlo(1-ng:nx+ng, 1-ng:ny+ng, 1-ng:nz+ng) &
                  + jz_D(1-ng:nx+ng, 1-ng:ny+ng, 1-ng:nz+ng) &
                  + jz_mpi(1-ng:nx+ng, 1-ng:ny+ng, 1-ng:nz+ng)

    jy(1-ng:nx+ng, 1-ng:ny+ng, 1-ng:nz+ng) = &
                  jy(1-ng:nx+ng, 1-ng:ny+ng, 1-ng:nz+ng) &
                  + jy_nlo(1-ng:nx+ng, 1-ng:ny+ng, 1-ng:nz+ng) &
                  + jy_D(1-ng:nx+ng, 1-ng:ny+ng, 1-ng:nz+ng) &
                  + jy_mpi(1-ng:nx+ng, 1-ng:ny+ng, 1-ng:nz+ng)

    jx(1-ng:nx+ng, 1-ng:ny+ng, 1-ng:nz+ng) = &
                 jx(1-ng:nx+ng, 1-ng:ny+ng, 1-ng:nz+ng) &
                   + jx_nlo(1-ng:nx+ng, 1-ng:ny+ng, 1-ng:nz+ng) &
                   + jx_D(1-ng:nx+ng, 1-ng:ny+ng, 1-ng:nz+ng) &
                   + jx_mpi(1-ng:nx+ng, 1-ng:ny+ng, 1-ng:nz+ng)

  END SUBROUTINE update_medium_polarisation


  SUBROUTINE update_electron_density_drude()

    REAL(num) :: gamma_keldysh, E_mag,omega, rate
    INTEGER :: ix,iy,iz,ig,ik
    TYPE(laser_block), POINTER :: current_laser
    REAL(num) :: mre_one_photon_absorption
    REAL(num) :: saturation,rate_col

    current_laser => laser_x_min

    DO iz = 1, nz
      DO iy = 1, ny
        DO ix = 1, nx

          ! Calculate the e field at every point
          E_mag = SQRT(ex(ix,iy,iz)**2 &
            +ey(ix,iy,iz)**2+ez(ix,iy,iz)**2)

          ! Then the keldysh parameter
          gamma_keldysh = SQRT(m0*bandGap_drude*reduced_mass) &
            *current_laser%omega/q0/E_mag

          ! Then the ionization rate from the ionRate_table
          ig = (gamma_keldysh - 0.01)/0.0005
          IF (gamma_keldysh .GT. 50) THEN
            rate = 0.0_num
          ELSE
            rate  = ionRate_table(1,ig) + (ionRate_table(1,ig+1) - &
              ionRate_table(1,ig))*(gamma_keldysh-ig*0.0005 - &
                0.01)/0.0005 * medium_mask(ix,iy,iz)
          ENDIF
          ! Saturation
          saturation = (atom_density - electron_density_Drude(ix,iy,iz))/atom_density
          rate = rate*saturation*atom_density

          ! We also add the collisionnal ionization rate
          ! Then electron populations are updated
          IF (use_mre_avalanche) THEN

            drude_cross_section = q0**2.0 &
              /(dynamic_gamma_drude(ix,iy,iz)*m0*reduced_mass &
                *(1.0_num+(current_laser%omega/dynamic_gamma_drude(ix,iy,iz))**2.0))

            ! Effective alpha value for the SRE limit
            ! mre_avalanche_factor = drude_cross_section &
            !  *amp_laser**2.0_num /(0.69315_num*mre_critical_energy)
            ! Alternative way to get alpha is to use 
            ! the electron-neutral atom collision rate
            IF (use_dynamic_gamma_drude) THEN
              ! We don't add the saturation factor, because it comes within g_en
              mre_avalanche_factor = g_en(ix,iy,iz)
            ELSE
              ! If gamma_drude is constant, we add the saturation factor
              mre_avalanche_factor = dynamic_gamma_drude(ix,iy,iz)*saturation
            ENDIF
            
            ! Inverse bremsstrahlung rate
            mre_one_photon_absorption = drude_cross_section &
              *E_mag**2.0_num/(0.69315_num*mre_critical_energy) &
              /(2.0_num**(1.0_num/mre_nb_levels)-1.0_num)
            ! We take the harmonic mean with the (g_en+g_ei) collision rate
            ! because an electron needs to be in the vinicity of a massive
            ! particle in order to asborb vis IB.
            IF (use_dynamic_gamma_drude) THEN
              mre_one_photon_absorption = 1.0_num/(1.0_num/ &
                (g_ei(ix,iy,iz)+g_en(ix,iy,iz)) &
                +1.0_num/mre_one_photon_absorption)
            ELSE
              mre_one_photon_absorption = 1.0_num/(1.0_num/(dynamic_gamma_drude(ix,iy,iz)) &
                +1.0_num/mre_one_photon_absorption)
            ENDIF

            ! track sfi contribution
            electron_density_sfi(ix,iy,iz) = electron_density_sfi(ix,iy,iz) + dt*rate
            ! track collisionnal contribution
            electron_density_col(ix,iy,iz) = electron_density_col(ix,iy,iz) + dt &
              *mre_avalanche_factor*mre_rho(mre_nb_levels,ix,iy,iz)

            ! Populations evolution in the CB
            mre_rho(0,ix,iy,iz) = mre_rho(0,ix,iy,iz) + dt &
              *(rate + 2.0_num*mre_avalanche_factor*mre_rho(mre_nb_levels,ix,iy,iz) &
                -(mre_one_photon_absorption + recombination_rate)*mre_rho(0,ix,iy,iz))

            DO ik=1,mre_nb_levels-1
              mre_rho(ik,ix,iy,iz) = mre_rho(ik,ix,iy,iz) + dt &
                *(mre_one_photon_absorption*mre_rho(ik-1,ix,iy,iz) &
                  -(mre_one_photon_absorption + recombination_rate) &
                    *mre_rho(ik,ix,iy,iz))
            ENDDO

            mre_rho(mre_nb_levels,ix,iy,iz) = mre_rho(mre_nb_levels,ix,iy,iz) + dt &
              *(mre_one_photon_absorption*mre_rho(mre_nb_levels-1,ix,iy,iz) &
                -(recombination_rate + mre_avalanche_factor) &
                  *mre_rho(mre_nb_levels,ix,iy,iz))

            ! !Plasma drift (still dirty and unstable)
            ! electron_density_Drude(ix,iy,iz) = SUM(mre_rho(:,ix,iy,iz))
            ! DO ik=0,mre_nb_levels
            !   mre_rho(ik,ix,iy,iz) = mre_rho(ik,ix,iy,iz) - dt/q0 &
            !     *((jx_D(ix+1,iy,iz)-jx_D(ix-1,iy,iz))/(2.0*dx) &
            !       +(jy_D(ix,iy+1,iz)-jy_D(ix,iy-1,iz))/(2.0*dy) &
            !       +(jz_D(ix,iy,iz+1)-jz_D(ix,iy,iz-1))/(2.0*dz)) &
            !       *(mre_rho(ik,ix,iy,iz)/electron_density_Drude(ix,iy,iz))
            ! ENDDO

            electron_density_Drude(ix,iy,iz) = SUM(mre_rho(:,ix,iy,iz))

            electron_temperature(ix,iy,iz) = room_temperature*medium_mask(ix,iy,iz)
            IF (electron_density_Drude(ix,iy,iz) .GT. 0.0_num) THEN
              DO ik=1,mre_nb_levels
                electron_temperature(ix,iy,iz) = electron_temperature(ix,iy,iz) &
                     + mre_rho(ik,ix,iy,iz)*ik * 2.0_num*h_bar*current_laser%omega &
                  /(3.0_num*kb*electron_density_Drude(ix,iy,iz))
              ENDDO
            ENDIF
          
          ELSE !SRE
            ! track sfi contribution
            electron_density_sfi(ix,iy,iz) = electron_density_sfi(ix,iy,iz) + dt*rate

            rate_col = avalanche_factor*intensity(ix,iy,iz) &
              *SQRT(1.0_num+chi1)*epsilon0*c/2.0_num &
                *electron_density_Drude(ix,iy,iz)*saturation

            ! track collisionnal contribution
            electron_density_col(ix,iy,iz) = electron_density_col(ix,iy,iz) + dt*rate_col

            ! Recombination rate
            electron_density_Drude(ix,iy,iz) = &
              electron_density_Drude(ix,iy,iz) - dt*recombination_rate &
                *electron_density_Drude(ix,iy,iz)

            ! Finally, the growth of electron_density_drude
            electron_density_Drude(ix,iy,iz) = &
              electron_density_Drude(ix,iy,iz)+dt*(rate+rate_col)
          ENDIF

          ! We also calculate the losses due to mpi absorption
          IF (E_mag .GT. 0.0_num) THEN
            jx_mpi(ix,iy,iz) = bandGap_drude*rate/(E_mag**2.0_num)*ex(ix,iy,iz)
            jy_mpi(ix,iy,iz) = bandGap_drude*rate/(E_mag**2.0_num)*ey(ix,iy,iz)
            jz_mpi(ix,iy,iz) = bandGap_drude*rate/(E_mag**2.0_num)*ez(ix,iy,iz)
          ELSE
            jx_mpi(ix,iy,iz) = 0.0_num
            jy_mpi(ix,iy,iz) = 0.0_num
            jz_mpi(ix,iy,iz) = 0.0_num
          ENDIF
        ENDDO
      ENDDO
    ENDDO

  END SUBROUTINE update_electron_density_drude

  SUBROUTINE update_gamma_drude()

    INTEGER :: ix,iy,iz,ik
    TYPE(laser_block), POINTER :: current_laser
    REAL(num) :: a,b,c,d,e
    REAL(num) :: Coulomb_log,mass_ratio
    REAL(num) :: d2,d3

    a = 4.0_num*pi*epsilon0/(q0**2.0_num) &
      *(6.0_num/(reduced_mass*m0))**0.5_num *kb**1.5_num
    b = q0**2.0*kb*room_temperature*reduced_mass*m0 &
      /(4.0_num*pi*epsilon0*h_bar**3.0_num) &
        *(3.0_num*pi**2.0_num)**(-1.0_num/3.0_num)
    c = atom_cross_section*(3.0_num*kb/(reduced_mass*m0))**0.5_num
    ! d = q0**4.0_num *(3.0_num*kb)**(-1.5_num) &
    !   /(4.0_num*pi*epsilon0**2.0_num*(reduced_mass*m0)**0.5_num)
    d2 = (3.0_num*kb)**2.0_num/(q0**4.0_num*(4.0_num*pi/3.0_num)**(4.0_num/3.0_num))
    d3 = (1.0_num/3.0_num)*(reduced_mass*m0/(6.0_num*pi*kb))**1.5_num
    e = (3.0_num*kb/(reduced_mass*m0))**0.5_num &
      *(4.0_num*pi*atom_density/3.0_num)**(1.0_num/3.0_num)

    DO iz = 1, nz
      DO iy = 1, ny
        DO ix = 1, nx
          ! electron-electron scattering, see [christensen2009]
          g_ee(ix,iy,iz) = a * electron_temperature(ix,iy,iz)**1.5_num
          ! electron-phonon scattering, see [bailling2013]
          g_ep(ix,iy,iz) = b * electron_density_Drude(ix,iy,iz)**(-1.0_num/3.0_num)
          ! electron-neutral atom scattering, see [bailling2013]
          g_en(ix,iy,iz) = c * (atom_density - electron_density_Drude(ix,iy,iz)) &
            * electron_temperature(ix,iy,iz)**0.5_num
          ! electron-ion scattering, see [gallais2015]
          ! g_ei(ix,iy,iz) = d * electron_temperature(ix,iy,iz)**(-1.5_num) &
          !   *(electron_density_Drude(ix,iy,iz)+1.0_num)
          ! electron-ion scattering, spitzel model
          Coulomb_log = 0.5_num*log(1.0_num+(d2*electron_temperature(ix,iy,iz)**2.0_num &
            /((electron_density_Drude(ix,iy,iz)+1.0_num)**(4.0_num/3.0_num))))
          g_ei(ix,iy,iz) = d3*omega_p(ix,iy,iz)**4.0_num*Coulomb_log &
            /(1.0_num+electron_density_Drude(ix,iy,iz)) &
            /(electron_temperature(ix,iy,iz)+1.0_num)**1.5_num
          ! maximal scattering, see [bailling2013]
          g_max(ix,iy,iz) = e * electron_temperature(ix,iy,iz)**0.5_num

          ! harmonic mean of the scattering rates, see [chimier2011]
          dynamic_gamma_drude(ix,iy,iz) = 1.0_num+ &
            ((g_ee(ix,iy,iz) + g_ep(ix,iy,iz))**(-1.0_num) &
            + (g_en(ix,iy,iz) + g_ei(ix,iy,iz))**(-1.0_num) &
            + g_max(ix,iy,iz)**(-1.0_num))**(-1.0_num)
        ENDDO
      ENDDO
    ENDDO

  END SUBROUTINE update_gamma_drude

#endif

END MODULE fields
