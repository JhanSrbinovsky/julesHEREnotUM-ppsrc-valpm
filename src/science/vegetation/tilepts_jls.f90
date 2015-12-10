
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! SUBROUTINE tilepts
!
! Purpose:
! Counts the number of points containing each surface type and creates
! a TILE_INDEX array specifying the location of these points on the land
! grid.

! Subroutine Interface:

SUBROUTINE tilepts(land_pts,frac,tile_pts,tile_index)

USE switches, ONLY :                                              &
 all_tiles

USE nstypes, ONLY :                                               &
!      imported scalars with intent(in)
  ntype

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
IMPLICIT NONE

INTEGER, INTENT(IN) :: land_pts  !num land points to process

REAL, INTENT(IN) :: frac(land_pts,ntype)  !fractions of surface types

INTEGER, INTENT(OUT) ::  tile_pts(ntype)                      &
                                          ! Number of land points which
                                          ! include the nth surface type
,                        tile_index(land_pts,ntype)
                                          ! Indices of land points which
                                          ! include the nth surface type

INTEGER :: n,l,c   !local counters: type, land pts, count

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook('TILEPTS',zhook_in,zhook_handle)

!-----------------------------------------------------------------------
! Create the TILE_INDEX array of land points with each surface type
!-----------------------------------------------------------------------
tile_pts (:) = 0
tile_index(:,:) = 0

DO n=1,ntype
  c=0
!CDIR NODEP
  DO l=1,land_pts
    IF ( ((all_tiles == 0).AND.                                  &
          (frac(l,n) >  0.0))                                    &
        .OR.                                                     &
         ((all_tiles == 1).AND.                                  &
          ( (n .lt. ntype .AND. frac(l,ntype) .lt. 0.5)          &
           .OR.                                                  &
            (n .eq. ntype .AND. frac(l,ntype) .gt. 0.5))) )THEN
      c = c + 1
      tile_index(c,n) = l
    END IF
  END DO
  tile_pts(n) = c
END DO

IF (lhook) CALL dr_hook('TILEPTS',zhook_out,zhook_handle)
RETURN
END SUBROUTINE tilepts
