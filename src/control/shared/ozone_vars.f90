! File containing variables for ozone implementation

MODULE ozone_vars







! Ozone diagnostics
  REAL, ALLOCATABLE ::         &
    flux_o3_ft(:,:)            &
                 ! Flux of O3 to stomata (nmol O3/m2/s).
   ,fo3_ft(:,:)  ! Ozone exposure factor.

END MODULE
