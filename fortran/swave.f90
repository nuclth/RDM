
MODULE swave
  IMPLICIT NONE

  TYPE h1
    REAL(8), ALLOCATABLE :: one_body (:,:)
  ENDTYPE h1

  TYPE h2
    REAL(8), ALLOCATABLE :: two_body (:,:,:,:)
  ENDTYPE h2

  CONTAINS

SUBROUTINE populate_h1

SUBROUTINE populate_h2
