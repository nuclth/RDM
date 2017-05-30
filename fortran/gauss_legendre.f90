
MODULE gauss_legendre
  IMPLICIT NONE
  REAL(8), PARAMETER :: pi  = 3.1415926535897932d0
  REAL(8), PARAMETER :: eps = 3.0D-14

  TYPE gauss_mesh
   INTEGER npts, job
   REAL(8) :: a, b
   REAL(8), ALLOCATABLE :: x, weights
  ENDTYPE

  CONTAINS

SUBROUTINE create_mesh (gmesh)
 
  TYPE gauss_mesh :: gmesh

  INTEGER :: i, m

  REAL (8) :: t, p1, p2, p3

  INTENT(IN)  :: gmesh
  INTENT(OUT) :: gmesh

  m = (npts + 1) / 2


  DO i = 1, m
    
    t = cos (pi * (i - 0.25) / (npts + 0.5))
    t1 = 1.0

    WHILE (abs(t-t1) .GE. eps)
      p1 = 1.0
      p2 = 0.0

      DO j = 1, npts
        p3 = p2
        p2 = p1
        p1 = ((2.0 *j - 1.0) * t * p2 - (j - 1.0) * p3) / j
      ENDDO
      pp = npts * (t * p1 - p2) / (t * t - 1.);
      t1 = t;
      t = t1 - p1 / pp;
    ENDWHILE   
 

    x(i - 1) = -t;
    x(npts - i) = t;
    weights(i - 1) = 2.0 / ((1. - t * t) * pp * pp);
    weights(npts - i) = weights(i - 1);
  



  ENDDO

nt i = 1; i <= m; i++)
  {
    t = cos (pi * (i - 0.25) / (npts + 0.5));
    t1 = 1.;
    while ((fabs (t - t1)) >= eps)
    {
      p1 = 1.0;
      p2 = 0.0;
      for (int j = 1; j <= npts; j++)
      {
        p3 = p2;
        p2 = p1;
        p1 = ((2. * j - 1.) * t * p2 - (j - 1.) * p3) / j;
      }
      pp = npts * (t * p1 - p2) / (t * t - 1.);
      t1 = t;
      t = t1 - p1 / pp;
    }
    xpts[i - 1] = -t;
    xpts[npts - i] = t;
    weights[i - 1] = 2.0 / ((1. - t * t) * pp * pp);
    weights[npts - i] = weights[i - 1];
  }


  ! Rescale from [a,inf] with 50% of points in [a,b+2a]
  IF job .EQ. 2 THEN

   DO i = 1, npts
     t = 1.0 - x[i]
     x[i] = (b * x[i] + b + a + a) / t
     w[i] = w[i] * 2.0 * (a + b) / (t * t)
   ENDDO
  ENDIF

END SUBROUTINE create_mesh

END MODULE gauss_legendre  
  
