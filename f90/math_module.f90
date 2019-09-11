MODULE math_module

 IMPLICIT NONE

 CONTAINS

 SUBROUTINE linear_interpolation(N, N0, x, y, x0, y0)
 !// Linear interpolation program. Provide number N of x and y values, 
 !// the new x-values x0[N0] values and an empty y0[N0]. y0 will then 
 !// recieve the interpolated y-values. // JSV
  USE precision_type
  IMPLICIT NONE
  integer   :: N, N0, i, j, k, k2
  real(dp)  :: x(N), y(N), x0(N0), y0(N0), dist, a, b
 !// LOOP OVER ELEMENTS IN X0 ARRAY
 
  DO j=1,N0
     !// Locate x0(i)
     k = 1
     dist = abs(x(1)-x0(j))
     DO i=1,N
        IF (abs(x(i) - x0(j)) < dist) THEN
           k = i
           dist = abs(x(i) - x0(j))
        END IF
     ENDDO
     k2 = k
     IF (x0(j) < x(k)) THEN
        k2 = k-1
     END IF
     IF (k2 == -1) THEN 
        k2 = 1
        write(*,*) "Error, interpolation x-value out of range"
     END IF
     IF (x0(j) > x(N)) THEN
        k2 = N-1
        write(*,*) "Error, interpolation x-value out of range"
        stop
     END IF
     a = (y(k2+1)- y(k2)) / (x(k2+1) - x(k2))
     b = -a*x(k2) + y(k2)
     y0(j) = a * x0(j) + b
  ENDDO
 END SUBROUTINE linear_interpolation
 
 SUBROUTINE cubic_interpolation_2d(Nx, Ny, Nx0, Ny0, x, y, x0, y0, matrix, matrix0, silentFLAG)
 USE precision_type
 IMPLICIT NONE
 !//------------------------------------------------------------------------//
 !// Code to perform linear interpolation from one matrix of a regular grid 
 !// to another one. Note that "matrix" and "matrix0" must be defined in 
 !// mother program as (not sure if N is rows or columns, should be checked 
 !// if not symmetric):
 !//
 !//    old_matrix[N][N];    //matrix of data to be interpolated from
 !//    old_matrix0[N0][N0]; // matrix representing the new grid
 !//    double *matrix[N], *matrix0[N0];
 !//    for(int i=0;i<N-1;i++)  matrix[i]  = old_matrix[i]; 
 !//    for(int i=0;i<N0-1;i++) matrix0[i] = old_matrix0[i]; 
 !//    cubic_interpolation_2d( ... , matrix, matrix0);
 !//
 !// JS
 !//------------------------------------------------------------------------//
 
  integer :: Nx, Ny, Nx0, Ny0, silentFLAG
  integer :: i, j, k, kx, ky, x0i(Nx0), y0i(Ny0)
  real(dp) :: u, t, dist
  real(dp) :: x(Nx), y(Ny), x0(Nx0), y0(Ny0), matrix(Nx,Ny), matrix0(Nx0,Ny0)
 
 ! CHECK INTERPOLATION RANGE
  IF (silentFLAG .ne. 1) THEN
     DO i=1,Nx0 !// LOOP OVER ELEMENTS IN X0 ARRAY
        IF (x0(i) < x(1)) THEN
           write(*,*) "Error, interpolation out of range"
           STOP
        END IF
        IF (x0(i) > x(Nx)) THEN
           write(*,*) "Error, interpolation out of range"
           STOP
        END IF
     ENDDO
     DO i=1,Ny0 !// LOOP OVER ELEMENTS IN Y0 ARRAY
        IF (y0(i) < y(1)) THEN 
           write(*,*) "Error, interpolation out of range"
           STOP
        END IF
        IF (y0(i) > y(Ny)) THEN
           write(*,*) "Error, interpolation out of range"
           STOP
        END IF
     ENDDO
  END IF
 
 !// LOCATE X-VALUES PLACEMENT IN GRID
  DO j=1,Nx0
     k = 1
     dist = abs(x(1)-x0(j))
     DO i=1,Nx
        IF (abs(x(i) - x0(j)) < dist) THEN
           k=i
           dist = abs(x(i) - x0(j))
        END IF
     ENDDO
     IF (x0(j) < x(k)) THEN 
        k = k-1
     END IF
     x0i(j) = k
  ENDDO 
 
 !// LOCATE Y-VALUES PLACEMENT IN GRID
  DO j=1,Ny0 
     k = 1
     dist = abs(y(1)-y0(j))
     DO i=1,Ny
        IF (abs(y(i) - y0(j)) < dist) THEN
           k=i
           dist = abs(y(i) - y0(j))
        END IF
     ENDDO
     IF (y0(j) < y(k)) THEN
        k = k-1
     END IF
     y0i(j) = k
  ENDDO
 
 !// INTERPOLATE
  DO i=1,Nx0
     DO j=1,Ny0
        kx = x0i(i)
        ky = y0i(j)
        t = (x0(i)-x(kx)) / (x(kx+1)-x(kx))
        u = (y0(j)-y(ky)) / (y(ky+1)-y(ky))
        matrix0(i,j) = (1.0-t)*(1.0-u)*matrix(kx,ky) + t*(1.0-u)*matrix(kx+1,ky) + t*u*matrix(kx+1,ky+1) + (1.0-t)*u*matrix(kx,ky+1)
     ENDDO
  ENDDO
 END SUBROUTINE CUBIC_INTERPOLATION_2D 
 
 
END MODULE math_module
