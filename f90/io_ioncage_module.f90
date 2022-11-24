MODULE io_ioncage_module

 IMPLICIT NONE

 CONTAINS
 
 SUBROUTINE read_coef_matrix(filename, matrix, d, silentFLAG, Ntot, NH)
    ! Imports coagulation coefficient matrix from file "filename" 
    ! and interpolates to grid provided by "d" into "matrix".
    !
    ! INPUT:
    ! Ntot : Number of diameter nodes in simulation for aerosols
    ! filename : Name of file to read coefficients from
    ! d : diameter of the Ntot nodes in the simulation. "d" is the grid that
    !     the coefficents from "filename" is interpolated to.
    !
    ! NH : number of header lines to skip in "filename".
    !
    ! OUTPUT: 
    ! matrix : An Ntot x Ntot size matrix containing interpolated coefficents
    !
    ! REQUIREMENTS FOR "filename":
    ! "filename" may have a header, which is skipped by providing number of
    ! header lines NH. Otherwise code assumes d-header, d-array, coef-header, 
    ! and coef-matrix. There must not be lines after the coef-matrix.
    !
    ! Example for Ntot = 3, and NH = 2:
    !  
    !   >HEADER HEADER HEADER HEADER
    !   >HEADER HEADER HEADER HEADER 
    !   >D-array [m]:
    !   >        7.00000E-10         8.91788E-10         1.13612E-09
    !   >Coagulation-matrix (s m^3):
    !   >        4.81600E-16         6.95071E-16         1.39363E-15
    !   >        6.95071E-16         6.34285E-16         9.15315E-16
    !   >        1.39363E-15         9.15315E-16         8.35062E-16
    !   

    USE precision_type
    USE math_module
    USE io_general_module
    IMPLICIT NONE
    integer           :: i, silentFLAG, lunit, Np, Ntot, x, N0, NH
    real(dp),allocatable :: w_d(:), w_matrix(:,:)
    real(dp)  :: matrix(Ntot,Ntot), d(Ntot)
    character(len=50) :: buffer
    character(*)      :: filename

   ! Count number of lines in file 
    CALL file_lines(filename,Np)

   !// Open file for reading
    lunit = 89
    OPEN(lunit, file = filename)

   ! Calculate Ntot from header and number of lines in file
    Np = Np - NH - 3 ! Np is being repurposed as number of diameters here

    !// Skip header and read d-array
    IF (NH >= 1) THEN 
       DO i=1,NH !// Skip header
          READ(lunit,'(A)') buffer
       END DO
    ENDIF
    READ(lunit,'(A)') buffer

    ALLOCATE(w_d(Np))
    read(lunit,*) w_d(:)
   
   !// Read w-matrix
    ALLOCATE(w_matrix(Np,Np))
    read(lunit,'(A)') buffer
   
    DO x=1,Np
       read(lunit,*) w_matrix(x,:) 
    ENDDO
    CLOSE(lunit)
   
   !// INTERPOLATE TO d[Ntot] GRID
    N0 = Ntot
        CALL cubic_interpolation_2d(Np, Np, N0, N0, w_d, w_d, d, d, w_matrix, matrix, silentFLAG)
    DEALLOCATE(w_d)
    DEALLOCATE(w_matrix)
   
 END SUBROUTINE read_coef_matrix

 SUBROUTINE read_coef_array(filename, array, d,di_in,silentFLAG, Ntot, read_type,NH)
    ! Imports coagulation coefficient matrix from file "filename"  
    ! and interpolates to grid provided by "d" into "array".
    !
    ! INPUT:
    ! Ntot : Number of diameter nodes in simulation for aerosols
    ! filename : Name of file to read coefficients from
    ! d : diameter of the Ntot nodes in the simulation. "d" is the grid that
    !     the coefficents from "filename" is interpolated to.
    !            1: Read a Ntot size column from the coefficient matrix
    !            2: Read a Ntot size row from the coefficient matrix
    !
    ! NH : number of header lines to skip in "filename".
    !
    ! OUTPUT: 
    ! arrat : An Ntot size array containing interpolated coefficents
    !
    ! REQUIREMENTS FOR "filename":
    ! "filename" may have a header, which is skipped by providing number of
    ! header lines NH. Otherwise code assumes d-header, d-array, coef-header, 
    ! and coef-matrix. There must not be lines after the coef-matrix.
    !
    ! Example for Ntot = 3, and NH = 2:
    !  
    !   >HEADER HEADER HEADER HEADER
    !   >HEADER HEADER HEADER HEADER 
    !   >D-array [m]:
    !   >        7.00000E-10         8.91788E-10         1.13612E-09
    !   >Coagulation-matrix (s m^3):
    !   >        4.81600E-16         6.95071E-16         1.39363E-15
    !   >        6.95071E-16         6.34285E-16         9.15315E-16
    !   >        1.39363E-15         9.15315E-16         8.35062E-16
    !   

    USE precision_type
    USE math_module
    USE io_general_module
    IMPLICIT NONE
    integer           :: i, silentFLAG, lunit, Np, Ntot, x, N0, NH, read_type
    real(dp),allocatable :: w_d(:), w_matrix(:,:)
    real(dp)  :: array(Ntot), d(Ntot),di_in,di(1)
    character(len=50) :: buffer
    character(*)      :: filename

   ! Count number of lines in file 
    CALL file_lines(filename,Np)

   !// Open file for reading
    lunit = 89
    OPEN(lunit, file = filename)

   ! Calculate Ntot from header and number of lines in file
    Np = Np - NH - 3 ! Np is being repurposed as number of diameters here

    !// Skip header and read d-array
    IF (NH >= 1) THEN 
       DO i=1,NH !// Skip header
          READ(lunit,'(A)') buffer
       END DO
    ENDIF
    READ(lunit,'(A)') buffer

    ALLOCATE(w_d(Np))
    read(lunit,*) w_d(:)
   
   !// Read w-matrix
    ALLOCATE(w_matrix(Np,Np))
    read(lunit,'(A)') buffer
   
    DO x=1,Np
       read(lunit,*) w_matrix(x,:) 
    ENDDO
    CLOSE(lunit)
   
    di(1) = di_in

   !// INTERPOLATE TO d[Ntot] GRID
    N0 = Ntot
    IF (read_type == 1) THEN 
        CALL cubic_interpolation_2d(Np, Np, 1,  N0, w_d, w_d, di, d, w_matrix, array, silentFLAG)
    ENDIF
    IF (read_type == 2) THEN
        CALL cubic_interpolation_2d(Np, Np, N0, 1,  w_d, w_d, d, di, w_matrix, array, silentFLAG)
    ENDIF
    DEALLOCATE(w_d)
    DEALLOCATE(w_matrix)
   
 END SUBROUTINE read_coef_array

 SUBROUTINE read_coef_twocolumn(filename, coef, d, Ntot, NH)
    ! Imports condensation coefficients from file "filename"
    ! and interpolates to grid provided by "d" into "coef".
    !
    ! INPUT:
    ! Ntot : Number of diameter nodes in simulation for aerosols
    ! filename : Name of file to read coefficients from
    ! d : diameter of the Ntot nodes in the simulation. "d" is the grid that
    !     the coefficents from "filename" is interpolated to.
    ! NH : number of header lines to skip in "filename".
    !
    ! OUTPUT: 
    ! coef : An Ntot size array containing interpolated coefficents
    !
    ! REQUIREMENTS FOR "filename":
    ! "filename" may have a header, which is skipped by providing number of
    ! header lines NH. Otherwise code assumes columns two columns - diameter
    ! and coefficient. Here is an example with NH = 3 and Ntot = 3:
    !
    !   >HEADER HEADER HEADER HEADER
    !   >HEADER HEADER HEADER HEADER 
    !   >HEADER  Diameter [m]        Coagulation coefficent [m^3/s]
    !   >        7.00000E-10         4.81600E-16
    !   >        8.91788E-10         6.95071E-16
    !   >        1.13612E-09         1.39363E-15

    USE precision_type
    USE math_module
    USE io_general_module
    IMPLICIT NONE
    integer              :: NL, i, Ntot, lunit, NH
    character(len=300)   :: filename 
    character(len=300)   :: buffer
    real(dp)             :: coef(Ntot), d(Ntot)
    real(dp),allocatable :: d_grid(:), coef_grid(:)
   
    !// Find number of lines in file
    call file_lines(filename,NL)
   
    !// Open file for reading
    lunit = 88
    OPEN(lunit, file = filename, status='OLD')
       !// Allocate arrays for reading data
       ALLOCATE(d_grid(NL-NH))
       ALLOCATE(coef_grid(NL-NH))
       IF (NH > 0) THEN
          DO i=1,NH ! Skip header
             READ(lunit,'(A)') buffer
          ENDDO
       ENDIF
       DO i=1,NL-NH
          READ(lunit,*) d_grid(i),coef_grid(i)
       ENDDO
    CLOSE(lunit)

    call linear_interpolation( NL-NH, Ntot, d_grid, coef_grid, d, coef)
    
    DEALLOCATE(d_grid)
    DEALLOCATE(coef_grid)
 END SUBROUTINE read_coef_twocolumn
 
 SUBROUTINE read_coef_twocolumn_singlevalue(filename, coef_out, d_in, NH)
    ! Imports condensation coefficients from file "filename"
    ! and interpolates to grid provided by "d" into "coef".
    !
    ! INPUT:
    ! filename : Name of file to read coefficients from
    ! d_in : diameter to extract condensation coefficient for. "d" is the value 
    !        that the coefficents from "filename" is interpolated to.
    ! NH : number of header lines to skip in "filename".
    !
    ! OUTPUT: 
    ! coef_out : An scalar containing interpolated coefficent
    !
    ! REQUIREMENTS FOR "filename":
    ! "filename" may have a header, which is skipped by providing number of
    ! header lines NH. Otherwise code assumes columns two columns - diameter
    ! and coefficient. Here is an example with NH = 3 and Ntot = 3:
    !
    !   >HEADER HEADER HEADER HEADER
    !   >HEADER HEADER HEADER HEADER 
    !   >HEADER  Diameter [m]        Condensation coefficent [m^3/s]
    !   >        7.00000E-10         4.81600E-16
    !   >        8.91788E-10         6.95071E-16
    !   >        1.13612E-09         1.39363E-15

    USE precision_type
    USE math_module
    USE io_general_module
    IMPLICIT NONE
    integer              :: NL, i, Ntot, lunit, NH
    character(len=300)   :: filename 
    character(len=300)   :: buffer
    real(dp)             :: coef(1), d_in,d(1),coef_out
    real(dp),allocatable :: d_grid(:), coef_grid(:)
   
    !// Find number of lines in file
    call file_lines(filename,NL)
   
    !// Open file for reading
    lunit = 88
    OPEN(lunit, file = filename, status='OLD')
       !// Allocate arrays for reading data
       ALLOCATE(d_grid(NL-NH))
       ALLOCATE(coef_grid(NL-NH))
       IF (NH > 0) THEN
          DO i=1,NH ! Skip header
             READ(lunit,'(A)') buffer
          ENDDO
       ENDIF
       DO i=1,NL-NH
          READ(lunit,*) d_grid(i),coef_grid(i)
       ENDDO
    CLOSE(lunit)

    d(1) = d_in
    call linear_interpolation( NL-NH, 1, d_grid, coef_grid, d, coef)
    coef_out = coef(1)
    DEALLOCATE(d_grid)
    DEALLOCATE(coef_grid)
 END SUBROUTINE read_coef_twocolumn_singlevalue

 SUBROUTINE read_controlfile(charFLAG, coagFLAG, nuclFLAG, condFLAG, lossFLAG, loadFLAG, &
                             ioncFLAG, prodFLAG, P0, Jn0, Jnp, Jnm, q, alpha, n00, np0, nm0, t0, &
                             dt, ts, te, dmin, dmax, input_file, load_state_file, &
                             fname_Kklp0, fname_Kklpm, fname_Kkl00, fname_Kkl0m, fname_bkp0, &
                             fname_bkpm, fname_bkmp, fname_bkm0, fname_bk0p, fname_bk00, &
                             fname_bk0m, head_Kklp0, head_Kklpm, head_Kkl00, head_Kkl0m, & 
                             head_bkp0, head_bkpm, head_bkmp, head_bkm0, head_bk0p, &
                             head_bk00, head_bk0m, rho_p, ggamma, lambda, mc0, &
                             binFLAG, mcp, mcm, rhocp, rhocm, Ntot, np_fix_flag, nm_fix_flag, &
                             n0_fix_flag,d0_crit,dp_crit,dn_crit, d_loss, NL, d_NL,relerr,abserr,flag,t_snap)
 ! Reads parameters set in controlfile
 
   USE precision_type
   USE io_general_module
   IMPLICIT NONE
   ! Input related variables
   character(len=100) :: buffer, label
   integer :: pos
   integer, parameter :: fh = 15
   integer :: ios = 0
   integer :: line = 0
 
   ! Control file variables
   integer :: Ntot, charFLAG, coagFLAG, nuclFLAG, condFLAG, lossFLAG, loadFLAG, np_fix_flag, nm_fix_flag, n0_fix_flag 
   integer :: head_Kklp0, head_Kklpm, head_Kkl00, head_Kkl0m, head_bkp0, binFLAG, ioncFLAG, prodFLAG,flag
   integer :: head_bkpm, head_bkmp, head_bkm0, head_bk0p, head_bk00, head_bk0m
   real(dp) :: T, P0, Jn0, Jnp, Jnm, q, alpha, n00, np0, nm0, t0, dt, ts, te, mcm, mcp, rhocp, rhocm, d_loss
   real(dp) :: dmin, dmax, rho_p, ggamma, lambda, mc0, d0_crit, dp_crit, dn_crit, NL, d_NL, relerr,abserr, t_snap
 
   character(LEN=300) :: input_file, load_state_file,fname_Kklp0, fname_Kklpm, fname_Kkl00, fname_Kkl0m, fname_bkp0
   character(LEN=300) :: fname_bkpm, fname_bkmp, fname_bkm0, fname_bk0p, fname_bk00, fname_bk0m
   LOGICAL :: file_ex
   integer :: l
   INQUIRE(FILE=input_file, EXIST=file_ex)
   IF (.not. file_ex) THEN 
      write(*,*) "Control-file not found"
      RETURN
   ENDIF
 
   open(fh, file=input_file)
 
   ! ios is negative if an end of record condition is encountered or if
   ! an endfile condition was detected.  It is positive if an error was
   ! detected.  ios is zero otherwise.
 
   do while (ios == 0)
      read(fh, '(A)', iostat=ios) buffer
      if (ios == 0) then
         line = line + 1
         ! Find the first instance of whitespace.  Split label and data.
         pos = scan(buffer, '    ')
         label = buffer(1:pos)
         buffer = buffer(pos+1:)
         ! LIST OF CASES! For adding new parameters to the inputfile, the
         ! cases should be defined here!
         select case (label)
         case ('ntot')
            read(buffer, *, iostat=ios) Ntot
         case ('charge_module')
            read(buffer, *, iostat=ios) charFLAG
         case ('coagulation_module')
            read(buffer, *, iostat=ios) coagFLAG
         case ('nucleation_module')
            read(buffer, *, iostat=ios) nuclFLAG
         case ('condensation_module')
            read(buffer, *, iostat=ios) condFLAG
         case ('ion_condensation_module')
            read(buffer, *, iostat=ios) ioncFLAG
         case ('loss_module')
            read(buffer, *, iostat=ios) lossFLAG
         case ('monomer_production_module')
            read(buffer, *, iostat=ios) prodFLAG
         case ('temperature')
            read(buffer, *, iostat=ios) T 
         case ('nucleation_rate_neutral')
            read(buffer, *, iostat=ios) Jn0
         case ('nucleation_rate_positive')
            read(buffer, *, iostat=ios) Jnp
         case ('nucleation_rate_negative')
            read(buffer, *, iostat=ios) Jnm
         case ('neutral_production_rate')
            read(buffer, *, iostat=ios) P0
         case ('ion_pair_production_rate')
            read(buffer, *, iostat=ios) q
         case ('recombination_coefficient')
            read(buffer, *, iostat=ios) alpha
         case ('density_of_aerosols')
            read(buffer, *, iostat=ios) rho_p
         case ('loss_term_gamma')
            read(buffer, *, iostat=ios) ggamma
         case ('loss_term_lambda')
            read(buffer, *, iostat=ios) lambda
         case ('loss_term_d_loss')
            read(buffer, *, iostat=ios) d_loss
         case ('neutral_mass')
            read(buffer, *, iostat=ios) mc0 
         case ('initial_ion(+)_concentr')
            read(buffer, *, iostat=ios) np0 
         case ('initial_ion(-)_concentr')
            read(buffer, *, iostat=ios) nm0
         case ('initial_neutral_concentr')
            read(buffer, *, iostat=ios) n00
         case('fix_initials')
            read(buffer, *, iostat=ios) np_fix_flag,nm_fix_flag,n0_fix_flag
         case ('start_time')
            read(buffer, *, iostat=ios) t0
         case ('end_time')
            read(buffer, *, iostat=ios) te
!        case ('time_step_length')
!           read(buffer, *, iostat=ios) dt
         case ('make_snapshot_every')
            read(buffer, *, iostat=ios) t_snap
         case ('integrator_step_length')
            read(buffer, *, iostat=ios) ts
         case ('minimum_aerosol_diameter')
            read(buffer, *, iostat=ios) dmin
         case ('maximum_aerosol_diameter')
            read(buffer, *, iostat=ios) dmax
         case ('positive_ion_mass')
            read(buffer, *, iostat=ios) mcp
         case ('negative_ion_mass')
            read(buffer, *, iostat=ios) mcm
         case ('positive_ion_density')
            read(buffer, *, iostat=ios) rhocp
         case ('negative_ion_density')
            read(buffer, *, iostat=ios) rhocm
         case ('neutral_critical_diameter')
            read(buffer, *, iostat=ios) d0_crit
         case ('positive_critical_diameter')
            read(buffer, *, iostat=ios) dp_crit
         case ('negative_critical_diameter')
            read(buffer, *, iostat=ios) dn_crit
         case ('lin_log_bin_flag') 
            read(buffer, *, iostat=ios) binFLAG
         case ('diameter_large_aerosol_loss') 
            read(buffer, *, iostat=ios) d_NL
         case ('concentr_large_aerosol_loss') 
            read(buffer, *, iostat=ios) NL
         case ('relerr') 
            read(buffer, *, iostat=ios) relerr
         case ('abserr') 
            read(buffer, *, iostat=ios) abserr
         case ('flag') 
            read(buffer, *, iostat=ios) flag
         case ('n_header_Kklp0')
            read(buffer, *, iostat=ios) head_Kklp0
         case ('filename_Kklp0')
            read(buffer, '(A)', iostat=ios) fname_Kklp0        ! Read buffer string into variable string
            fname_Kklp0 = adjustl(fname_Kklp0)                 ! Remove space before path
            l = findnl(fname_Kklp0)                            ! Locate placement of "newline" character in string
            if (l .le. len(fname_Kklp0)) fname_Kklp0(l:) = " " ! Substitute from "newline" to end of string with " " 
         case ('n_header_Kklpm')
            read(buffer, *, iostat=ios) head_Kklpm
         case ('filename_Kklpm')
            read(buffer, '(A)', iostat=ios) fname_Kklpm
            fname_Kklpm = adjustl(fname_Kklpm)
            l = findnl(fname_Kklpm)
            if (l .le. len(fname_Kklpm)) fname_Kklpm(l:) = " " 
         case ('n_header_Kkl00')
            read(buffer, *, iostat=ios) head_Kkl00
         case ('filename_Kkl00')
            read(buffer, '(A)', iostat=ios) fname_Kkl00
            fname_Kkl00 = adjustl(fname_Kkl00)
            l = findnl(fname_Kkl00)             
            if (l .le. len(fname_Kkl00)) fname_Kkl00(l:) = " "      
         case ('n_header_Kkl0m')
            read(buffer, *, iostat=ios) head_Kkl0m
         case ('filename_Kkl0m')
            read(buffer, '(A)', iostat=ios) fname_Kkl0m
            fname_Kkl0m = adjustl(fname_Kkl0m)
            l = findnl(fname_Kkl0m)
            if (l .le. len(fname_Kkl0m)) fname_Kkl0m(l:) = " "      
         case ('n_header_bkp0')
            read(buffer, *, iostat=ios) head_bkp0
         case ('filename_bkp0')
            read(buffer, '(A)', iostat=ios) fname_bkp0
            fname_bkp0 = adjustl(fname_bkp0)
            l = findnl(fname_bkp0)
            if (l .le. len(fname_bkp0)) fname_bkp0(l:) = " "      
         case ('n_header_bkpm')
            read(buffer, *, iostat=ios) head_bkpm
         case ('filename_bkpm')
            read(buffer, '(A)', iostat=ios) fname_bkpm
            fname_bkpm = adjustl(fname_bkpm)
            l = findnl(fname_bkpm)
            if (l .le. len(fname_bkpm)) fname_bkpm(l:) = " "      
         case ('n_header_bkmp')
            read(buffer, *, iostat=ios) head_bkmp
         case ('filename_bkmp')
            read(buffer, '(A)', iostat=ios) fname_bkmp
            fname_bkmp = adjustl(fname_bkmp)
            l = findnl(fname_bkmp)
            if (l .le. len(fname_bkmp)) fname_bkmp(l:) = " "      
         case ('n_header_bkm0')
            read(buffer, *, iostat=ios) head_bkm0
         case ('filename_bkm0')
            read(buffer, '(A)', iostat=ios) fname_bkm0
            fname_bkm0 = adjustl(fname_bkm0)
            l = findnl(fname_bkm0)             
            if (l .le. len(fname_bkm0)) fname_bkm0(l:) = " "      
         case ('n_header_bk0p')
            read(buffer, *, iostat=ios) head_bk0p
         case ('filename_bk0p')
            read(buffer, '(A)', iostat=ios) fname_bk0p 
            fname_bk0p = adjustl(fname_bk0p)
            l = findnl(fname_bk0p)
            if (l .le. len(fname_bk0p)) fname_bk0p(l:) = " "      
         case ('n_header_bk00')
            read(buffer, *, iostat=ios) head_bk00
         case ('filename_bk00')
            read(buffer, '(A)', iostat=ios) fname_bk00
            fname_bk00 = adjustl(fname_bk00)
            l = findnl(fname_bk00)             
            if (l .le. len(fname_bk00)) fname_bk00(l:) = " "      
         case ('n_header_bk0m')
            read(buffer, *, iostat=ios) head_bk0m
         case ('filename_bk0m')
            read(buffer, '(A)', iostat=ios) fname_bk0m
            fname_bk0m = adjustl(fname_bk0m)
            l = findnl(fname_bk0m)             
            if (l .le. len(fname_bk0m)) fname_bk0m(l:) = " "      
         case ('load_file_flag')
            read(buffer, *, iostat=ios) loadFLAG 
         case ('filename_load')
            read(buffer, '(A)', iostat=ios) load_state_file
            load_state_file = adjustl(load_state_file)
            l = findnl(load_state_file)
            if (l .le. len(load_state_file)) load_state_file(l:) = " "      
         case default
         end select
      end if
   end do
 
  close(fh)
 
 !/////////////////////////
 !// DO UNIT CONVERSIONS //
 !/////////////////////////
 ! The code uses SI units, so controlfile
 ! units are here converted where needed.
  P0    = P0    * 1.0e6
  Jn0   = Jn0   * 1.0e6
  Jnp   = Jnp   * 1.0e6
  Jnm   = Jnm   * 1.0e6
  q     = q     * 1.0e6
  alpha = alpha * 1.0e-6
  n00   = n00   * 1.0e6
  np0   = np0   * 1.0e6
  nm0   = nm0   * 1.0e6
  mc0   = mc0   * 1.66053904e-27
  mcp   = mcp   * 1.66053904e-27
  mcm   = mcm   * 1.66053904e-27
  d_NL  = d_NL  * 1.0e-9
  NL    = NL    * 1.0e6

 
 END SUBROUTINE read_controlfile
 
 SUBROUTINE print_controlfile_to_output(input_file)
   ! Prints contents of inputfile to terminal or outputfile
    IMPLICIT NONE
    CHARACTER(*) :: input_file
    CHARACTER(LEN=150) :: buffer
    integer :: ios = 0
    integer :: lunit = 99
    OPEN (lunit, file = input_file)                    ! Open file for reading 
    DO                                                 ! Loop until end of file in READ
       READ (lunit,'(A)', iostat=ios, END=10) buffer   ! Read line into buffer
       write(*,'(A)') buffer                           ! Print line 
       buffer(:) = " "                                 ! Clear buffer
    END DO
    10 CLOSE (lunit)
   
 END SUBROUTINE print_controlfile_to_output
 
 SUBROUTINE print_snapshot( Nk0, Nkp, Nkm, n0, np, nm, t, Ntot)
   ! Prints a snapshot to terminal or outputfile
    USE precision_type
    IMPLICIT NONE
    integer :: i, Ntot
    real(dp) :: Nk0(Ntot), Nkp(Ntot), Nkm(Ntot), n0, np, nm, t
    write(*,*) "--------------------------------------------" 
    write(*,*) " "
    write(*,*) "Time [s]" 
    write(*,"(F22.9)") t 
    write(*,*) " " 
    write(*,*) "n0" 
    write(*,"(ES16.6)") n0 
    write(*,*) " " 
    write(*,*) "np              nm"
    write(*,"(2ES16.6)") np, nm
    write(*,*) " "
    write(*,*) "Nk0             Nkp             Nkm"
    DO i=1,Ntot 
       write(*,"(3ES16.6E3)") Nk0(i), Nkp(i), Nkm(i)
    ENDDO
    FLUSH(6)
 END SUBROUTINE print_snapshot

 SUBROUTINE load_state(filename, vmin, vmax, Nk0, Nkp, Nkm, n0, np, nm, Ntot)
 ! Loads the latest snapshot of an outputfile from a previous run, and continues
 ! simulation from that state.
  USE precision_type
  USE io_general_module
  IMPLICIT NONE
  integer                    :: Ntot,lunit,NP_file,Nsnaps,i,pos,number_of_lines
  character(len=300)         :: filename,buffer
  character(len=32)          :: smallbuf,label
  real(dp)                   :: vmin, vmax, n0, np, nm, pi,dmin_file,dmax_file
  real(dp),dimension(Ntot)   :: Nk0, Nkp, Nkm
 
  pi = 3.141592653589793
 
  !// Open file for reading and check for size compatibility //
  lunit = 29
  OPEN(lunit, file = filename)
     DO i=1,5 ! Skip first 5 lines
        read(lunit,"(A)") buffer
     ENDDO
 
     ! Read Ntot from old inputfile    
     read(lunit,"(A)") buffer
     pos = scan(buffer, '    ')
     label = buffer(1:pos)
     buffer = buffer(pos+1:)
     read(buffer, *) Np_file
     DO i=1,47 ! Skip lines
        read(lunit,"(A)") buffer
     ENDDO
 
     ! Read d_min from old inputfile    
     read(lunit,"(A)") buffer
     pos = scan(buffer, '    ')
     label = buffer(1:pos)
     buffer = buffer(pos+1:)
     read(buffer, *) dmin_file
     ! Read d_max from old inputfile    
     read(lunit,"(A)") buffer
     pos = scan(buffer, '    ')
     label = buffer(1:pos)
     buffer = buffer(pos+1:)
     read(buffer, *) dmax_file
 
     IF (Np_file .ne. Ntot) write(*,*) "load_state: Error: Size array not compatible"
     IF (abs((4.0d0*pi*(1.0d-9*dmin_file/2.0d0)**3)/3.0d0 - vmin) > 0.001d0*vmin) &
        write(*,*) 'load_state: Error: vmin not compatible'
     IF (abs((4.0d0*pi*(1.0d-9*dmax_file/2.0d0)**3)/3.0d0 - vmax) > 0.001d0*vmax) &
        write(*,*) 'load_state: Error: vmax not compatible'
 
  CLOSE(lunit)

  call file_lines(filename,number_of_lines)
 
  lunit = 39
 !/ Read final snapshot into arrays
  OPEN(lunit, file = filename)
     Nsnaps = (number_of_lines-(89+Ntot+2))/(Ntot+12)
     DO i=1, 93+(2+Ntot)+(Nsnaps-1)*(Ntot+12)
        read(lunit,"(A)") buffer !// Skip to final snapshot
     ENDDO
     DO i=1,7
        read(lunit,"(A)") buffer !// Skip header of final snapshot
     ENDDO
     read(lunit,*) n0            !// Read n0
     read(lunit,"(A)") buffer    !// Skip line
     read(lunit,"(A)") buffer    !// Skip line
     read(lunit,*) np,nm         !// Read np
     read(lunit,"(A)") buffer    !// Skip line
     read(lunit,"(A)") buffer    !// Skip line
 
     DO i=1,Ntot
        read(lunit,*) Nk0(i), Nkp(i),Nkm(i)
     ENDDO
  CLOSE(lunit)
 END SUBROUTINE load_state
 
END MODULE io_ioncage_module
