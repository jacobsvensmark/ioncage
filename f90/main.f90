
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!! FORTRAN version v0.9 of the ION-CAGE code.  !!
!! // Jacob Svensmark                          !!
!! In case of using the following code, please !!
!! add a reference to                          !!
!!    arXiv:1909.12784                         !!
!! or its peer reviewed counterpart as it is   !!
!! accepted for publication                    !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

PROGRAM main

 USE precision_type
 USE io_general_module
 USE io_ioncage_module
 USE math_module
 USE calc_module
 USE dNdt_module
 USE shared_data
 USE rk45_mod

 IMPLICIT NONE

!//////////////////////////
!/ Variable Declarations //
!//////////////////////////
 integer :: i, j, neqn, flag
 integer :: head_Kklp0,head_Kklpm,head_Kkl00,head_Kkl0m, head_bkp0, head_bkpm 
 integer :: head_bkmp, head_bkm0, head_bk0p, head_bk00, head_bk0m, binFLAG, snapFLAG
 real(dp) :: d0_crit, dp_crit, dm_crit
 real(dp) :: v0_crit, vp_crit, vm_crit
 real(dp) :: vmin, vmax, n0, np, nm, dn0dt, dnpdt, dnmdt, tt, ts_old, t_snap
 real(dp) :: mcp, mcm, rhocp, rhocm, d_NL, abserr, relerr
 real(dp) :: n00, np0, nm0, t0, dt, ts, te, pi
 real(dp) :: dmin, dmax, Jn0, rho_p, mc0
 real(dp),allocatable :: Nk0(:),Nkp(:),Nkm(:),dNk0dt(:),dNkpdt(:),dNkmdt(:), y(:), yp(:)
 character(LEN=300) :: load_state_file, fname_Kklp0, fname_Kklpm, fname_Kkl00, fname_Kkl0m, fname_bkp0 
 character(LEN=300) :: fname_bkpm, fname_bkmp, fname_bkm0, fname_bk0p, fname_bk00, fname_bk0m
 real(dp) :: fac, Dp0, sigma, d1,d2,v1,v2, NN0
 external calculate_derivatives
!///////////////////////////
!/     Read input-file     //
!////////////////////////////
 character(LEN=300) :: input_file = 'controlfile.txt' !/ Name of input-file
 call read_controlfile(charFLAG, coagFLAG, nuclFLAG, condFLAG, lossFLAG, loadFLAG, ioncFLAG, & 
                       prodFLAG, P0, Jn0, Jnp, Jnm, q, alpha, n00,  np0, nm0, t0, dt, ts, te,  &
                       dmin, dmax, input_file, load_state_file, & 
                       fname_Kklp0, fname_Kklpm, fname_Kkl00, fname_Kkl0m, fname_bkp0, & 
                       fname_bkpm, fname_bkmp, fname_bkm0, fname_bk0p, fname_bk00, fname_bk0m,&
                       head_Kklp0, head_Kklpm, head_Kkl00, head_Kkl0m, head_bkp0,&
                       head_bkpm, head_bkmp, head_bkm0, head_bk0p, head_bk00, &
                       head_bk0m, rho_p, ggamma, lambda, mc0 , binFLAG, &
                       mcp, mcm, rhocp, rhocm, Ntot, np_fix_flag, nm_fix_flag, n0_fix_flag, &
                       d0_crit, dp_crit, dm_crit, d_loss, NL, d_NL, relerr, abserr, flag, t_snap)

!/////////////////////////////
!/  Alocate Arrays to Ntot  //
!/////////////////////////////
 allocate(Vij(Ntot,Ntot))
 allocate(Nkp(Ntot),Nk0(Ntot),Nkm(Ntot),v(Ntot),d(Ntot),bkp0(Ntot),bkpm(Ntot))
 allocate(dNk0dt(Ntot),dNkpdt(Ntot),dNkmdt(Ntot))
 allocate(bk0p(Ntot),bk00(Ntot),bk0m(Ntot),bkmp(Ntot),bkm0(Ntot))
 allocate(kL0(Ntot),kLp(Ntot),kLm(Ntot))
 allocate(Kklp0(Ntot,Ntot),Kklpm(Ntot,Ntot),Kkl00(Ntot,Ntot),Kkl0m(Ntot,Ntot),S(Ntot,Ntot))
 allocate(y(3*Ntot+3),yp(3*Ntot+3))

!/****************************************************************************\\
!/****************************************************************************\\
!/****************************************************************************\\

!////////////////////////////
!/ CREATE LOG-SPACED NODES //
!////////////////////////////
 pi = 3.141592653589793d0
 vmin = (dmin/2.0d0*1.0d-9)**3 * (4.0d0/3.0d0*pi)
 vmax = (dmax/2.0d0*1.0d-9)**3 * (4.0d0/3.0d0*pi)
 binFLAG = 0
 IF (binFLAG == 0) THEN  ! Nodes logarithmically spaced in v
    DO j = 1, Ntot
      v(j) = dble(10.)**(dlog10(vmin) + dble(j-1)*(dlog10(vmax)-dlog10(vmin))/dble(Ntot-1))
      d(j) = 2.0d0*(v(j)*3.0d0/(4.0d0*pi))**(1.0d0/3.0d0)
    ENDDO
 END IF
 IF (binFLAG == 1) THEN ! Nodes linearly spaced in d
    DO j = 1, Ntot
       d(j) = 1.0d-9*(dmin+dble(j-1)*(dmax-dmin)/dble(Ntot))
       v(j) = (4.d0/3.d0)*pi*(d(j)/2.D0)**3.0D0
    ENDDO
 END IF

 vcm = mcm/rhocm
 vcp = mcp/rhocp
 v0 = mc0/rho_p
 d_n0 = 2.0d0*(v0*3.0d0/(4.0d0*pi))**(1.0d0/3.0d0)
 d_np = 2.0d0*(vcp*3.0d0/(4.0d0*pi))**(1.0d0/3.0d0)
 d_nm = 2.0d0*(vcm*3.0d0/(4.0d0*pi))**(1.0d0/3.0d0)

! Check size nucleating particles (SECTION NOT DONE YET!!!)

! Neutral nucleation
 v0_crit = (d0_crit/2.*1.0d-9)**3 * (4.0d0/3.0d0*pi)
 IF (v0_crit <= v(1)) THEN
    Jn0_scaled = v0_crit/v(1) * Jn0
    i_Jn0 = 1
 ELSE
    DO i=1,NTOT-1
       IF ((v0_crit > v(i)) .AND. (v0_crit <= v(i+1))) THEN
          Jn0_scaled = v0_crit/v(i+1) * Jn0
          i_Jn0 = i+1
       END IF 
    ENDDO 
 END IF

! Positive nucleation
 vp_crit = (dp_crit/2.*1.0d-9)**3 * (4.0d0/3.0d0*pi)
 IF (vp_crit <= v(1)) THEN
    Jnp_scaled = vp_crit/v(1) * Jnp
    i_Jnp = 1
 ELSE
    DO i=1,NTOT-1
       IF ((vp_crit > v(i)) .AND. (vp_crit <= v(i+1))) THEN
          Jnp_scaled = vp_crit/v(i+1) * Jnp
          i_Jnp = i+1
       END IF 
    ENDDO 
 END IF

! Negative nucleation
 vm_crit = (dm_crit/2.*1.0d-9)**3 * (4.0d0/3.0d0*pi)
 IF (vm_crit <= v(1)) THEN
    Jnm_scaled = vm_crit/v(1) * Jnm
    i_Jnm = 1
 ELSE
    DO i=1,NTOT-1
       IF ((vm_crit > v(i)) .AND. (vm_crit <= v(i+1))) THEN
          Jnm_scaled = vm_crit/v(i+1) * Jnm
          i_Jnm = i+1
       END IF 
    ENDDO 
 END IF

! Check if volume of ions or neutral monomer is larger than smallest bin spacing

!IF ((v(2) - v(1) < vcm) .OR. (v(2) - v(1) < vcp) .OR. (v(2) - v(1) < v0)) THEN
!   write(*,*) "*******************************************************************"
!   write(*,*) "   Error, smallest bin spacing has volume less than monomers,"
!   write(*,*) "   either ion volume or sulfuric acid. Change vmin, vmax or Ntot"
!   write(*,*) "   such that v(2)-v(1) is larger than monomer volumes."
!   write(*,*) "*******************************************************************"
!   GOTO 90
!END IF 

IF (ioncFLAG == 1) THEN
   IF (charFLAG == 0) THEN
       write(*,*) '***************************************'
       write(*,*) ' WARNING: ioncFLAG = 1 requires that'
       write(*,*) '          charFLAG = 1. Comment out '
       write(*,*) '          this requirement in main.f90' 
       write(*,*) '          to ignore.'
       write(*,*) '***************************************'
       GOTO 90
   END IF
END IF

!/////////////////////////
!/ Initial distribution // 
!/////////////////////////

! Clear arrays
 DO i = 1,Ntot
    Nk0(i) = 0
    Nkp(i) = 0
    Nkm(i) = 0
 ENDDO

 DO i=1,Ntot
    dNk0dt(i) = 0
    dNkpdt(i) = 0
    dNkmdt(i) = 0
 ENDDO
 dn0dt = 0
 dnpdt = 0
 dnmdt = 0

 neqn = 3*Ntot+3
 DO i=1,neqn
    yp(i) = 0
 ENDDO 
 
!// Read datafile, continue from end of run and overwrite initial distributions set above if loadFLAG=1 
 if (loadFLAG == 1) call load_state(load_state_file, vmin, vmax, Nk0, Nkp, Nkm, n0, np, nm, Ntot)

! Set initial values for n0, np and nm as specified in controlfile
 n0 = n00
 np = np0
 nm = nm0

!///////////////////////////////////////////////////
!/ Make coagulation volume and splitting matrices // 
!///////////////////////////////////////////////////
 DO i = 1, Ntot 
    DO j = 1, Ntot
       Vij(i,j) = floor(((dlog10(v(i)+v(j))-dlog10(vmin))*dble(Ntot-1))/(dlog10(vmax)-dlog10(vmin))) + 1
       S(i,j) = (v(Vij(i,j)+1)-(v(i)+v(j))) / (v(Vij(i,j)+1)-v(Vij(i,j)))
       ! Handle splitting when aerosols coagulate out of diameter node range
       IF (Vij(i,j) == Ntot ) THEN 
          S(i,j) = ( v(Ntot)  * (v(Ntot) / v(Ntot-1))-(v(i)+v(j))) / ( v(Ntot)  *(v(Ntot)  / v(Ntot-1))-v(Vij(i,j)) )
       ELSE IF (Vij(i,j) > Ntot) THEN
          S(i,j) = 0.0d0
          Vij(i,j) = Ntot
       ENDIF
    ENDDO
 ENDDO

!////////////////////////////////////////////////////////
!/ Load coagulation and condensation interaction terms //
!////////////////////////////////////////////////////////
! Load condensation coefficients 
 call read_coef_twocolumn( fname_bk00, bk00, d, Ntot, head_bk00 )
 call read_coef_twocolumn( fname_bkp0, bkp0, d, Ntot, head_bkp0 )
 call read_coef_twocolumn( fname_bkm0, bkm0, d, Ntot, head_bkm0 )
 call read_coef_twocolumn( fname_bk0p, bk0p, d, Ntot, head_bk0p )
 call read_coef_twocolumn( fname_bk0m, bk0m, d, Ntot, head_bk0m )
 call read_coef_twocolumn( fname_bkpm, bkpm, d, Ntot, head_bkpm )
 call read_coef_twocolumn( fname_bkmp, bkmp, d, Ntot, head_bkmp )
! Load coagulation coefficients
 call read_coef_matrix( fname_Kklp0, Kklp0, d, 0, Ntot, head_Kklp0 )
 call read_coef_matrix( fname_Kklpm, Kklpm, d, 0, Ntot, head_Kklpm )
 call read_coef_matrix( fname_Kkl00, Kkl00, d, 0, Ntot, head_Kkl00 )
 call read_coef_matrix( fname_Kkl0m, Kkl0m, d, 0, Ntot, head_Kkl0m )
! Calculate interaction terms for the large mode particles:
 call read_coef_twocolumn_singlevalue( fname_bk00, bL0, d_NL, head_bk00 )
 call read_coef_twocolumn_singlevalue( fname_bkp0, bLp, d_NL, head_bkp0 )
 call read_coef_twocolumn_singlevalue( fname_bkm0, bLm, d_NL, head_bkm0 )
 call read_coef_array( fname_Kkl00, kL0, d, d_NL, 0, Ntot, 1, head_Kkl00 )
 call read_coef_array( fname_Kklp0, kLp, d, d_NL, 0, Ntot, 1, head_Kklp0 )
 call read_coef_array( fname_Kklp0, kLm, d, d_NL, 0, Ntot, 1, head_Kklp0 ) ! Note assume symmetric for pos/neg

!//////////////////////////////////////////////////////////////////////////////////////////////////////////
!/ TIME INTEGRATION LOOP //                             START                                            //
!//////////////////////////                                                                              //

 write(*,*) "//*************************** INPUT FILE *************************//"
 call print_controlfile_to_output(input_file)
 write(*,*) "//*************************** INPUT FILE *************************//"
 write(*,*) "//********************* Diameter/Volume nodes ********************//"
 write(*,*) "Diameter nodes [m], Volume nodes [m^3]"
 DO i=1,Ntot
    write(*,"(2ES16.6)") d(i), v(i) 
 ENDDO
 write(*,*) "//********************* Diameter/Volume nodes ********************//"
 write(*,*) " "
 
 DO i=1,Ntot
    y(i)        = Nk0(i)
    y(i+Ntot)   = Nkp(i)
    y(i+2*Ntot) = Nkm(i)
 ENDDO
 y(3*Ntot+1) = n0
 y(3*Ntot+2) = np
 y(3*Ntot+3) = nm

 tt = t0 
 ts_old = ts 
 snapFLAG = 0
 call print_snapshot( Nk0, Nkp, Nkm, n0, np, nm, tt, Ntot)
 DO WHILE (tt<te)
    call r8_rkf45 ( calculate_derivatives, neqn, y, yp, tt, tt+ts, relerr, abserr, flag )
    IF (flag == 2) THEN ! Normal operation
        GOTO 30
    ELSE IF (flag == 3) THEN ! RELERROR ADJUSTED UP
        flag = 2
    ELSE IF (flag == 4) THEN ! ALLOW FOR MORE INTEGRATION STEPS
        GOTO 30
    ELSE IF (flag == 5) THEN ! SOLUTION VANISHED
        write(*,*) "rk45 flag=5: Solution vanished, making error check impossible. Consider single step mode." 
        GOTO 90
    ELSE IF (flag == 6) THEN
        write(*,*) "rk45 flag=6: Integration was not completed because the requested accuracy could not be achieved"
        GOTO 90
    ELSE IF (flag == 7) THEN
        flag = 2 ! Replace "flag = 2" with "GOTO 90" to not insist on not using single step mode.
    ELSE IF (flag == 8) THEN
        write(*,*) "rk45 flag=8: Invalid input."
        GOTO 90
    END IF

    ! Check if snapshot is expected from next rk45 call, and adjust step length accordingly
30  IF (mod(tt+ts,t_snap) < ts) THEN
        IF (snapFLAG == 0) THEN
            ts_old = ts
            snapFLAG = 1 ! snapFLAG = 1: Snapshot expected after next rkf45 call
        ELSE 
            snapFLAG = 2 ! snapFLAG = 2; Snapshot expected, but ts_old is not updated
        END IF
        ts = t_snap - mod(tt,t_snap)
    ELSE IF (snapFLAG > 0) THEN !/  PRINT SNAPSHOT
        DO i=1,Ntot
           Nk0(i) = y(i)
           Nkp(i) = y(i+Ntot)
           Nkm(i) = y(i+2*Ntot)
        ENDDO
        n0 = y(3*Ntot+1)
        np = y(3*Ntot+2)
        nm = y(3*Ntot+3)
        call print_snapshot( Nk0, Nkp, Nkm, n0, np, nm, tt, Ntot)
        snapFLAG = 0
        ts = ts_old
    END IF
    
 ENDDO
!//////////////////////////                                                                              //
!/ TIME INTEGRATION LOOP //                              END                                             //
!//////////////////////////////////////////////////////////////////////////////////////////////////////////

 90 CONTINUE

! Deallocate arrays
 deallocate(Vij)
 deallocate(Nkp, Nk0, Nkm, v, d, bkp0, bkpm, bk0p, bk00, bk0m, bkmp, bkm0)
 deallocate(kL0,kLp,kLm)
 deallocate(S, Kklp0,  Kklpm, Kkl00, Kkl0m)
 deallocate(y,yp)

END


!////////////////////////////////////
!/ SUBROUTINE FOR TOTAL DERIVATIVE //
!////////////////////////////////////

SUBROUTINE calculate_derivatives(t,y,yp)
   USE precision_type
   USE dNdt_module
   USE shared_data, ONLY: Ntot,charFLAG,coagFLAG,nuclFLAG,condFLAG,lossFLAG,loadFLAG,ioncFLAG,prodFLAG,&
                          np_fix_flag, nm_fix_flag, n0_fix_flag 

   IMPLICIT NONE

   integer :: i
   real(dp) :: n0, np, nm, dn0dt, dnpdt, dnmdt, t
   real(dp), dimension(Ntot) :: Nk0, Nkp, Nkm, dNk0dt, dNkpdt, dNkmdt
   real(dp), dimension(3*Ntot+3) :: y, yp
   DO i=1,Ntot
      Nk0(i) = y(i)
      Nkp(i) = y(i+Ntot)
      Nkm(i) = y(i+2*Ntot)
      dNk0dt(i) = 0
      dNkpdt(i) = 0
      dNkmdt(i) = 0
   ENDDO
   n0 = y(3*Ntot+1)
   np = y(3*Ntot+2)
   nm = y(3*Ntot+3)
   dn0dt = 0
   dnpdt = 0
   dnmdt = 0
   !/ Calculate desired derivatives
   IF (nuclFLAG==1) call dNdt_nucl( dNk0dt, dNkpdt, dNkmdt, dnmdt, dnpdt)
   IF (coagFLAG==1) call dNdt_coag( Nk0, Nkp, Nkm, dNk0dt, dNkpdt, dNkmdt)
   IF (condFLAG==1) call dNdt_cond( Nk0, Nkp, Nkm, dNk0dt, dNkpdt, dNkmdt, n0, dn0dt)
   IF (ioncFLAG==1) call dNdt_ionc( Nk0, Nkp, Nkm, dNk0dt, dNkpdt, dNkmdt, nm, np)
   IF (charFLAG==1) call dNdt_char( Nk0, Nkp, Nkm, dNk0dt, dNkpdt, dNkmdt, nm, np, dnmdt, dnpdt)
   IF (lossFLAG==1) call dNdt_loss( Nk0, Nkp, Nkm, dNk0dt, dNkpdt, dNkmdt, n0, nm, np, dn0dt, dnpdt, dnmdt)
   IF (prodFLAG==1) call dNdt_prod( dn0dt, dnpdt, dnmdt)
   !/ Fix monomers if fixflags are set
   IF (n0_fix_flag == 1) dn0dt = 0
   IF (np_fix_flag == 1) dnpdt = 0
   IF (nm_fix_flag == 1) dnmdt = 0
   DO i=1,Ntot
      yp(i)        = dNk0dt(i)
      yp(i+Ntot)   = dNkpdt(i)
      yp(i+Ntot*2) = dNkmdt(i)
   ENDDO
   yp(3*Ntot+1) = dn0dt
   yp(3*Ntot+2) = dnpdt
   yp(3*Ntot+3) = dnmdt
END SUBROUTINE calculate_derivatives
