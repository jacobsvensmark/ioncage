!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! FORTRAN version v0.9 of the ION-CAGE code.  !!
!! // Jacob Svensmark                          !!
!! In case of using the following code, please !!
!! add a reference to the work                 !!
!! XXXX Article_reference XXXX                 !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

PROGRAM main

 USE precision_type
 USE io_general_module
 USE io_ioncage_module
 USE math_module
 USE calc_module
 USE dNdt_module

 IMPLICIT NONE

!//////////////////////////
!/ Variable Declarations //
!//////////////////////////
 integer :: i,j,rk,Ntot, np_fix_flag, nm_fix_flag, n0_fix_flag 
 integer :: charFLAG,coagFLAG,nuclFLAG,condFLAG,lossFLAG,loadFLAG,ioncFLAG,prodFLAG
 integer :: head_Kklp0,head_Kklpm,head_Kkl00,head_Kkl0m, head_bkp0, head_bkpm 
 integer :: head_bkmp, head_bkm0, head_bk0p, head_bk00, head_bk0m, binFLAG
 integer :: i_Jn0, i_Jnp, i_Jnm
 integer,allocatable :: Vij(:,:) 
 real(dp) :: d0_crit, dp_crit, dm_crit, Jn0_scaled, Jnp_scaled, Jnm_scaled
 real(dp) :: v0_crit, vp_crit, vm_crit, d_loss
 real(dp) :: rk_coef(4) = (/1.0d0,2.0d0,2.0d0,1.0d0/)
 real(dp) :: vmin, vmax,n0, np, nm, n0_rk,  np_rk,  nm_rk, tt, v0
 real(dp) :: mcp, mcm, rhocp, rhocm, vcm,vcp, NL, d_NL, bL0, bLp, bLm
 real(dp) :: q, alpha, n00, np0, nm0, t0, dt, ts, te, ggamma, lambda, pi
 real(dp) :: dmin, dmax, d_np, d_nm, d_n0, P0, Jn0, Jnp, Jnm, rho_p, mc0, Nx, Vx
 real(dp),dimension(4) :: fn0, fnp, fnm
 real(dp),allocatable :: fNk0(:,:),fNkp(:,:),fNkm(:,:) 
 real(dp),allocatable :: v(:),d(:),bkp0(:),bkpm(:),bk0p(:),bk00(:),bk0m(:),bkmp(:),bkm0(:)
 real(dp),allocatable :: Nk0(:),Nkp(:),Nkm(:),Nk0_rk(:),Nkp_rk(:),Nkm_rk(:)
 real(dp),allocatable :: S(:,:),Kklp0(:,:),Kklpm(:,:),Kkl00(:,:),Kkl0m(:,:), kL0(:), kLp(:), kLm(:)
 character(LEN=300) :: load_state_file, fname_Kklp0, fname_Kklpm, fname_Kkl00, fname_Kkl0m, fname_bkp0 
 character(LEN=300) :: fname_bkpm, fname_bkmp, fname_bkm0, fname_bk0p, fname_bk00, fname_bk0m
 real(dp) :: fac, Dp0, sigma, d1,d2,v1,v2, NN0
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
                       d0_crit, dp_crit, dm_crit, d_loss, NL, d_NL)

!/////////////////////////////
!/  Alocate Arrays to Ntot  //
!/////////////////////////////
 allocate(Vij(Ntot,Ntot))
 allocate(fNk0(4,Ntot),fNkp(4,Ntot),fNkm(4,Ntot))
 allocate(Nkp(Ntot),Nk0(Ntot),Nkm(Ntot),v(Ntot),d(Ntot),bkp0(Ntot),bkpm(Ntot))
 allocate(bk0p(Ntot),bk00(Ntot),bk0m(Ntot),bkmp(Ntot),bkm0(Ntot))
 allocate(Nk0_rk(Ntot),Nkp_rk(Ntot),Nkm_rk(Ntot),kL0(Ntot),kLp(Ntot),kLm(Ntot))
 allocate(Kklp0(Ntot,Ntot),Kklpm(Ntot,Ntot),Kkl00(Ntot,Ntot),Kkl0m(Ntot,Ntot),S(Ntot,Ntot))

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

 IF ((v(2) - v(1) < vcm) .OR. (v(2) - v(1) < vcp) .OR. (v(2) - v(1) < v0)) THEN
    write(*,*) "*******************************************************************"
    write(*,*) "   Error, smallest bin spacing has volume less than monomers,"
    write(*,*) "   either ion volume or sulfuric acid. Change vmin, vmax or Ntot"
    write(*,*) "   such that v(2)-v(1) is larger than monomer volumes."
    write(*,*) "*******************************************************************"
    GOTO 90
 END IF 

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
 DO j = 1,Ntot
    Nk0(j) = 0.0d0 
    Nkp(j) = 0.0d0
    Nkm(j) = 0.0d0
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

 tt = t0 
 call print_snapshot( Nk0, Nkp, Nkm, n0, np, nm, tt, Ntot)
 DO WHILE (tt<te)
    tt = tt + dt
    DO i=1,Ntot !/ Set RK-N's to current value 
       Nk0_rk(i) = Nk0(i)
       Nkp_rk(i) = Nkp(i)
       Nkm_rk(i) = Nkm(i)
    ENDDO
    np_rk = np
    nm_rk = nm
    n0_rk = n0
    DO rk=1,4 !/ RungeKutta4-loop (RK) fNk0 etc are arrays of size (4,Ntot) containing RK-derivatives
       DO i=1,Ntot       !/ Clear derivative functions from previous run
          fNk0(rk,i)=0.0d0    !/ and set dummy variable N's to current state
          fNkp(rk,i)=0.0d0
          fNkm(rk,i)=0.0d0
       ENDDO
       fn0(rk) = 0.0d0
       fnp(rk) = 0.0d0
       fnm(rk) = 0.0d0
       !/ Calculate desired derivatives
       if (nuclFLAG==1)  call dNdt_nucl( fNk0(rk,:), fNkp(rk,:), fNkm(rk,:), fnm(rk), fnp(rk), &
                                         Jn0_scaled, Jnp_scaled, Jnm_scaled, i_Jn0, i_Jnp, i_Jnm, &
                                         Jnp, Jnm, Ntot)
       if (coagFLAG==1)  call dNdt_coag( Nk0_rk, Nkp_rk, Nkm_rk, fNk0(rk,:), fNkp(rk,:), fNkm(rk,:), &
                                         S, Vij, Kklp0, Kklpm, Kkl00, Kkl0m, Ntot)
       if (condFLAG==1)  call dNdt_cond( Nk0_rk, Nkp_rk, Nkm_rk, fNk0(rk,:), fNkp(rk,:), fNkm(rk,:), &
                                         v, v0, n0_rk, bk00, bk0p, bk0m, fn0(rk), Ntot)
       if (ioncFLAG==1)  call dNdt_ionc( Nk0_rk, Nkp_rk, Nkm_rk, fNk0(rk,:), fNkp(rk,:), fNkm(rk,:), &
                                         v, bkm0, bkp0, bkpm, bkmp, vcm, vcp, nm, np, Ntot)
       if (charFLAG==1)  call dNdt_char( Nk0_rk, Nkp_rk, Nkm_rk, fNk0(rk,:), fNkp(rk,:), fNkm(rk,:), &
                                         nm_rk, np_rk, fnm(rk), fnp(rk), bkm0, bkmp, bkp0, bkpm,  &
                                         Ntot)
       if (lossFLAG==1)  call dNdt_loss( Nk0_rk, Nkp_rk, Nkm_rk, fNk0(rk,:), fNkp(rk,:), fNkm(rk,:), &
                                         n0_rk, nm_rk, np_rk, fn0(rk), fnp(rk), fnm(rk), d, d_np,    &
                                         d_nm, d_n0, ggamma, lambda, d_loss, alpha,NL, bL0, bLp, bLm,& 
                                         kL0, kLp, kLm, Ntot)
       if (prodFLAG==1)  call dNdt_prod( fn0(rk), fnp(rk), fnm(rk), n0_rk, np_rk, nm_rk, P0, q, Ntot)
      !/ Calculate aerosol concentrations N's at relevant RK-timestep, and load in to RK-concentrations
       IF (rk < 4) THEN
          DO i=1,Ntot
             Nk0_rk(i) = Nk0(i) + (dt/rk_coef(rk)) * fNk0(rk,i)
             Nkp_rk(i) = Nkp(i) + (dt/rk_coef(rk)) * fNkp(rk,i)
             Nkm_rk(i) = Nkm(i) + (dt/rk_coef(rk)) * fNkm(rk,i)
          ENDDO
          !/ Either update or keep constant monomer and ion concentrations in RK depending on controlfile
          IF (n0_fix_flag == 0) THEN
             n0_rk = n0 + (dt/rk_coef(rk)) * fn0(rk)
          END IF
          IF (np_fix_flag == 0) THEN
             np_rk  = np  + (dt/rk_coef(rk)) * fnp(rk)  
          END IF 
          IF (nm_fix_flag == 0) THEN
             nm_rk  = nm  + (dt/rk_coef(rk)) * fnm(rk)
          END IF
       END IF
    ENDDO
    DO i=1,Ntot !/ Calculate final new aerosol concentrations	
       Nk0(i) = Nk0(i) + dt * ( fNk0(1,i) + 2.0d0*fNk0(2,i) + 2.0d0*fNk0(3,i) + fNk0(4,i) ) / 6.0d0
       Nkp(i) = Nkp(i) + dt * ( fNkp(1,i) + 2.0d0*fNkp(2,i) + 2.0d0*fNkp(3,i) + fNkp(4,i) ) / 6.0d0
       Nkm(i) = Nkm(i) + dt * ( fNkm(1,i) + 2.0d0*fNkm(2,i) + 2.0d0*fNkm(3,i) + fNkm(4,i) ) / 6.0d0
    ENDDO
    !/ Either update or keep constan monomert and ion concentrations depending on controlfile
    IF (n0_fix_flag == 0) THEN 
       n0 = n0 + dt * ( fn0(1) + 2.0d0*fn0(2) + 2.0d0*fn0(3) + fn0(4))  / 6.0d0
    END IF
    IF (np_fix_flag == 0) THEN
       np  = np +  dt * ( fnp(1)  + 2.0d0*fnp(2)  + 2.0d0*fnp(3)  + fnp(4) )  / 6.0d0
    END IF
    IF (nm_fix_flag == 0) THEN 
       nm  = nm +  dt * ( fnm(1)  + 2.0d0*fnm(2)  + 2.0d0*fnm(3)  + fnm(4) )  / 6.0d0
    END IF
 
!/ PRINT SNAPSHOT
   IF (mod(tt,ts) < dt) call print_snapshot( Nk0, Nkp, Nkm, n0, np, nm, tt, Ntot)

 ENDDO
!//////////////////////////                                                                              //
!/ TIME INTEGRATION LOOP //                              END                                             //
!//////////////////////////////////////////////////////////////////////////////////////////////////////////

 90 CONTINUE

! Deallocate arrays
 deallocate(Vij)
 deallocate(fNk0,fNkp,fNkm)
 deallocate(Nkp, Nk0, Nkm, v, d, bkp0, bkpm, bk0p, bk00, bk0m, bkmp, bkm0)
 deallocate(Nk0_rk, Nkp_rk, Nkm_rk, kL0,kLp,kLm)
 deallocate(S, Kklp0,  Kklpm, Kkl00, Kkl0m)

END
