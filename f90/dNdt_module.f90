MODULE dNdt_module
!///////////////////////////////////////////
!/ Ions are treated as clusters with mass //
!///////////////////////////////////////////

 IMPLICIT NONE

 CONTAINS

!//////////////////////
!/ NUCLEATION MODULE //
!//////////////////////

 SUBROUTINE dNdt_nucl( dNk0dt, dNkpdt, dNkmdt, dnmdt, dnpdt)

    USE precision_type
    USE shared_data, ONLY : i_Jn0, i_Jnp, i_Jnm, Jn0_scaled, Jnp_scaled, Jnm_scaled, Jnp, Jnm, Ntot 

    real(dp) :: dnpdt, dnmdt
    real(dp),dimension(Ntot) :: dNk0dt, dNkpdt, dNkmdt

    dNk0dt(i_Jn0) = dNk0dt(i_Jn0) + Jn0_scaled
    dNkpdt(i_Jnp) = dNkpdt(i_Jnp) + Jnp_scaled
    dNkmdt(i_Jnm) = dNkmdt(i_Jnm) + Jnm_scaled
    dnpdt = dnpdt - Jnp
    dnmdt = dnmdt - Jnm

 END SUBROUTINE dNdt_nucl

!///////////////////////
!/ COAGULATION MODULE //
!///////////////////////

 SUBROUTINE dNdt_coag( Nk0, Nkp, Nkm, dNk0dt, dNkpdt, dNkmdt)
    USE precision_type
    USE shared_data, ONLY          : S, Vij, Kklp0, Kklpm, Kkl00, Kkl0m, Ntot
    integer                       :: i, j
    real(dp)                      :: Ip0, Ipm, I00, I0m
    real(dp),dimension(Ntot)      :: Nk0, Nkp, Nkm, dNk0dt, dNkpdt,dNkmdt

    DO i=1,Ntot
       DO j=1,Ntot 
          Ip0 = Kklp0(i,j) * Nkp(i) * Nk0(j)
          Ipm = Kklpm(i,j) * Nkp(i) * Nkm(j)
          I00 = Kkl00(i,j) * Nk0(i) * Nk0(j)
          I0m = Kkl0m(i,j) * Nk0(i) * Nkm(j) 
          dNk0dt(Vij(i,j)) = dNk0dt(Vij(i,j)) + S(i,j) * (I00/2.0d0+Ipm)
          dNkpdt(Vij(i,j)) = dNkpdt(Vij(i,j)) + S(i,j) * Ip0
          dNkmdt(Vij(i,j)) = dNkmdt(Vij(i,j)) + S(i,j) * I0m
          IF (Vij(i,j) < Ntot) THEN
             dNk0dt(Vij(i,j)+1) = dNk0dt(Vij(i,j)+1) + (1.0d0-S(i,j))*(I00/2.0d0+Ipm)
             dNkpdt(Vij(i,j)+1) = dNkpdt(Vij(i,j)+1) + (1.0d0-S(i,j))*Ip0
             dNkmdt(Vij(i,j)+1) = dNkmdt(Vij(i,j)+1) + (1.0d0-S(i,j))*I0m
          ENDIF
          dNk0dt(i) = dNk0dt(i) - I00/2.0d0 - I0m
          dNkpdt(i) = dNkpdt(i) - Ip0       - Ipm
          dNk0dt(j) = dNk0dt(j) - Ip0       - I00/2.0d0
          dNkmdt(j) = dNkmdt(j) - Ipm       - I0m 
      ENDDO
    ENDDO
 END SUBROUTINE dNdt_coag

!////////////////////////
!/ CONDENSATION MODULE //
!////////////////////////

 SUBROUTINE dNdt_cond( Nk0, Nkp, Nkm, dNk0dt, dNkpdt, dNkmdt, n0, dn0dt)
    USE precision_type
    USE shared_data, ONLY     : bk00, bk0p, bk0m, v, v0, Ntot
    integer                  :: i
    real(dp)                 :: n0, dn0dt, temp, vmaxp1
    real(dp),dimension(Ntot) :: Nk0, Nkp, Nkm, dNk0dt, dNkpdt,dNkmdt

    dNk0dt(1) = dNk0dt(1) - n0 * bk00(1)*Nk0(1)/(v(2)-v(1)) * v0
    dNkpdt(1) = dNkpdt(1) - n0 * bk0p(1)*Nkp(1)/(v(2)-v(1)) * v0
    dNkmdt(1) = dNkmdt(1) - n0 * bk0m(1)*Nkm(1)/(v(2)-v(1)) * v0
    DO i=2,Ntot-1
       dNk0dt(i) = dNk0dt(i) + n0 * (bk00(i-1)*Nk0(i-1)/(v(i)-v(i-1)) - bk00(i)*Nk0(i)/(v(i+1)-v(i))) * v0
       dNkpdt(i) = dNkpdt(i) + n0 * (bk0p(i-1)*Nkp(i-1)/(v(i)-v(i-1)) - bk0p(i)*Nkp(i)/(v(i+1)-v(i))) * v0
       dNkmdt(i) = dNkmdt(i) + n0 * (bk0m(i-1)*Nkm(i-1)/(v(i)-v(i-1)) - bk0m(i)*Nkm(i)/(v(i+1)-v(i))) * v0
    ENDDO
    vmaxp1 = v(Ntot)*(v(Ntot)/v(Ntot-1))
 
    dNk0dt(Ntot) = dNk0dt(Ntot) + n0 * (bk00(Ntot-1)*Nk0(Ntot-1)/(v(Ntot)  -v(Ntot-1)) &
                                - bk00(Ntot)  *Nk0(Ntot)  /(vmaxp1-v(Ntot)))   * v0  
    dNkpdt(Ntot) = dNkpdt(Ntot) + n0 * (bk0p(Ntot-1)*Nkp(Ntot-1)/(v(Ntot)  -v(Ntot-1)) &
                                - bk0p(Ntot)  *Nkp(Ntot)  /(vmaxp1-v(Ntot)))   * v0  
    dNkmdt(Ntot) = dNkmdt(Ntot) + n0 * (bk0m(Ntot-1)*Nkm(Ntot-1)/(v(Ntot)  -v(Ntot-1)) &
                                - bk0m(Ntot)  *Nkm(Ntot)  /(vmaxp1-v(Ntot)))   * v0  

    temp = 0.0d0
    DO i=1,Ntot
       temp = temp + Nk0(i)*bk00(i) + Nkp(i)*bk0p(i) + Nkm(i)*bk0m(i)
    ENDDO
    dn0dt = dn0dt - n0*temp

 END SUBROUTINE dNdt_cond


!////////////////////////////////////
!/ ION-INDUCED CONDENSATION MODULE //
!////////////////////////////////////

 SUBROUTINE dNdt_ionc( Nk0, Nkp, Nkm, dNk0dt, dNkpdt, dNkmdt, nm, np)
    USE precision_type
    USE shared_data, ONLY     : P0, v0, vcm, vcp, v, bkm0, bkp0, bkpm, bkmp, Ntot
    integer                  :: i
    real(dp)                 :: np, nm, temp, vm, vp, vmaxp1
    real(dp),dimension(Ntot) :: Nk0, Nkp, Nkm, dNk0dt, dNkpdt, dNkmdt
    vm = vcm
    vp = vcp

    dNk0dt(1) = dNk0dt(1) - nm * bkmp(1)*Nkp(1)/(v(2)-v(1)) * vm   - np * bkpm(1)*Nkm(1)/(v(2)-v(1)) * vp
    dNkpdt(1) = dNkpdt(1) - np * bkp0(1)*Nk0(1)/(v(2)-v(1)) * vp
    dNkmdt(1) = dNkmdt(1) - nm * bkm0(1)*Nk0(1)/(v(2)-v(1)) * vm
    DO i=2,Ntot-1
       dNk0dt(i) = dNk0dt(i) +  nm * (bkmp(i-1)*Nkp(i-1)/(v(i)-v(i-1)) - bkmp(i)*Nkp(i)/(v(i+1)-v(i))) * vm   &
                             +  np * (bkpm(i-1)*Nkm(i-1)/(v(i)-v(i-1)) - bkpm(i)*Nkm(i)/(v(i+1)-v(i))) * vp
       dNkpdt(i) = dNkpdt(i) +  np * (bkp0(i-1)*Nk0(i-1)/(v(i)-v(i-1)) - bkp0(i)*Nk0(i)/(v(i+1)-v(i))) * vp
       dNkmdt(i) = dNkmdt(i) +  nm * (bkm0(i-1)*Nk0(i-1)/(v(i)-v(i-1)) - bkm0(i)*Nk0(i)/(v(i+1)-v(i))) * vm
    ENDDO
    vmaxp1 = v(Ntot)*(v(Ntot)/v(Ntot-1))

    dNk0dt(Ntot) = dNk0dt(Ntot) +  nm * (bkmp(Ntot-1)*Nkp(Ntot-1)/(v(Ntot)  -v(Ntot-1)) &
                                - bkmp(Ntot)  *Nkp(Ntot)  /(vmaxp1-v(Ntot)))   * vm   &
                                +  np * (bkpm(Ntot-1)*Nkm(Ntot-1)/(v(Ntot)  -v(Ntot-1)) &
                                - bkpm(Ntot)  *Nkm(Ntot)  /(vmaxp1-v(Ntot)))   * vp
    dNkpdt(Ntot) = dNkpdt(Ntot) +  np * (bkp0(Ntot-1)*Nk0(Ntot-1)/(v(Ntot)  -v(Ntot-1)) &
                                - bkp0(Ntot)  *Nk0(Ntot)  /(vmaxp1-v(Ntot)))   * vp
    dNkmdt(Ntot) = dNkmdt(Ntot) +  nm * (bkm0(Ntot-1)*Nk0(Ntot-1)/(v(Ntot)  -v(Ntot-1)) &
                                - bkm0(Ntot)  *Nk0(Ntot)  /(vmaxp1-v(Ntot)))   * vm
    ! Note, change in np and nm is calculated in the dNdt_char subroutine, where also
    ! the switching of charge is taken into account. If the flag for this subroutine is 
    ! switched on in the controlfile, then dNdt_char flag must also be on.

 END SUBROUTINE dNdt_ionc

!///////////////////////////
!/ CHARGE EXCHANGE MODULE //
!///////////////////////////

 SUBROUTINE dNdt_char( Nk0, Nkp, Nkm, dNk0dt, dNkpdt, dNkmdt, nm, np, dnmdt,dnpdt)
    USE precision_type
    USE shared_data, ONLY     : bkm0, bkmp, bkp0, bkpm, Ntot
    integer                  :: i
    real(dp)                 :: nm, np, dnmdt, dnpdt, sum1, sum2
    real(dp),dimension(Ntot) :: Nk0,Nkp,Nkm,dNk0dt,dNkpdt,dNkmdt

    sum1 = 0.0d0
    sum2 = 0.0d0
    DO i=1,Ntot
       sum1 = sum1 + bkm0(i)*Nk0(i) + bkmp(i)*Nkp(i)
       sum2 = sum2 + bkp0(i)*Nk0(i) + bkpm(i)*Nkm(i)
    ENDDO
    dnmdt = dnmdt - nm*sum1
    dnpdt = dnpdt - np*sum2
    DO i=1,Ntot
       dNk0dt(i) = dNk0dt(i) + np*bkpm(i)*Nkm(i) + nm*bkmp(i)*Nkp(i) &
                             - np*bkp0(i)*Nk0(i) - nm*bkm0(i)*Nk0(i) 
       dNkpdt(i) = dNkpdt(i) + np*bkp0(i)*Nk0(i) - nm*bkmp(i)*Nkp(i)
       dNkmdt(i) = dNkmdt(i) + nm*bkm0(i)*Nk0(i) - np*bkpm(i)*Nkm(i)
    ENDDO
 END SUBROUTINE dNdt_char


!///////////////////////////
!/      LOSS MODULE       //
!///////////////////////////

 SUBROUTINE dNdt_loss( Nk0, Nkp, Nkm,  dNk0dt, dNkpdt, dNkmdt, n0,  nm, np, dn0dt, dnpdt, dnmdt)

    USE precision_type
    USE shared_data, ONLY     : d_np, d_nm, d_n0, ggamma, lambda, d_loss, alpha,d, kL0, &
                                kLp, kLm,d_np, NL, bL0, bLp, bLm, Ntot
    integer                  :: i
    real(dp)                 :: n0, nm, np, dn0dt, dnpdt, dnmdt
    real(dp),dimension(Ntot) :: Nk0, Nkp, Nkm, dNk0dt, dNkpdt, dNkmdt 

    ! Wall loss
    DO i=1,Ntot
       dNk0dt(i) = dNk0dt(i) - lambda * ( 1.0d9*d(i) / d_loss )**(-ggamma) * Nk0(i)
       dNkmdt(i) = dNkmdt(i) - lambda * ( 1.0d9*d(i) / d_loss )**(-ggamma) * Nkm(i)
       dNkpdt(i) = dNkpdt(i) - lambda * ( 1.0d9*d(i) / d_loss )**(-ggamma) * Nkp(i)
    ENDDO
    dn0dt = dn0dt - lambda * (d_n0 * 1.0d9 / d_loss)**(-ggamma) * n0 
    dnmdt = dnmdt - lambda * (d_np * 1.0d9 / d_loss)**(-ggamma) * nm - alpha*nm*np
    dnpdt = dnpdt - lambda * (d_nm * 1.0d9 / d_loss)**(-ggamma) * np - alpha*nm*np

    ! Loss to large particles
    DO i=1,Ntot
       dNk0dt(i) = dNk0dt(i) - Nk0(i) * kL0(i) * NL
       dNkmdt(i) = dNkmdt(i) - Nkm(i) * kLm(i) * NL
       dNkpdt(i) = dNkpdt(i) - Nkp(i) * kLp(i) * NL
    ENDDO
    dn0dt = dn0dt - n0 * bL0 * NL
    dnmdt = dnmdt - nm * bLm * NL
    dnpdt = dnpdt - np * bLp * NL

 END SUBROUTINE dNdt_loss


!///////////////////////////////
!/ MONOMER PRODUCTION MODULE  //
!///////////////////////////////

 SUBROUTINE dNdt_prod(dn0dt, dnpdt, dnmdt)

    USE precision_type
    USE shared_data, ONLY : P0, q

    real(dp) :: dn0dt, dnpdt, dnmdt
    dn0dt = dn0dt + P0 
    dnpdt = dnpdt + q
    dnmdt = dnmdt + q
    
 END SUBROUTINE dNdt_prod


END MODULE dNdt_module
