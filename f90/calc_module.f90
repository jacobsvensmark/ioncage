MODULE calc_module

 IMPLICIT NONE

 CONTAINS

 double precision FUNCTION brownian_kernel( ri, rj, rho_p, T, eta, lambda_a)
 !/////////////////////////////////////////////////////////////////////////////////////////////////////
 !/ Function returning the brownian coagulation coeff. Eq. 1 as                                      //
 !/ given in Jacobsen and Seinfeld 2004 Atmospheric Environment 38 (2004) 1839–1850 (doi:10.1016).   // 
 !/////////////////////////////////////////////////////////////////////////////////////////////////////
 !/ Input(Radius_i, Radius_j,Aerosol Density, Temperature, Viscosity, Mean Free Path of Air)         //
 !/////////////////////////////////////////////////////////////////////////////////////////////////////
    USE precision_type
    IMPLICIT NONE
    real(dp) :: ri,rj,rho_p,T,eta,lambda_a
    real(dp) :: pi,k,Mi,Mj,vi,vj,Kni,Knj,Am,Bm,Cm,Dmi,Dmj,lambdai,lambdaj,deltai,deltaj, T1,T2,T3,beta
    pi = 3.141592653589793
    k  = 1.380648813d-23
    Mi = 4.0d0 * pi * (ri**3) * rho_p / 3.0d0
    Mj = 4.0d0 * pi * (rj**3) * rho_p / 3.0d0
    vi = sqrt(8.0d0*k*T/(pi*Mi))
    vj = sqrt(8.0d0*k*T/(pi*Mj))
    Kni = lambda_a / ri
    Knj = lambda_a / rj
    write(*,"(7ES16.6E3)") Mi,Mj,vi,vj,Kni,Knj,T
    Am = 0.864d0
    Bm = 0.29d0
    Cm = 1.25d0
    Dmi = k*T/(6.0d0*pi*ri * eta) * (1.0d0+Kni*(Am+Bm*exp(-Cm/Kni)))
    Dmj = k*T/(6.0d0*pi*rj * eta) * (1.0d0+Knj*(Am+Bm*exp(-Cm/Knj)))
    lambdai = 8.0d0*Dmi/(pi*vi)
    lambdaj = 8.0d0*Dmj/(pi*vj)
    deltai = ((2.0d0*ri+lambdai)**3 - (4.0d0*(ri**2) + (lambdai**2))**(3.0d0/2.0d0) ) / (6.0d0*ri*lambdai) - 2.0d0 * ri
    deltaj = ((2.0d0*rj+lambdaj)**3 - (4.0d0*(rj**2) + (lambdaj**2))**(3.0d0/2.0d0) ) / (6.0d0*rj*lambdaj) - 2.0d0 * rj
    T1 = 4.0d0 * pi * (ri+rj) * (Dmi + Dmj)
    T3 = 4.0d0 * (Dmi + Dmj)/(sqrt(vi*vi + vj*vj) * (ri+rj))
    T2 = (ri+rj) / (ri+rj+sqrt((deltai**2)+(deltaj**2)))
    beta  = T1/(T2+T3)
    brownian_kernel = beta

    RETURN
 END FUNCTION brownian_kernel

END MODULE calc_module
