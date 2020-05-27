MODULE shared_data
     use precision_type
     integer,allocatable  :: Vij(:,:) 
     integer              :: charFLAG,coagFLAG,nuclFLAG,condFLAG,lossFLAG,loadFLAG,ioncFLAG,prodFLAG
     integer              :: np_fix_flag, nm_fix_flag, n0_fix_flag, i_Jn0, i_Jnp, i_Jnm, Ntot
     real(dp)             :: d_np, d_nm, d_n0, NL, d_loss, bL0, bLp, bLm, Jn0_scaled, Jnp_scaled, Jnm_scaled
     real(dp)             :: ggamma, lambda, v0, vcm, vcp, q, alpha, P0, Jnp, Jnm
     real(dp),allocatable :: v(:),d(:),bkp0(:),bkpm(:),bk0p(:),bk00(:),bk0m(:),bkmp(:),bkm0(:)
     real(dp),allocatable :: S(:,:),Kklp0(:,:),Kklpm(:,:),Kkl00(:,:),Kkl0m(:,:), kL0(:), kLp(:), kLm(:)
END MODULE shared_data
