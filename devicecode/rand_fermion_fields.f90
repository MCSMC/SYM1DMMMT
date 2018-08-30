!@TODO use here a better CURAND random number generator.
module rand_fermion_fields
   implicit none
   contains
SUBROUTINE generate_rand_vect(pf)

   use compiletimeconstants
   implicit none

    double complex pf(1:nmat,1:nmat,1:nspin,-(nmargin-1):nsite+nmargin,1:npf)
    !$acc declare device_resident(pf)
    double complex pfhost(1:nmat,1:nmat,1:nspin,-(nmargin-1):nsite+nmargin,1:npf)
    double precision r1,r2
    integer ipf,isite,ispin,jmat,imat

    do ipf=1,npf
        do isite=1,nsite
            do ispin=1,nspin
                do jmat=1,nmat
                    do imat=1,nmat
                        call BoxMuller(r1,r2)
                        pfhost(imat,jmat,ispin,isite,ipf)=&
                            &(dcmplx(r1)+dcmplx(r2)*(0D0,1D0))/dcmplx(dsqrt(2d0))
                    end do
                end do
            end do
        end do
    end do
    !$acc data copyin(pfhost)
    !$acc kernels
    pf=pfhost
  !$acc end kernels
  !$acc end data
END  SUBROUTINE generate_rand_vect

SUBROUTINE Adjust_margin_and_bc_pf_device(pf,nbc)

    use compiletimeconstants
    use dirac_operator
    implicit none

    !***** input *****
    integer nbc
    !***** input & output *****
    double complex pf(1:nmat,1:nmat,1:nspin,-(nmargin-1):nsite+nmargin,1:npf)
    !$acc declare present(nbc)
    !$acc declare device_resident(pf)
    !*******************
    integer isite,imat,jmat,ispin,ipf

    do ipf=1,npf
        call set_boundary_device(nbc,pf(:,:,:,:,ipf))
    end do

    return

END SUBROUTINE Adjust_margin_and_bc_pf_device


!pseudo fermion
!phi is generated by Gaussian weight e^{-Tr(phi^dag*phi)}, 
!pf=D^{-1/8}*phi
!traceless condition sum_i phi_ii= sum_i pf_ii=0 is taken into account. 
SUBROUTINE generate_pseudo_fermion_SUN_device(pf,xmat,alpha,phase,&
    &Gam123,acoeff_pf,bcoeff_pf,max_err,max_iteration,iteration,&
    &nbc,nbmn,temperature,flux,info)
  
    use cgm_solver
    !use mtmod !Mersenne twistor
    implicit none
  
  include '../staticparameters.f90'
    !***** input *****
    integer nbc,nbmn,info
    double precision temperature,flux
    double complex xmat(1:nmat,1:nmat,1:ndim,-(nmargin-1):nsite+nmargin)
    double precision alpha(1:nmat)
    integer max_iteration
    double precision acoeff_pf(0:nremez_pf)
    double precision bcoeff_pf(1:nremez_pf)!bcoeff(1) is the smallest.
    double precision max_err
    double complex :: Gam123(1:nspin,1:nspin)
    double complex :: phase(1:nmat,1:nmat,1:2)
    !***** output *****
    double complex pf(1:nmat,1:nmat,1:nspin,-(nmargin-1):nsite+nmargin,1:npf)
    integer iteration
    !*********************
    double complex phi(1:nmat,1:nmat,1:nspin,-(nmargin-1):nsite+nmargin,1:npf)
    double complex chi(1:nmat,1:nmat,1:nspin,&
        &-(nmargin-1):nsite+nmargin,1:nremez_pf,1:npf)
    !$acc declare present(xmat,alpha,nbc,nbmn,acoeff_pf,bcoeff_pf,flux,temperature)
    !$acc declare device_resident(phi,chi,Gam123,phase,pf)
    double complex trace,tmp
    integer imat,jmat
    integer ispin
    integer isite
    integer iremez
    integer i
    integer ipf

    !**************************************************************************
    !**** we must be careful about the normalization of the Gaussian term. ****
    !**************************************************************************
    call generate_rand_vect(phi)

    !$acc kernels
    do ipf=1,npf
        !****************************
        !*** traceless projection ***
        !****************************
        do isite=1,nsite
            do ispin=1,nspin
                trace=(0d0,0d0)
                do imat=1,nmat
                    trace=trace+phi(imat,imat,ispin,isite,ipf)
                end do
                trace=trace/dcmplx(nmat)
                do imat=1,nmat
                    phi(imat,imat,ispin,isite,ipf)=phi(imat,imat,ispin,isite,ipf)-trace
                end do
            end do
        end do
    end do
    !$acc end kernels
    !*********************************
    !*** adjust the margin and b.c.***
    !*********************************
    call Adjust_margin_and_bc_pf_device(phi,nbc)
  
    !**************************
    !*** pf = D^{-1/8}*phi  ***
    !**************************

    call cgm_solver_device(nremez_pf,bcoeff_pf,nbmn,nbc,temperature,&
        max_err,max_iteration,xmat,phase,Gam123,phi,chi,info,iteration)

    if(rhmc_verbose.EQ.1) then
        !$acc kernels
        tmp=Sum(chi)
        !$acc end kernels
        print*, "pseudoferm generation  dev chi ", tmp
           !$acc kernels
        tmp=Sum(phi)
        !$acc end kernels
        print*, "pseudoferm generation  dev phi ", tmp
    end if

    !$acc kernels
    do ipf=1,npf
        do isite=-(nmargin-1),nsite+nmargin
            do ispin=1,nspin
                do jmat=1,nmat
                    do imat=1,nmat
                        pf(imat,jmat,ispin,isite,ipf)=&
                            dcmplx(acoeff_pf(0))*phi(imat,jmat,ispin,isite,ipf)
                    end do
                end do
            end do
        end do
        do iremez=1,nremez_pf
            do isite=-(nmargin-1),nsite+nmargin
                do ispin=1,nspin
                    do jmat=1,nmat
                        do imat=1,nmat
                    
                            pf(imat,jmat,ispin,isite,ipf)=pf(imat,jmat,ispin,isite,ipf)&
                                &+dcmplx(acoeff_pf(iremez))&
                                &*chi(imat,jmat,ispin,isite,iremez,ipf)
                        end do
                    end do
                end do
            end do
        end do
    end do
    !$acc end kernels

    return
  
END SUBROUTINE Generate_pseudo_fermion_SUN_device

! At the moment this just redirects to host code.
! It is not easy to implement the parallel rnd generator.
! @TODO: Add efficient cuda random number generator.
SUBROUTINE Generate_Momenta_device(P_xmat,P_alpha)

  use compiletimeconstants
  implicit none

  integer imat,jmat,idim,isite

  double complex P_xmat(1:nmat,1:nmat,1:ndim,1:nsite)
  double precision P_alpha(1:nmat)
  !$acc declare device_resident(P_xmat,P_alpha)
  double complex P_xmat_host(1:nmat,1:nmat,1:ndim,1:nsite)
  double precision P_alpha_host(1:nmat)

  call Generate_P_xmat(P_xmat_host)
  call Generate_P_alpha(P_alpha_host)
  !$acc data copyin(P_xmat_host,P_alpha_host)
  !$acc kernels
  P_xmat=P_xmat_host
  P_alpha=P_alpha_host
  !$acc end kernels
  !$acc end data

  return

END SUBROUTINE Generate_Momenta_device
end module rand_fermion_fields
