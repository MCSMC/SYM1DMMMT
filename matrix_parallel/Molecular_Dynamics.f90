! alpha -> alpha + P_alpha*dtau_alpha
! xmat -> xmat + P_xmat*dtau_xmat
! P_alpha -> P_alpha - delh_alpha*dtau_alpha
! P_xmat -> P_xmat - delh_xmat*dtau_xmat
! delh_xmat(imat,jmat)=dS/dxmat(jmat,imat)

subroutine Molecular_Dynamics(nbc,temperature,&
    &ntau,dtau_xmat,dtau_alpha,xmat,alpha,P_xmat,P_alpha,&
    &acoeff_md,bcoeff_md,pf,max_iteration,max_err,iteration,&
    &gamma10d,g_alpha,g_R,RCUT,acceleration,&
    &nbmn,flux,info,nsmear,s_smear,ngauge,purebosonic)

    implicit none

  include 'mpif.h'
  include 'size_parallel.h'
  include '../Fourier.inc'
  include '../unit_number.inc'
    integer nbc,nbmn,nsmear,ngauge,purebosonic
    integer max_iteration
    double precision max_err
    double precision g_alpha
    integer ntau
    doubleprecision temperature,flux
    double precision acceleration(1:nsite_local)
    double precision dtau_xmat,dtau_alpha
    double precision g_R,RCUT
    double complex pf(1:nmat_block,1:nmat_block,&
        1:nspin,-(nmargin-1):nsite_local+nmargin)
    double complex GAMMA10d(1:ndim,1:nspin,1:nspin)
    double precision s_smear
    !***** input & output *****
    double complex xmat(1:nmat_block,1:nmat_block,1:ndim,&
        &-(nmargin-1):nsite_local+nmargin)
    double precision alpha(1:nmat_block*nblock)
    double complex P_xmat(1:nmat_block,1:nmat_block,1:ndim,1:nsite_local)
    double precision P_alpha(1:nmat_block*nblock)
    double precision acoeff_md(0:nremez_md),bcoeff_md(1:nremez_md)
    !***** output *****
    integer info!0 -> OK, 1 -> err in CG solver.
    integer iteration
    !**************************
    double complex xmat_smeared(1:nmat_block,1:nmat_block,1:ndim,&
        &-(nmargin-1):nsite_local+nmargin)
    double complex xmat_mom(1:nmat_block,1:nmat_block,1:ndim,1:nsite_local)
    double complex P_xmat_mom(1:nmat_block,1:nmat_block,1:ndim,1:nsite_local)
    double complex delh_xmat(1:nmat_block,1:nmat_block,1:ndim,1:nsite_local)
    double complex delh_xmat_mom(1:nmat_block,1:nmat_block,1:ndim,1:nsite_local)
    double precision delh_alpha(1:nmat_block*nblock)
    double precision delh_alpha_local(1:nmat_block*nblock)
    integer imat,jmat
    integer imom
    integer idim
    integer step
    double complex Chi(1:nremez_md,1:nmat_block,1:nmat_block,&
        1:nspin,-(nmargin-1):nsite_local+nmargin)
    !***** for MPI *****
    integer ierr,myrank
    integer iblock,jblock,isublat
   
    call MPI_COMM_RANK(MPI_COMM_WORLD,MYRANK, IERR)
    call who_am_i(myrank,isublat,iblock,jblock)
    !*******************************
    !*** first step of leap frog ***
    !*******************************
    !Fourier transform from P_xmat to P_xmat_mom
    call Fourier_transform_P_xmat(P_xmat,P_xmat_mom,myrank,x2p)
    !Fourier transform from xmat to xmat_mom
    call Fourier_transform_xmat(xmat,xmat_mom,myrank,x2p)
    !move xmat_mom and alpha (1/2)-step forward.
    do imom=1,nsite_local
        do idim=1,ndim
            !$omp parallel
            !$omp do
            do imat=1,nmat_block
                do jmat=1,nmat_block
                    xmat_mom(imat,jmat,idim,imom)=&
                        &xmat_mom(imat,jmat,idim,imom)&
                        &+P_xmat_mom(imat,jmat,idim,imom)&
                        &*dcmplx(0.5d0*dtau_xmat*acceleration(imom))
                end do
            end do
        !$omp end do
        !$omp end parallel
        end do
    end do
    alpha=alpha+P_alpha*0.5d0*dtau_alpha
    !Fourier transform from xmat_mom to xmat
    call Fourier_transform_xmat(xmat,xmat_mom,myrank,p2x)

    !************************************
    !*** second step,...,Ntau-th step ***
    !************************************
    info=0
    !info=0 -> CG solver correctly worked.
    !info=1 -> CG solver returned error.
    step=1
    do while ((info.EQ.0).AND.(step.LT.ntau))
        step=step+1
        !adjust the margin.
        call Adjust_margin_xmat(xmat,myrank)
        !smear
        if(purebosonic.eq.0) then
            call smearing_xmat(xmat,xmat_smeared,s_smear,myrank,nsmear)
            !Calculate Chi_k = (D+bcoeff_md(k))^{-1}*pf by using multi-mass CG-solver
            call solver_biCGm(nbc,nremez_md,bcoeff_md,temperature,&
                &xmat_smeared,alpha,pf,chi,GAMMA10d,max_err,max_iteration,iteration,myrank,&
                &nbmn,flux,info)
        end if
        if(myrank.EQ.0)then
            !Take CG_log
            write(unit_CG_log,*)"molecular evolution",iteration
        end if

        !calculate the force term in the coordinate space.
        !delh_xmat=dH/dX, delh_alpha=dH/(d alpha)
        call Calc_Force(delh_xmat,delh_alpha_local,xmat,alpha,&
            &chi,GAMMA10d,g_alpha,g_R,RCUT,nbmn,temperature,flux,acoeff_md,&
            &xmat_smeared,nsmear,s_smear,ngauge,purebosonic)
        call MPI_Allreduce(delh_alpha_local(1),delh_alpha(1),nmat_block*nblock,&
            &MPI_DOUBLE_PRECISION,&
            &MPI_SUM,MPI_COMM_WORLD,IERR)
        !Fourier transform from delh_xmat to delh_xmat_mom,
        !i.e. convert the force(delh_xmat) to the Fourier mode (delh_xmat_mom).
        call Fourier_transform_P_xmat(delh_xmat,delh_xmat_mom,myrank,x2p)
        !move P_xmat_mom and P_alpha one step forward.
        P_alpha=P_alpha-delh_alpha*dtau_alpha
        do imom=1,nsite_local
            do idim=1,ndim
                !$omp parallel
                !$omp do
                do imat=1,nmat_block
                    do jmat=1,nmat_block
                        P_xmat_mom(imat,jmat,idim,imom)=&
                            &P_xmat_mom(imat,jmat,idim,imom)&
                            &-delh_xmat_mom(imat,jmat,idim,imom)&
                            &*dcmplx(dtau_xmat*acceleration(imom))
                    end do
                end do
            !$omp end do
            !$omp end parallel
            end do
        end do
        !move xmat_mom and alpha one step forward.
        do imom=1,nsite_local
            do idim=1,ndim
                !$omp parallel
                !$omp do
                do imat=1,nmat_block
                    do jmat=1,nmat_block
                        xmat_mom(imat,jmat,idim,imom)=&
                            &xmat_mom(imat,jmat,idim,imom)&
                            &+P_xmat_mom(imat,jmat,idim,imom)&
                            &*dcmplx(dtau_xmat*acceleration(imom))
                    end do
                end do
            !$omp end do
            !$omp end parallel
            end do
        end do
        alpha=alpha+P_alpha*dtau_alpha
        !Fourier transform from xmat_mom to xmat
        call Fourier_transform_xmat(xmat,xmat_mom,myrank,p2x)
    end do
    !*****************
    !*** last step ***
    !*****************
    !info=0 -> so far so good, let's go to the last step.
    !info=1 -> error happened somewhere,
    !      just skip this last step and return info=1.
    if(info.EQ.0)then
        !adjust margin.
        call Adjust_margin_xmat(xmat,myrank)
        !smear
        if(purebosonic.eq.0) then
            call smearing_xmat(xmat,xmat_smeared,s_smear,myrank,nsmear)
            !Calculate Chi_k = (D+bcoeff_md(k))^{-1}*pf by using multi-mass CG-solver
            call solver_biCGm(nbc,nremez_md,bcoeff_md,temperature,&
                &xmat_smeared,alpha,pf,chi,GAMMA10d,max_err,max_iteration,iteration,myrank,&
                &nbmn,flux,info)
        end if
        if(myrank.EQ.0)then
            !Take CG_log
            write(unit_CG_log,*)"molecular evolution",iteration
        end if
        !calculate the force term in the coordinate space.
        !delh_xmat=dH/dX, delh_alpha=dH/(d alpha)
        call Calc_Force(delh_xmat,delh_alpha_local,xmat,alpha,&
            &chi,GAMMA10d,g_alpha,g_R,RCUT,nbmn,temperature,flux,acoeff_md,&
            &xmat_smeared,nsmear,s_smear,ngauge,purebosonic)
        call MPI_Allreduce(delh_alpha_local(1),delh_alpha(1),nmat_block*nblock,MPI_DOUBLE_PRECISION,&
            &MPI_SUM,MPI_COMM_WORLD,IERR)
        !Fourier transform from delh_xmat to delh_xmat_mom,
        !i.e. convert the force(delh_xmat) to the Fourier mode (delh_xmat_mom).
        call Fourier_transform_P_xmat(delh_xmat,delh_xmat_mom,myrank,x2p)
        !move P_xmat_mom and P_alpha one step forward.
        P_alpha=P_alpha-delh_alpha*dtau_alpha
        do imom=1,nsite_local
            do idim=1,ndim
                !$omp parallel
                !$omp do
                do imat=1,nmat_block
                    do jmat=1,nmat_block
                        P_xmat_mom(imat,jmat,idim,imom)=&
                            &P_xmat_mom(imat,jmat,idim,imom)&
                            &-delh_xmat_mom(imat,jmat,idim,imom)&
                            &*dcmplx(dtau_xmat*acceleration(imom))
                    end do
                end do
            !$omp end do
            !$omp end parallel
            end do
        end do
        !move xmat_mom and alpha (1/2)-step forward.
        do imom=1,nsite_local
            do idim=1,ndim
                !$omp parallel
                !$omp do
                do imat=1,nmat_block
                    do jmat=1,nmat_block
                        xmat_mom(imat,jmat,idim,imom)=&
                            &xmat_mom(imat,jmat,idim,imom)&
                            &+P_xmat_mom(imat,jmat,idim,imom)&
                            &*dcmplx(0.5d0*dtau_xmat*acceleration(imom))
                    end do
                end do
            !$omp end do
            !$omp end parallel
            end do
        end do
        alpha=alpha+P_alpha*0.5d0*dtau_alpha
        !Fourier transform from xmat_mom to xmat
        call Fourier_transform_xmat(xmat,xmat_mom,myrank,p2x)
        !adjust the margin.
        call Adjust_margin_xmat(xmat,myrank)
        !Fourier transform from P_xmat_mom to P_xmat
        call Fourier_transform_P_xmat(P_xmat,P_xmat_mom,myrank,p2x)

    end if

    return

END subroutine Molecular_Dynamics
