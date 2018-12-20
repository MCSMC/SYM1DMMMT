SUBROUTINE RHMC_evolution(xmat,alpha,ncv,n_bad_CG,nacceptance,nbc,nbmn,&
    &temperature,flux,GAMMA10d,ntau,dtau_xmat,dtau_alpha,&
    &acceleration,g_alpha,g_R,RCUT,&
    &acoeff_md,bcoeff_md,acoeff_pf,bcoeff_pf,&
    &max_err,max_iteration,iteration,&
    &ham_init,ham_fin,ntrial,imetropolis,nsmear,s_smear,ngauge,purebosonic)

    use mtmod !Mersenne twistor
    implicit none

  include 'mpif.h'
  include 'size_parallel.h'
  include '../Fourier.inc'
  include '../unit_number.inc'
    !input
    integer nbc,nbmn,nsmear,ngauge,purebosonic
    double precision temperature,flux,s_smear
    double complex GAMMA10d(1:ndim,1:nspin,1:nspin)
    double precision max_err
    integer max_iteration
    double precision acoeff_md(0:nremez_md),bcoeff_md(1:nremez_md)
    double precision acoeff_pf(0:nremez_pf),bcoeff_pf(1:nremez_pf)
    double precision g_alpha,g_R,RCUT
    integer ntau
    integer imetropolis!1-> no Metropolis test
    double precision acceleration(1:nsite_local),dtau_alpha,dtau_xmat
    !input & output
    double complex xmat(1:nmat_block,1:nmat_block,1:ndim,-(nmargin-1):nsite_local+nmargin)
    double precision alpha(1:nmat_block*nblock)
    integer ncv,n_bad_CG,nacceptance,ntrial
    double precision ham_init,ham_fin
    !output
    integer iteration


    double precision metropolis
    integer info_pf,info_mol,info_CG_init,info_CG_fin,info_alpha,info_accept,info
    double complex backup_xmat(1:nmat_block,1:nmat_block,1:ndim,-(nmargin-1):nsite_local+nmargin)
    double precision backup_alpha(1:nmat_block*nblock)
    double complex P_xmat(1:nmat_block,1:nmat_block,1:ndim,1:nsite_local)
    double precision P_alpha(1:nmat_block*nblock)
    double complex pf(1:nmat_block,1:nmat_block,&
        1:nspin,-(nmargin-1):nsite_local+nmargin)
    double complex Chi(1:nremez_md,1:nmat_block,1:nmat_block,&
        1:nspin,-(nmargin-1):nsite_local+nmargin)
    double precision ham_init_local,ham_fin_local
    integer IERR,myrank
    double complex xmat_smeared(1:nmat_block,1:nmat_block,1:ndim,-(nmargin-1):nsite_local+nmargin)
  
    call MPI_COMM_RANK(MPI_COMM_WORLD,MYRANK, IERR)



    !initialize info flugs
    info_pf=0
    info_mol=0
    info_CG_init=0
    info_CG_fin=0
    iteration=0
    !**********************************
    !**** Generate pseudo fermion. ****
    !**********************************
    !smear
    if(purebosonic.eq.0) then
        call smearing_xmat(xmat,xmat_smeared,s_smear,myrank,nsmear)
  

        call generate_pseudo_fermion_SUN(nbc,acoeff_pf,bcoeff_pf,temperature,&
            &xmat_smeared,alpha,pf,GAMMA10d,max_err,max_iteration,iteration,myrank,&
            &nbmn,flux,info_pf)
    end if
    !info_pf=0 -> OK (CG solver converged)
    !info_pf=1 -> error (CG solver did not converge)

    if(myrank.EQ.0)then
        !Take CG_log
        write(unit_CG_log,*)"pseudo-fermion generation",iteration
    end if
    !******************************
    !**** Molecular Evolution. ****
    !******************************
    !save the old config.
    backup_xmat=xmat
    backup_alpha=alpha
    !generate auxiliary momenta.
    call Generate_P_xmat(P_xmat,myrank)
    if(ngauge.EQ.0)then
        !gauged
        call Generate_P_alpha(nmat_block*nblock,P_alpha)
    else if(ngauge.EQ.1)then
        !ungauged
        P_alpha=0d0
    end if
    !Calculate ham_init
    if(purebosonic.eq.0) then
        call solver_biCGm(nbc,nremez_md,bcoeff_md,temperature,&
            &xmat_smeared,alpha,pf,chi,GAMMA10d,max_err,max_iteration,iteration,myrank,&
            &nbmn,flux,info_CG_init)
    end if
    !info_CG_init=0 -> OK (CG solver converged)
    !info_CG_init=1 -> error (CG solver did not converge)

    if(myrank.EQ.0)then
        !Take CG_log
        write(unit_CG_log,*)"ham_init",iteration
    end if

    call Calc_Ham(temperature,xmat,alpha,&
        &P_xmat,P_alpha,ham_init_local,myrank,pf,chi,&
        &acoeff_md,g_R,RCUT,nbmn,flux,ngauge,purebosonic)
    !collect ham_int to myrank=0 and calculate total value
    call MPI_Reduce(ham_init_local,ham_init,1,MPI_DOUBLE_PRECISION,&
        &MPI_SUM,0,MPI_COMM_WORLD,IERR)
    !Molecular Evolution
    call Molecular_Dynamics(nbc,temperature,&
        &ntau,dtau_xmat,dtau_alpha,xmat,alpha,P_xmat,P_alpha,&
        &acoeff_md,bcoeff_md,pf,max_iteration,max_err,iteration,&
        &gamma10d,g_alpha,g_R,RCUT,acceleration,&
        &nbmn,flux,info_mol,nsmear,s_smear,ngauge,purebosonic)
    !info_mol=0 -> OK (CG solver converged)
    !info_mol=1 -> error (CG solver did not converge)
    if(info_mol.EQ.0)then
        if(purebosonic.eq.0) then
            !smear
            call smearing_xmat(xmat,xmat_smeared,s_smear,myrank,nsmear)
            !calculate ham_fin
            call solver_biCGm(nbc,nremez_md,bcoeff_md,temperature,&
                &xmat_smeared,alpha,pf,chi,GAMMA10d,max_err,max_iteration,iteration,myrank,&
                &nbmn,flux,info_CG_fin)
        !info_CG_fin=0 -> OK (CG solver converged)
        !info_CG_fin=1 -> error (CG solver did not converge)
        end if
        if(myrank.EQ.0)then
            !Take CG_log
            write(unit_CG_log,*)"ham_fin",iteration
        end if
     
        call Calc_Ham(temperature,&
            &xmat,alpha,P_xmat,P_alpha,ham_fin_local,myrank,&
            &pf,chi,acoeff_md,g_R,RCUT,nbmn,flux,ngauge,purebosonic)
        !collect ham_fin to myrank=0 and calculate total value
        call MPI_Reduce(ham_fin_local,ham_fin,1,MPI_DOUBLE_PRECISION,&
            &MPI_SUM,0,MPI_COMM_WORLD,IERR)
    end if
    !#######################################
    !#### Did CG solver work correctly? ####
    !#######################################
    if((info_pf.EQ.0).AND.(info_mol.EQ.0).AND.(info_CG_init.EQ.0).AND.(info_CG_fin.EQ.0))then
        info=0!CG solver worked fine everywhere.
    else
        info=1!CG solver did not converge at least one place.
    end if
    !Count how many times the CG solver failed to converge.
    if(info.EQ.1)then
        n_bad_CG=n_bad_CG+1
    end if
    !#########################
    !#### metropolis test ####
    !#########################
    ntrial=ntrial+1
    metropolis=grnd()
    !"metropolis=grnd()" must be put outside "if(myrank.EQ.0)then"
    !so that random numbers synchronize.
    if(imetropolis.EQ.1)then
        !No Metropolis test
        metropolis=-1d0 !Then, it is always accepted.
    end if
    if(myrank.EQ.0)then

        !####################################
        !#### check constraint for alpha ####
        !####################################
        call check_alpha_constraint(nmat_block*nblock,alpha,info_alpha)
        !Count how many times the constraint is violated.
        if(info_alpha.EQ.1)then
            ncv=ncv+1
            info=1
        end if
        if(info.EQ.0)then
            !alpha satisfies the constraint.
            !CG converged everywhere during Molecular evolution.
        
            if(dexp(ham_init-ham_fin) > metropolis)THEN
                !accept
                info_accept=0
            else
                !reject
                info_accept=1
            end if
        else!automatic reject, before the Metropolis test.
            info_accept=1
        end if
    end if
    !"info_accept" is broadcasted to all nodes.
    !if(ngauge.EQ.0)then
    call MPI_Bcast(info_accept,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
    !end if
    ! accept/reject at each node.
    if(info_accept.EQ.0)then
        !accept
        nacceptance=nacceptance+1
        ! Project out the trace parts of X and alpha.
        ! For X, we remove only the zero mode of X,
        !        i.e. (¥int dt X) is made traceless.


        call subtract_U1(xmat,alpha,myrank)

    else
        !reject
        xmat=backup_xmat
        alpha=backup_alpha
    end if
  
    return

END SUBROUTINE RHMC_evolution
