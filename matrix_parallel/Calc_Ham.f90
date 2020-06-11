subroutine Calc_Ham(temperature,xmat,alpha,&
    &P_xmat,P_alpha,ham,myrank,pf,chi,acoeff_md,g_R,RCUT,nbmn,flux,ngauge,purebosonic,&
    &polfix,g_pol,pol_fix_width,myersfix,g_myers,myers_fix_width)

    implicit none
  include 'mpif.h'
  include 'size_parallel.h'
    !***** input *****
    integer nbmn,ngauge,purebosonic
    integer myrank
    double complex xmat(1:nmat_block,1:nmat_block,1:ndim,&
        &-(nmargin-1):nsite_local+nmargin)
    double complex P_xmat(1:nmat_block,1:nmat_block,1:ndim,1:nsite_local)
    double precision alpha(1:nmat_block*nblock)
    double precision P_alpha(1:nmat_block*nblock)
    double precision temperature,flux
    double precision g_R,RCUT
    double complex pf(1:nmat_block,1:nmat_block,1:nspin,&
        &-(nmargin-1):nsite_local+nmargin)
    double complex Chi(1:nremez_md,1:nmat_block,1:nmat_block,1:nspin,&
        &-(nmargin-1):nsite_local+nmargin)
    double precision acoeff_md(0:nremez_md)
    !***** output *****
    double precision ham
    !******************
    double precision alpha_max,alpha_min
    double precision pi
    double precision trx2
    integer isite
    integer idim
    integer imat,jmat,kmat
    integer iremez,ispin
    double complex x23(1:nmat_block,1:nmat_block),&
        x32(1:nmat_block,1:nmat_block),&
        trx123,trx132
    double precision trx2_123,trx2_456789,lattice_spacing

 

    double complex xmat_row(1:nmat_block,1:nmat_block*nblock,1:ndim,1:nsite_local)
    double complex xmat_column(1:nmat_block*nblock,1:nmat_block,&
        &1:ndim,1:nsite_local)

    double precision g_pol,pol_fix_width,pol,polfix
    double precision g_myers,myers_fix_width,myers,myersfix
    ! double precision polfix,g_pol,pol_fix_width,myersfix,g_myers,myers_fix_width

    pi=2d0*dasin(1d0)
    lattice_spacing=1d0/temperature/dble(nsite_local*nsublat)

    !move i-th row and j-th row of xmat to (i,j)-th node.
    call mpi_xmat_row(xmat,xmat_row,myrank)
    call mpi_xmat_column(xmat,xmat_column,myrank)

    call Calc_action(temperature,xmat,xmat_row,xmat_column,&
        &alpha,ham,myrank,ngauge)

    do isite=1,nsite_local
        do idim=1,ndim
            do imat=1,nmat_block
                do jmat=1,nmat_block
                    ham=ham&
                        +dble(P_xmat(imat,jmat,idim,isite)&
                        *dconjg(P_xmat(imat,jmat,idim,isite)))*0.5d0
                end do
            end do
        end do
    end do
    if(ngauge.eq.0)then
        if(myrank.EQ.0)then
            do imat=1,nmat_block*nblock
                ham=ham+dble(P_alpha(imat)*P_alpha(imat))*0.5d0
            end do
        end if
    end if
    !**********************
    !*** pseudo-fermion ***
    !**********************
    if(purebosonic.eq.0) then
        do iremez=1,nremez_md
            do isite=1,nsite_local
                do ispin=1,nspin
                    do imat=1,nmat_block
                        do jmat=1,nmat_block
                            ham=ham+acoeff_md(iremez)*&
                                dble(Chi(iremez,imat,jmat,ispin,isite)*&
                                dconjg(pf(imat,jmat,ispin,isite)))
                        end do
                    end do
                end do
            end do
        end do
    end if
    !****************************
    !*** constraint for alpha ***
    !****************************
 if(ngauge.eq.0)then
    if(myrank.EQ.0)then
        alpha_max=alpha(1)
        alpha_min=alpha(1)
        do imat=2,nmat_block*nblock
            if(alpha(imat).GT.alpha_max)then
                alpha_max=alpha(imat)
            else if(alpha(imat).LT.alpha_min)then
                alpha_min=alpha(imat)
            end if
        end do
        if(alpha_max-alpha_min.LT.2d0*pi)then
            ham=ham-dlog(2d0*pi-(alpha_max-alpha_min))
        end if
    end if
  end if
    !****************************
    !*** constraint for TrX^2 ***
    !****************************
    call Calc_TrX2(xmat,trx2,myrank)
    if(myrank.EQ.0)then
        if(trx2.GE.RCUT)then
            ham=ham+g_R*(trx2-RCUT)*dble(nmat_block*nblock)
        end if
    end if

  !***************************
  !*** fixing polakov loop ***
  !***************************
  if((ngauge.EQ.0).AND.(g_pol.NE.0))then
     call Calc_Polyakov(nmat_block*nblock,alpha,Pol)
     if(myrank.EQ.0)then
        if(pol.GT.(polfix+0.5d0*pol_fix_width))then
           ham=ham+0.5d0*g_pol*(pol-(polfix+0.5d0*pol_fix_width))**2d0
        else if(pol.LT.(polfix-0.5d0*pol_fix_width))then
           ham=ham+0.5d0*g_pol*(pol-(polfix-0.5d0*pol_fix_width))**2d0
        end if
     end if
  end if
    !****************************
    !*** constraint for Myers observable ***
    !****************************
  if(g_myers.NE.0)then
     call Calc_Myers(xmat,myers,myrank)
     if(myrank.EQ.0)then
        if(myers.GT.(myersfix+0.5d0*myers_fix_width))then
           ham=ham+0.5d0*g_myers*(myers-(myersfix+0.5d0*myers_fix_width))**2d0
        else if(myers.LT.(myersfix-0.5d0*myers_fix_width))then
           ham=ham+0.5d0*g_myers*(myers-(myersfix-0.5d0*myers_fix_width))**2d0
        end if
     end if
  end if
    !******************************
    !*** Plane wave deformation ***
    !******************************
    if(nbmn.EQ.1)then
        !*****************
        !*** mass term ***
        !*****************
        trx2_123=0d0
        trx2_456789=0d0
        do isite=1,nsite_local
            do imat=1,nmat_block
                do jmat=1,nmat_block
                    do idim=1,3
                        trx2_123=trx2_123&
                            +dble(xmat(imat,jmat,idim,isite)&
                            *dconjg(xmat(imat,jmat,idim,isite)))
                    end do
                    do idim=4,9
                        trx2_456789=trx2_456789&
                            +dble(xmat(imat,jmat,idim,isite)&
                            *dconjg(xmat(imat,jmat,idim,isite)))
                    end do
                end do
            end do
        end do
        ham=ham+flux*flux*(0.5d0*trx2_123+0.125d0*trx2_456789)&
            *lattice_spacing*dble(nmat_block*nblock)
        !******************
        !*** cubic term ***
        !******************
        trx123=(0d0,0d0)
        trx132=(0d0,0d0)
        do isite=1,nsite_local
            x23=(0d0,0d0)
            x32=(0d0,0d0)
            !$omp parallel
            !$omp do
            do imat=1,nmat_block
                do jmat=1,nmat_block
                    do kmat=1,nmat_block*nblock
                        x23(imat,jmat)=x23(imat,jmat)&
                            +xmat_row(imat,kmat,2,isite)&
                            *xmat_column(kmat,jmat,3,isite)
                        x32(imat,jmat)=x32(imat,jmat)&
                            +xmat_row(imat,kmat,3,isite)&
                            *xmat_column(kmat,jmat,2,isite)
                    end do
                end do
            end do
            !$omp end do
            !$omp end parallel
            do imat=1,nmat_block
                do jmat=1,nmat_block
                    trx123=trx123+x23(imat,jmat)*dconjg(xmat(imat,jmat,1,isite))
                    trx132=trx132+x32(imat,jmat)*dconjg(xmat(imat,jmat,1,isite))
                end do
            end do
        end do
        ham=ham+dble((0d0,3d0)*(trx123-trx132))*flux*lattice_spacing*dble(nmat_block*nblock)
    end if

    return

END subroutine Calc_Ham
