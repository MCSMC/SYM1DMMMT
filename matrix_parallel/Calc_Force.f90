! alpha -> alpha + P_alpha*dtau_alpha
! xmat -> xmat + P_xmat*dtau_xmat
! P_alpha -> P_alpha - delh_alpha*dtau_alpha
! P_xmat -> P_xmat - delh_xmat*dtau_xmat
! delh_xmat(imat,jmat)=dS/dxmat(jmat,imat), be careful about the difference of 
! ordering of indices. 
SUBROUTINE Calc_Force(delh_xmat,delh_alpha,xmat,alpha,&
    chi,GAMMA10d,gcoeff_alpha,g_R,RCUT,nbmn,temperature,flux,acoeff_md,&
    &xmat_smeared,nsmear,s_smear,ngauge,purebosonic,&
    &polfix,g_pol,pol_fix_width,myersfix,g_myers,myers_fix_width)

    implicit none

  include 'mpif.h'
  include 'size_parallel.h'
    !***** input *****
    integer nbmn,myrank,nsmear,ngauge,purebosonic
    double precision temperature,flux,s_smear
    double precision gcoeff_alpha
    double precision g_R,trx2,RCUT
    double complex xmat(1:nmat_block,1:nmat_block,1:ndim,&
        &-(nmargin-1):nsite_local+nmargin)
    double complex xmat_smeared(1:nmat_block,1:nmat_block,1:ndim,&
        &-(nmargin-1):nsite_local+nmargin)
    double precision alpha(1:nmat_block*nblock),alpha_max,alpha_min
    double complex Chi(1:nremez_md,1:nmat_block,1:nmat_block,&
        1:nspin,-(nmargin-1):nsite_local+nmargin)
    double complex GAMMA10d(1:ndim,1:nspin,1:nspin)
    double precision acoeff_md(0:nremez_md)
    !***** output *****
    double complex delh_xmat(1:nmat_block,1:nmat_block,1:ndim,1:nsite_local)
    double precision delh_alpha(1:nmat_block*nblock)
    !******************
    double precision lattice_spacing
    integer imat,jmat,kmat,lmat,imax,imin
    integer idim,jdim
    integer isite
    double precision pi
    double complex commutator(1:nmat_block,1:nmat_block,1:ndim,1:ndim),&
        &uxdu(1:nmat_block,1:nmat_block),&
        &udxu(1:nmat_block,1:nmat_block),&
        &ei,ej
    double complex com_x(1:nmat_block,1:nmat_block,1:nsite_local,1:ndim,1:nblock)
    double complex x_com(1:nmat_block,1:nmat_block,1:nsite_local,1:ndim,1:nblock)
    double complex Deriv_Mchi_alpha(1:nremez_md,1:nmat_block*nblock,1:nmat_block,&
        &1:nmat_block,1:nspin,1:nsite_local)
    double complex MChi(1:nremez_md,1:nmat_block,1:nmat_block,&
        1:nspin,1:nsite_local)
    double complex Chi_column(1:nremez_md,1:nmat_block,1:nmat_block,&
        &1:nspin,1:nsite_local,1:nblock)
    double complex MChi_column(1:nremez_md,1:nmat_block,1:nmat_block,&
        &1:nspin,1:nsite_local,1:nblock)
    double complex Chi_row(1:nremez_md,1:nmat_block,1:nmat_block,&
        &1:nspin,1:nsite_local,1:nblock)
    double complex MChi_row(1:nremez_md,1:nmat_block,1:nmat_block,&
        &1:nspin,1:nsite_local,1:nblock)
    double complex delPF_xmat(1:nmat_block,1:nmat_block,1:ndim,1:nsite_local)
    double complex delPF_xmat_smeared(1:nmat_block,1:nmat_block,1:ndim,1:nsite_local)
    double precision delPF_alpha(1:nmat_block*nblock)
    integer iremez
    integer ispin,jspin
    double complex com12(1:nmat_block,1:nmat_block),&
        &com23(1:nmat_block,1:nmat_block),&
        &com31(1:nmat_block,1:nmat_block),temp
    double complex xmat_row(1:nmat_block,1:nmat_block*nblock,1:ndim,1:nsite_local)
    double complex xmat_column(1:nmat_block*nblock,1:nmat_block,1:ndim,1:nsite_local)
    double complex xmat_row_2(1:nmat_block,1:nmat_block*nblock,1:ndim,1:nsite_local)
    double complex xmat_column_2(1:nmat_block*nblock,1:nmat_block,1:ndim,1:nsite_local)
    double complex xmat_row_smeared(1:nmat_block,1:nmat_block*nblock,1:ndim,1:nsite_local)
    double complex xmat_column_smeared(1:nmat_block*nblock,1:nmat_block,1:ndim,1:nsite_local)
    double complex uxumx(1:nmat_block,1:nmat_block)
    !***** for MPI *****
    integer send_rank,receive_rank,ireq,ierr,tag
    integer status(MPI_STATUS_SIZE)
    double complex mat_rcv(1:nmat_block,1:nmat_block,1:nsite_local,1:ndim)
    double complex mat_send(1:nmat_block,1:nmat_block,1:nsite_local,1:ndim)
    double complex mat_rcv_2(1:nmat_block,1:nmat_block,1:ndim,1:nsite_local)
    integer iblock,jblock,kblock,isublat,ishift,kblock_send,kblock_rcv
    double precision polfix,g_pol,pol_fix_width,myersfix,g_myers,myers_fix_width
    double precision pol,re_pol,im_pol,myers

    call MPI_COMM_RANK(MPI_COMM_WORLD,MYRANK, IERR)
    call who_am_i(myrank,isublat,iblock,jblock)
    !move i-th row and j-th row of xmat to (i,j)-th node.
    call mpi_xmat_row(xmat,xmat_row,myrank)
    call mpi_xmat_column(xmat,xmat_column,myrank)
    call mpi_xmat_row(xmat_smeared,xmat_row_smeared,myrank)
    call mpi_xmat_column(xmat_smeared,xmat_column_smeared,myrank)

  
    do idim=1,ndim
        !$omp parallel
        !$omp do
        do imat=1,nmat_block
            do jmat=1,nmat_block*nblock
                do isite=1,nsite_local
                    xmat_row_2(imat,jmat,idim,isite)=&
                        &dconjg(xmat_column(jmat,imat,idim,isite))
                    xmat_column_2(jmat,imat,idim,isite)=&
                        &dconjg(xmat_row(imat,jmat,idim,isite))
                end do
            end do
        end do
    !$omp end do
    !$omp end parallel
    end do

    pi=2d0*dasin(1d0)
    lattice_spacing=1d0/temperature/dble(nsite_local*nsublat)
    !*****************************
    !*****************************
    !*** calculate delh_alpha ****
    !*****************************
    !*****************************
    delh_alpha=0d0
    !********************
    !*** kinetic term ***
    !********************
    if(nimprove.EQ.0)then
        !naive action.
        do isite=1,nsite_local
            !Neighboring MPI processes take care of margins.
            do idim=1,ndim
                !u(t)*x(t+a)*u^dagger(t) - x(t)
                do imat=1,nmat_block
                    do jmat=1,nmat_block
                        !exp(i*alpha_i)
                        ei=dcmplx(dcos(alpha(imat+(iblock-1)*nmat_block)&
                            &/dble(nsite_local*nsublat)))&
                            &+(0d0,1d0)*dcmplx(dsin(alpha(imat+(iblock-1)*nmat_block)&
                            &/dble(nsite_local*nsublat)))
                        !exp(-i*alpha_j)
                        ej=dcmplx(dcos(alpha(jmat+(jblock-1)*nmat_block)&
                            &/dble(nsite_local*nsublat)))&
                            -(0d0,1d0)*dcmplx(dsin(alpha(jmat+(jblock-1)*nmat_block)&
                            &/dble(nsite_local*nsublat)))
                 
                        delh_alpha(imat+(iblock-1)*nmat_block)&
                            &=delh_alpha(imat+(iblock-1)*nmat_block)&
                            &-dble(ei*xmat(imat,jmat,idim,isite+1)&
                            &*ej*dconjg(xmat(imat,jmat,idim,isite))*(0d0,1d0))
                 
                        delh_alpha(jmat+(jblock-1)*nmat_block)&
                            &=delh_alpha(jmat+(jblock-1)*nmat_block)&
                            &+dble(ei*xmat(imat,jmat,idim,isite+1)&
                            &*ej*dconjg(xmat(imat,jmat,idim,isite))*(0d0,1d0))
                 
                 
                    end do
                end do
            end do
        end do
    else if(nimprove.EQ.1)then
        !improved action.
        do isite=1,nsite_local
            !Neighboring MPI processes take care of margins.
            do idim=1,ndim
                !-0.5*u^2*x(t+2a)*(u^dagger)^2 + 2*u*x(t+a)*u^dagger - 1.5*x(t)
                do imat=1,nmat_block
                    do jmat=1,nmat_block
                        !exp(i*alpha_i)
                        ei=dcmplx(dcos(alpha(imat+(iblock-1)*nmat_block)&
                            &/dble(nsite_local*nsublat)))&
                            &+(0d0,1d0)*dcmplx(dsin(alpha(imat+(iblock-1)*nmat_block)&
                            &/dble(nsite_local*nsublat)))
                        !exp(-i*alpha_j)
                        ej=dcmplx(dcos(alpha(jmat+(jblock-1)*nmat_block)&
                            &/dble(nsite_local*nsublat)))&
                            -(0d0,1d0)*dcmplx(dsin(alpha(jmat+(jblock-1)*nmat_block)&
                            &/dble(nsite_local*nsublat)))
                        uxumx(imat,jmat)=&
                            -(0.5d0,0d0)*ei*ei*xmat(imat,jmat,idim,isite+2)*ej*ej&
                            +(2d0,0d0)*ei*xmat(imat,jmat,idim,isite+1)*ej&
                            -(1.5d0,0d0)*xmat(imat,jmat,idim,isite)
                 
                        delh_alpha(imat+(iblock-1)*nmat_block)&
                            &=delh_alpha(imat+(iblock-1)*nmat_block)&
                            &+dble(((2d0,0d0)*ei*xmat(imat,jmat,idim,isite+1)*ej&
                            &-ei*ei*xmat(imat,jmat,idim,isite+2)*ej*ej)&
                            *dconjg(uxumx(imat,jmat))*(0d0,1d0))
                 
                        delh_alpha(jmat+(jblock-1)*nmat_block)&
                            &=delh_alpha(jmat+(jblock-1)*nmat_block)&
                            &-dble(((2d0,0d0)*ei*xmat(imat,jmat,idim,isite+1)*ej&
                            &-ei*ei*xmat(imat,jmat,idim,isite+2)*ej*ej)&
                            *dconjg(uxumx(imat,jmat))*(0d0,1d0))

                    end do
                end do
            end do
        end do
    end if
    delh_alpha=delh_alpha*dble(nmat_block*nblock)&
        &/lattice_spacing/dble(nsite_local*nsublat)
    !*************************
    !*** gauge-fixing term ***
    !*************************
    if(ngauge.EQ.0)then
        if(myrank.EQ.0)then
            !$omp parallel
            !$omp do
            do imat=1,nmat_block*nblock
                do jmat=1,nmat_block*nblock
                    if(jmat.NE.imat)then
                        delh_alpha(imat)=delh_alpha(imat)&
                            -1d0/dtan(0.5d0*(alpha(imat)-alpha(jmat)))
                    end if
                end do
            end do
           !$omp end do
           !$omp end parallel
        end if
    end if
    !***************************
    !***************************
    !*** calculate delh_xmat ***
    !***************************
    !***************************
    delh_xmat=(0d0,0d0)
    !********************
    !*** kinetic term ***
    !********************
    if(nimprove.EQ.0)then
        !naive action.
        do isite=1,nsite_local
            !Neighboring MPI processes take care of margins.
            do idim=1,ndim
                !u(t)*x(t+a)*u(t)^¥dagger
                !u(t-a)^¥dagger*x(t-a)*u(t-a)
                !$omp parallel
                !$omp do
                do imat=1,nmat_block
                    do jmat=1,nmat_block
                 
                        !exp(i*alpha_i)
                        ei=dcmplx(dcos(alpha(imat+(iblock-1)*nmat_block)&
                            &/dble(nsite_local*nsublat)))&
                            &+(0d0,1d0)*dcmplx(dsin(alpha(imat+(iblock-1)*nmat_block)&
                            &/dble(nsite_local*nsublat)))
                        !exp(-i*alpha_j)
                        ej=dcmplx(dcos(alpha(jmat+(jblock-1)*nmat_block)&
                            &/dble(nsite_local*nsublat)))&
                            &-(0d0,1d0)*dcmplx(dsin(alpha(jmat+(jblock-1)*nmat_block)&
                            &/dble(nsite_local*nsublat)))
                 
                        uxdu(imat,jmat)=&
                            ei*xmat(imat,jmat,idim,isite+1)*ej
                        udxu(imat,jmat)=&
                            dconjg(ei)*xmat(imat,jmat,idim,isite-1)*dconjg(ej)
                    end do
                end do
                !$omp end do
                !$omp end parallel
                !$omp parallel
                !$omp do
                do imat=1,nmat_block
                    do jmat=1,nmat_block
                        delh_xmat(imat,jmat,idim,isite)=&
                            &delh_xmat(imat,jmat,idim,isite)&
                            &+dcmplx(2d0)*xmat(imat,jmat,idim,isite)
                 
                        delh_xmat(imat,jmat,idim,isite)=&
                            &delh_xmat(imat,jmat,idim,isite)&
                            &-uxdu(imat,jmat)&
                            &-udxu(imat,jmat)
                    end do
                end do
            !$omp end do
            !$omp end parallel
            end do
        end do
    else if(nimprove.EQ.1)then
        !improved action.
        do isite=1,nsite_local
            do idim=1,ndim
                !$omp parallel
                !$omp do
                do imat=1,nmat_block
                    do jmat=1,nmat_block
                        !exp(i*alpha_i)
                        ei=dcmplx(dcos(alpha(imat+(iblock-1)*nmat_block)&
                            &/dble(nsite_local*nsublat)))&
                            &+(0d0,1d0)*dcmplx(dsin(alpha(imat+(iblock-1)*nmat_block)&
                            &/dble(nsite_local*nsublat)))
                        !exp(-i*alpha_j)
                        ej=dcmplx(dcos(alpha(jmat+(jblock-1)*nmat_block)&
                            &/dble(nsite_local*nsublat)))&
                            &-(0d0,1d0)*dcmplx(dsin(alpha(jmat+(jblock-1)*nmat_block)&
                            &/dble(nsite_local*nsublat)))
                 
                        delh_xmat(imat,jmat,idim,isite)=&
                            &delh_xmat(imat,jmat,idim,isite)&
                            &+(6.5d0,0d0)*xmat(imat,jmat,idim,isite)&
                            &-(4d0,0d0)*ei*xmat(imat,jmat,idim,isite+1)*ej&
                            &-(4d0,0d0)*dconjg(ei)*xmat(imat,jmat,idim,isite-1)&
                            &*dconjg(ej)&
                            &+(0.75d0,0d0)*ei*ei*xmat(imat,jmat,idim,isite+2)*ej*ej&
                            &+(0.75d0,0d0)*dconjg(ei*ei)&
                            &*xmat(imat,jmat,idim,isite-2)*dconjg(ej*ej)
                 
                    end do
                end do
            !$omp end do
            !$omp end parallel
            end do
        end do
    end if
    delh_xmat=delh_xmat*dcmplx(nmat_block*nblock)/dcmplx(lattice_spacing)
    !***********************
    !*** commutator term ***
    !***********************
    com_x=(0d0,0d0)
    x_com=(0d0,0d0)
    !these are parts of com*X and x*com
    do isite=1,nsite_local
        commutator=(0d0,0d0)
        !Information of the (iblock,jblock)-component of the commutator is
        !stored at (iblock,jblock).
        do idim=1,ndim
            do jdim=1,ndim
                !$omp parallel
                !$omp do
                do imat=1,nmat_block
                    do jmat=1,nmat_block
                        do kmat=1,nmat_block*nblock
                            commutator(imat,jmat,idim,jdim)=&
                                commutator(imat,jmat,idim,jdim)&
                                &+xmat_row(imat,kmat,idim,isite)&
                                &*xmat_column(kmat,jmat,jdim,isite)&
                                &-xmat_row(imat,kmat,jdim,isite)&
                                &*xmat_column(kmat,jmat,idim,isite)
                        end do
                    end do
                end do
            !$omp end do
            !$omp end parallel
            end do
        end do
        !com_x(1:nmat_block,1:nmat_block,1:nsite_local,1:ndim,1:nblock)
        !x_com(1:nmat_block,1:nmat_block,1:nsite_local,1:ndim,1:nblock)
        do idim=1,ndim
            do jdim=1,ndim
                !$omp parallel
                !$omp do
                do imat=1,nmat_block
                    do jmat=1,nmat_block
                        do kblock=1,nblock
                            do kmat=1,nmat_block
                                !Information of the (kblock,jblock)-component of x_com is
                                !stored at (iblock,jblock). [kblock=1,2,...,nblock]
                                x_com(imat,jmat,isite,idim,kblock)=&
                                    &x_com(imat,jmat,isite,idim,kblock)&
                                    &+xmat_column_2(imat+(kblock-1)*nmat_block,&
                                    &kmat,jdim,isite)&
                                    &*commutator(kmat,jmat,idim,jdim)
                                !Information of the (iblock,kblock)-component of com_x is
                                !stored at (iblock,jblock). [kblock=1,2,...,nblock]
                                com_x(imat,jmat,isite,idim,kblock)=&
                                    &com_x(imat,jmat,isite,idim,kblock)&
                                    &+commutator(imat,kmat,idim,jdim)&
                                    &*xmat_row_2(kmat,jmat+(kblock-1)*nmat_block,&
                                    &jdim,isite)
                            end do
                        end do
                    end do
                end do
            !$omp end do
            !$omp end parallel
            end do
        end do
    end do
    !###############################################
    !## xxx=xxx+&                                 ##
    !## xmat(imat,kmat,jdim,isite)&               ##
    !## *commutator(kmat,jmat,idim,jdim)&         ##
    !## -commutator(imat,kmat,idim,jdim)&         ##
    !## *xmat(kmat,jmat,jdim,isite)               ##
    !##                                           ##
    !## delh_xmat(imat,jmat,idim,isite)=&         ##
    !## delh_xmat(imat,jmat,idim,isite)&          ##
    !## -dcmplx(nmat)*dcmplx(lattice_spacing)*xxx ##
    !###############################################
    !$omp parallel
    !$omp do
    do imat=1,nmat_block
        do jmat=1,nmat_block
            do idim=1,ndim
                do isite=1,nsite_local
                    delh_xmat(imat,jmat,idim,isite)=&
                        &delh_xmat(imat,jmat,idim,isite)&
                        &-dcmplx(nmat_block*nblock)*dcmplx(lattice_spacing)&
                        &*x_com(imat,jmat,isite,idim,iblock)
                end do
            end do
        end do
    end do
    !$omp end do
    !$omp end parallel
    !Information of the (kblock,jblock)-component of x_com is
    !stored at (iblock,jblock). [kblock=1,2,...,nblock]
    !Hence, we should send them to (kblock,jblock).
    kblock_send=iblock
    kblock_rcv=iblock
    do ishift=1,nblock-1
        kblock_send=kblock_send+1
        kblock_rcv=kblock_rcv-1
        if(kblock_send.EQ.nblock+1)then
            kblock_send=1
        end if
        if(kblock_rcv.EQ.0)then
            kblock_rcv=nblock
        end if
        send_rank=(isublat-1)*nblock*nblock+(kblock_send-1)*nblock+jblock-1
        receive_rank=(isublat-1)*nblock*nblock+(kblock_rcv-1)*nblock+jblock-1
        !double complex mat_rcv(1:nmat_block,1:nmat_block,1:nsite_local,1:ndim)
        !double complex mat_send(1:nmat_block,1:nmat_block,1:nsite_local,1:ndim)
        !$omp parallel
        !$omp do
        do imat=1,nmat_block
            do jmat=1,nmat_block
                do idim=1,ndim
                    do isite=1,nsite_local
                        mat_send(imat,jmat,isite,idim)=&
                            &x_com(imat,jmat,isite,idim,kblock_send)
                    end do
                end do
            end do
        end do
        !$omp end do
        !$omp end parallel
        tag=1
        call MPI_Isend(mat_send(1,1,1,1),&
            &nmat_block*nmat_block*nsite_local*ndim,&
            &MPI_DOUBLE_COMPLEX,&
            &send_rank,tag,MPI_COMM_WORLD,ireq,ierr)
        call MPI_Recv(mat_rcv(1,1,1,1),&
            &nmat_block*nmat_block*nsite_local*ndim,&
            &MPI_DOUBLE_COMPLEX,&
            &receive_rank,tag,MPI_COMM_WORLD,status,ierr)
        call MPI_Wait(ireq,status,ierr)
        !$omp parallel
        !$omp do
        do imat=1,nmat_block
            do jmat=1,nmat_block
                do idim=1,ndim
                    do isite=1,nsite_local
                        delh_xmat(imat,jmat,idim,isite)=&
                            &delh_xmat(imat,jmat,idim,isite)&
                            &-dcmplx(nmat_block*nblock)*dcmplx(lattice_spacing)&
                            &*mat_rcv(imat,jmat,isite,idim)
                    end do
                end do
            end do
        end do
    !$omp end do
    !$omp end parallel
    end do
    !$omp parallel
    !$omp do
    do imat=1,nmat_block
        do jmat=1,nmat_block
            do idim=1,ndim
                do isite=1,nsite_local
                    delh_xmat(imat,jmat,idim,isite)=&
                        &delh_xmat(imat,jmat,idim,isite)&
                        &+dcmplx(nmat_block*nblock)*dcmplx(lattice_spacing)&
                        &*com_x(imat,jmat,isite,idim,jblock)
                end do
            end do
        end do
    end do
    !$omp end do
    !$omp end parallel
    !Information of the (iblock,kblock)-component of com_x is
    !stored at (iblock,jblock). [kblock=1,2,...,nblock]
    !Hence, we should send them to (iblock,kblock).
    kblock_send=jblock
    kblock_rcv=jblock
    do ishift=1,nblock-1
        kblock_send=kblock_send+1
        kblock_rcv=kblock_rcv-1
        if(kblock_send.EQ.nblock+1)then
            kblock_send=1
        end if
        if(kblock_rcv.EQ.0)then
            kblock_rcv=nblock
        end if
        send_rank=(isublat-1)*nblock*nblock+(iblock-1)*nblock+kblock_send-1
        receive_rank=(isublat-1)*nblock*nblock+(iblock-1)*nblock+kblock_rcv-1
        !$omp parallel
        !$omp do
        do imat=1,nmat_block
            do jmat=1,nmat_block
                do idim=1,ndim
                    do isite=1,nsite_local
                        mat_send(imat,jmat,isite,idim)=&
                            &com_x(imat,jmat,isite,idim,kblock_send)
                    end do
                end do
            end do
        end do
        !$omp end do
        !$omp end parallel
        tag=2
        call MPI_Isend(mat_send(1,1,1,1),&
            nmat_block*nmat_block*nsite_local*ndim,&
            &MPI_DOUBLE_COMPLEX,&
            &send_rank,tag,MPI_COMM_WORLD,ireq,ierr)
        call MPI_Recv(mat_rcv(1,1,1,1),nmat_block*nmat_block*nsite_local*ndim,&
            &MPI_DOUBLE_COMPLEX,&
            &receive_rank,tag,MPI_COMM_WORLD,status,ierr)
        call MPI_Wait(ireq,status,ierr)
        !$omp parallel
        !$omp do
        do imat=1,nmat_block
            do jmat=1,nmat_block
                do idim=1,ndim
                    do isite=1,nsite_local
                        delh_xmat(imat,jmat,idim,isite)=&
                            &delh_xmat(imat,jmat,idim,isite)&
                            &+dcmplx(nmat_block*nblock)*dcmplx(lattice_spacing)&
                            &*mat_rcv(imat,jmat,isite,idim)
                    end do
                end do
            end do
        end do
    !$omp end do
    !$omp end parallel
    end do
    !*******************************
    !*******************************
    !***** pseudo-fermion part *****
    !*******************************
    !*******************************
    !** xmat has to be replaced with the smeared one. ***
    !###########################################################################
    !## xmat-derivative of the pseudo-fermion action can easily be calculated ##
    !## from Mchi and chi.                                                    ##
    !## So, firstly we gather                                                 ##
    !## Mchi(:,jblock), Mchi(iblock,:), Chi(:,iblock),Chi(jblock,:)           ##
    !## to (iblock,jblock)-process.                                           ##
    !##                                                                       ##
    !## In (iblock,jblock)-process, we use the following notation:            ##
    !## Mchi(kblock,jblock) = Mchi_column(kblock)                             ##
    !## Mchi(iblock,kblock) = Mchi_row(kblock)                                ##
    !## chi(kblock,iblock) = chi_column(kblock)                               ##
    !## chi(jblock,kblock) = chi_row(kblock)                                  ##
    !## Note that the ordering of the indices of chi are tricky.              ##
    !## It makes the MPI communication a bit tricky too.                      ##
    !##                                                                       ##
    !## The format is                                                         ##
    !##                                                                       ##
    !## Chi_column(1:nremez_md,1:nmat_block,1:nmat_block,&                    ##
    !##           &1:nspin,0:nsite_local+1,1:nblock)                          ##
    !## MChi_column(1:nremez_md,1:nmat_block,1:nmat_block,&                   ##
    !##           &1:nspin,1:nsite_local,1:nblock)                            ##
    !## Chi_row(1:nremez_md,1:nmat_block,1:nmat_block,&                       ##
    !##           &1:nspin,0:nsite_local+1,1:nblock)                          ##
    !## MChi_row(1:nremez_md,1:nmat_block,1:nmat_block,&                      ##
    !##           &1:nspin,1:nsite_local,1:nblock)                            ##
    !###########################################################################
    !calculate Mchi.
    if(purebosonic.eq.0) then
        call Multiply_Dirac_to_chi(temperature,xmat_row_smeared,xmat_column_smeared,alpha,&
            chi,mchi,GAMMA10d,nbmn,flux,myrank)
        !*************************************************
        !*** Mchi(kblock,jblock) = Mchi_column(kblock) ***
        !*************************************************
        kblock_send=iblock
        kblock_rcv=iblock
        do ishift=1,nblock-1
            kblock_send=kblock_send+1
            kblock_rcv=kblock_rcv-1
            if(kblock_send.EQ.nblock+1)then
                kblock_send=1
            end if
            if(kblock_rcv.EQ.0)then
                kblock_rcv=nblock
            end if
            send_rank=(isublat-1)*nblock*nblock+(kblock_send-1)*nblock+jblock-1
            receive_rank=(isublat-1)*nblock*nblock+(kblock_rcv-1)*nblock+jblock-1
            tag=4
            call MPI_Isend(Mchi(1,1,1,1,1),&
                &nmat_block*nmat_block*nsite_local*nspin*nremez_md,&
                &MPI_DOUBLE_COMPLEX,&
                &send_rank,tag,MPI_COMM_WORLD,ireq,ierr)
            call MPI_Recv(Mchi_column(1,1,1,1,1,kblock_rcv),&
                &nmat_block*nmat_block*nsite_local*nspin*nremez_md,&
                &MPI_DOUBLE_COMPLEX,&
                &receive_rank,tag,MPI_COMM_WORLD,status,ierr)
            call MPI_Wait(ireq,status,ierr)
        end do
        do iremez=1,nremez_md
            !$omp parallel
            !$omp do
            do imat=1,nmat_block
                do jmat=1,nmat_block
                    do ispin=1,nspin
                        do isite=1,nsite_local
                            Mchi_column(iremez,imat,jmat,ispin,isite,iblock)&
                                &=Mchi(iremez,imat,jmat,ispin,isite)
                        end do
                    end do
                end do
            end do
        !$omp end do
        !$omp end parallel
        end do
        !**********************************************
        !*** Mchi(iblock,kblock) = Mchi_row(kblock) ***
        !**********************************************
        kblock_send=jblock
        kblock_rcv=jblock
        do ishift=1,nblock-1
            kblock_send=kblock_send+1
            kblock_rcv=kblock_rcv-1
            if(kblock_send.EQ.nblock+1)then
                kblock_send=1
            end if
            if(kblock_rcv.EQ.0)then
                kblock_rcv=nblock
            end if
            send_rank=(isublat-1)*nblock*nblock+(iblock-1)*nblock+kblock_send-1
            receive_rank=(isublat-1)*nblock*nblock+(iblock-1)*nblock+kblock_rcv-1
            tag=5
            call MPI_Isend(Mchi(1,1,1,1,1),&
                &nmat_block*nmat_block*nsite_local*nspin*nremez_md,&
                &MPI_DOUBLE_COMPLEX,&
                &send_rank,tag,MPI_COMM_WORLD,ireq,ierr)
            call MPI_Recv(Mchi_row(1,1,1,1,1,kblock_rcv),&
                &nmat_block*nmat_block*nsite_local*nspin*nremez_md,&
                &MPI_DOUBLE_COMPLEX,&
                &receive_rank,tag,MPI_COMM_WORLD,status,ierr)
            call MPI_Wait(ireq,status,ierr)
        end do
        do iremez=1,nremez_md
            !$omp parallel
            !$omp do
            do imat=1,nmat_block
                do jmat=1,nmat_block
                    do ispin=1,nspin
                        do isite=1,nsite_local
                            Mchi_row(iremez,imat,jmat,ispin,isite,jblock)&
                                &=Mchi(iremez,imat,jmat,ispin,isite)
                        end do
                    end do
                end do
            end do
        !$omp end do
        !$omp end parallel
        end do
        !***********************************************
        !*** chi(kblock,iblock) = chi_column(kblock) ***
        !***********************************************
        kblock_send=iblock
        kblock_rcv=jblock
        do ishift=1,nblock-1
            kblock_send=kblock_send+1
            kblock_rcv=kblock_rcv-1
            if(kblock_send.EQ.nblock+1)then
                kblock_send=1
            end if
            if(kblock_rcv.EQ.0)then
                kblock_rcv=nblock
            end if
            send_rank=(isublat-1)*nblock*nblock+(jblock-1)*nblock+kblock_send-1
            receive_rank=(isublat-1)*nblock*nblock+(kblock_rcv-1)*nblock+iblock-1
            tag=6
            call MPI_Isend(chi(1,1,1,1,1),&
                &nmat_block*nmat_block*nsite_local*nspin*nremez_md,&
                &MPI_DOUBLE_COMPLEX,&
                &send_rank,tag,MPI_COMM_WORLD,ireq,ierr)
            call MPI_Recv(chi_column(1,1,1,1,1,kblock_rcv),&
                &nmat_block*nmat_block*nsite_local*nspin*nremez_md,&
                &MPI_DOUBLE_COMPLEX,&
                &receive_rank,tag,MPI_COMM_WORLD,status,ierr)
            call MPI_Wait(ireq,status,ierr)
        end do
        if(iblock.NE.jblock)then
            send_rank=(isublat-1)*nblock*nblock+(jblock-1)*nblock+iblock-1
            receive_rank=(isublat-1)*nblock*nblock+(jblock-1)*nblock+iblock-1
            tag=7
            call MPI_Isend(chi(1,1,1,1,1),&
                &nmat_block*nmat_block*nsite_local*nspin*nremez_md,&
                &MPI_DOUBLE_COMPLEX,&
                &send_rank,tag,MPI_COMM_WORLD,ireq,ierr)
            call MPI_Recv(chi_column(1,1,1,1,1,jblock),&
                &nmat_block*nmat_block*nsite_local*nspin*nremez_md,&
                &MPI_DOUBLE_COMPLEX,&
                &receive_rank,tag,MPI_COMM_WORLD,status,ierr)
            call MPI_Wait(ireq,status,ierr)
        else
            do iremez=1,nremez_md
                !$omp parallel
                !$omp do
                do imat=1,nmat_block
                    do jmat=1,nmat_block
                        do ispin=1,nspin
                            do isite=1,nsite_local
                                chi_column(iremez,imat,jmat,ispin,isite,jblock)&
                                    &=chi(iremez,imat,jmat,ispin,isite)
                            end do
                        end do
                    end do
                end do
            !$omp end do
            !$omp end parallel
            end do
        end if
        !********************************************
        !*** chi(jblock,kblock) = chi_row(kblock) ***
        !********************************************
        kblock_send=jblock
        kblock_rcv=iblock
        do ishift=1,nblock-1
            kblock_send=kblock_send+1
            kblock_rcv=kblock_rcv-1
            if(kblock_send.EQ.nblock+1)then
                kblock_send=1
            end if
            if(kblock_rcv.EQ.0)then
                kblock_rcv=nblock
            end if
            send_rank=(isublat-1)*nblock*nblock+(kblock_send-1)*nblock+iblock-1
            receive_rank=(isublat-1)*nblock*nblock+(jblock-1)*nblock+kblock_rcv-1
            tag=8
            call MPI_Isend(chi(1,1,1,1,1),&
                &nmat_block*nmat_block*nsite_local*nspin*nremez_md,&
                &MPI_DOUBLE_COMPLEX,&
                &send_rank,tag,MPI_COMM_WORLD,ireq,ierr)
            call MPI_Recv(chi_row(1,1,1,1,1,kblock_rcv),&
                &nmat_block*nmat_block*nsite_local*nspin*nremez_md,&
                &MPI_DOUBLE_COMPLEX,&
                &receive_rank,tag,MPI_COMM_WORLD,status,ierr)
            call MPI_Wait(ireq,status,ierr)
        end do
        if(iblock.NE.jblock)then
            send_rank=(isublat-1)*nblock*nblock+(jblock-1)*nblock+iblock-1
            receive_rank=(isublat-1)*nblock*nblock+(jblock-1)*nblock+iblock-1
            tag=9
            call MPI_Isend(chi(1,1,1,1,1),&
                &nmat_block*nmat_block*nsite_local*nspin*nremez_md,&
                &MPI_DOUBLE_COMPLEX,&
                &send_rank,tag,MPI_COMM_WORLD,ireq,ierr)
            call MPI_Recv(chi_row(1,1,1,1,1,iblock),&
                &nmat_block*nmat_block*nsite_local*nspin*nremez_md,&
                &MPI_DOUBLE_COMPLEX,&
                &receive_rank,tag,MPI_COMM_WORLD,status,ierr)
            call MPI_Wait(ireq,status,ierr)
        else
            do iremez=1,nremez_md
                !$omp parallel
                !$omp do
                do imat=1,nmat_block
                    do jmat=1,nmat_block
                        do ispin=1,nspin
                            do isite=1,nsite_local
                                chi_row(iremez,imat,jmat,ispin,isite,iblock)&
                                    &=chi(iremez,imat,jmat,ispin,isite)
                            end do
                        end do
                    end do
                end do
            !$omp end do
            !$omp end parallel
            end do
        end if
        !##################################################################
        !### Now we calculate xmat-derivative of pseudo fermion action. ###
        !### Let us start with (M*chi)^dagger*(dM/dX)*chi.              ###
        !###     delPF_xmat = ((dM/dX)*chi)^dagger*(M*chi)              ###
        !##################################################################
        delPF_xmat_smeared=(0d0,0d0)
        !$omp parallel
        !$omp do
        do kmat=1,nmat_block
            do lmat=1,nmat_block
                do idim=1,ndim
                    do isite=1,nsite_local
                        do iremez=1,nremez_md
                            do imat=1,nmat_block
                                do ispin=1,nspin
                                    do jspin=1,nspin
                                        do kblock=1,nblock
                                            delPF_xmat_smeared(kmat,lmat,idim,isite)=&
                                                &delPF_xmat_smeared(kmat,lmat,idim,isite)&
                                                &-acoeff_md(iremez)*dcmplx(lattice_spacing)&
                                                &*GAMMA10d(idim,ispin,jspin)&
                                                &*(dconjg(Mchi_column(iremez,imat,lmat,ispin,isite,kblock))&
                                                &*chi_column(iremez,imat,kmat,jspin,isite,kblock)&
                                                &-dconjg(Mchi_row(iremez,kmat,imat,ispin,isite,kblock))&
                                                &*chi_row(iremez,lmat,imat,jspin,isite,kblock))

                                        end do
                                    end do
                                end do
                            end do
                        end do
                    end do
                end do
            end do
        end do
        !$omp end do
        !$omp end parallel
        !###################################################################
        !### Wow I calculated the complex conjugate by mistake, sorry... ###
        !###################################################################
        call smearing_delPF(delPF_xmat_smeared,delPF_xmat,s_smear,myrank,nsmear)
        delPF_xmat=dconjg(delPF_xmat)
        !##################################################################
        !## We have to take into account both                            ##
        !## (M*chi)^dagger*(dM/dX)*chi and ((dM/dX)*chi)^dagger*(M*chi). ##
        !##################################################################
        !**********************************
        !*** (M*chi)^dagger*(dM/dX)*chi ***
        !**********************************
        delh_xmat=delh_xmat+delPF_xmat
        !************************************
        !*** ((dM/dX)*chi)^dagger*(M*chi) ***
        !************************************
        if(iblock.EQ.jblock)then
            !$omp parallel
            !$omp do
            do kmat=1,nmat_block
                do lmat=1,nmat_block
                    do idim=1,ndim
                        do isite=1,nsite_local
                            delh_xmat(kmat,lmat,idim,isite)=&
                                delh_xmat(kmat,lmat,idim,isite)&
                                +dconjg(delPF_xmat(lmat,kmat,idim,isite))
                        end do
                    end do
                end do
            end do
        !$omp end do
        !$omp end parallel
        else
            send_rank=(isublat-1)*nblock*nblock+(jblock-1)*nblock+iblock-1
            receive_rank=(isublat-1)*nblock*nblock+(jblock-1)*nblock+iblock-1

            tag=10
            call MPI_Isend(delPF_xmat(1,1,1,1),&
                nmat_block*nmat_block*nsite_local*ndim,&
                &MPI_DOUBLE_COMPLEX,&
                &send_rank,tag,MPI_COMM_WORLD,ireq,ierr)
            call MPI_Recv(mat_rcv_2(1,1,1,1),&
                nmat_block*nmat_block*nsite_local*ndim,&
                &MPI_DOUBLE_COMPLEX,&
                &receive_rank,tag,MPI_COMM_WORLD,status,ierr)
            call MPI_Wait(ireq,status,ierr)
            !$omp parallel
            !$omp do
            do kmat=1,nmat_block
                do lmat=1,nmat_block
                    do idim=1,ndim
                        do isite=1,nsite_local
                            delh_xmat(kmat,lmat,idim,isite)=&
                                delh_xmat(kmat,lmat,idim,isite)&
                                +dconjg(mat_rcv_2(lmat,kmat,idim,isite))
                        end do
                    end do
                end do
            end do
        !$omp end do
        !$omp end parallel
        end if
        !######################################################
        !######################################################
        !## Now we have finished calculating xmat-derivarive ##
        !## of the pseudofermion action.                     ##
        !## Next we move to to the alpha-derivative of       ##
        !## the pseudofermion action.                        ##
        !######################################################
        !######################################################
        !(dM/d¥alpha)*chi
        if(ngauge.eq.0)then
            call Derivative_Dirac_alpha(alpha,chi,Deriv_Mchi_alpha,myrank)
            !-acoeff*(M*chi)^dagger*(dM/alpha)*chi
            delPF_alpha=0d0
            do iremez=1,nremez_md
                !$omp parallel
                !$omp do
                do kmat=(iblock-1)*nmat_block+1,iblock*nmat_block
                    do isite=1,nsite_local
                        do imat=1,nmat_block
                            do jmat=1,nmat_block
                                do ispin=1,nspin
                                    delPF_alpha(kmat)=delPF_alpha(kmat)&
                                        -acoeff_md(iremez)&
                                        *dble(dconjg(Mchi(iremez,imat,jmat,ispin,isite))&
                                        *Deriv_Mchi_alpha(iremez,kmat,imat,jmat,ispin,isite))
                                end do
                            end do
                        end do
                    end do
                end do
                !$omp end do
                !$omp end parallel
                if(jblock.NE.iblock)then
                    !$omp parallel
                    !$omp do
                    do kmat=(jblock-1)*nmat_block+1,jblock*nmat_block
                        do isite=1,nsite_local
                            do imat=1,nmat_block
                                do jmat=1,nmat_block
                                    do ispin=1,nspin
                                        delPF_alpha(kmat)=delPF_alpha(kmat)&
                                            -acoeff_md(iremez)&
                                            *dble(dconjg(Mchi(iremez,imat,jmat,ispin,isite))&
                                            *Deriv_Mchi_alpha(iremez,kmat,imat,jmat,ispin,isite))
                                    end do
                                end do
                            end do
                        end do
                    end do
                   !$omp end do
                   !$omp end parallel
                end if
            end do
            !take into account both
            !(M*chi)^dagger*(dM/d¥alpha)*chi and ((dM/d¥alpha)*chi)^dagger*(M*chi)
            !$omp parallel
            !$omp do
            do kmat=1,nmat_block*nblock
                delh_alpha(kmat)=delh_alpha(kmat)+delPF_alpha(kmat)*2d0
            end do
           !$omp end do
           !$omp end parallel
        end if
    ! if pruebosonic.eq.0
    end if
    !###########################################
    !### Pseudofermion part is finally over. ###
    !###########################################
    !****************************
    !****************************
    !*** constraint for alpha ***
    !****************************
    !****************************
    if(ngauge.eq.0)then
        if(myrank.EQ.0)then
            imax=1
            imin=1
            alpha_max=alpha(1)
            alpha_min=alpha(1)
            do imat=2,nmat_block*nblock
                if(alpha(imat).GT.alpha_max)then
                    imax=imat
                    alpha_max=alpha(imat)
                else if(alpha(imat).LT.alpha_min)then
                    imin=imat
                    alpha_min=alpha(imat)
                end if
            end do
            if(alpha_max-alpha_min.LT.2d0*pi)then
                delh_alpha(imax)=delh_alpha(imax)&
                    +1d0/(2d0*pi-(alpha_max-alpha_min)+1d0/gcoeff_alpha)
                delh_alpha(imin)=delh_alpha(imin)&
                    -1d0/(2d0*pi-(alpha_max-alpha_min)+1d0/gcoeff_alpha)
            else if(alpha_max-alpha_min.GE.2d0*pi)then
                delh_alpha(imax)=delh_alpha(imax)+gcoeff_alpha
                delh_alpha(imin)=delh_alpha(imin)-gcoeff_alpha
            end if

        end if
    end if
    !****************************
    !****************************
    !*** constraint for TrX^2 ***
    !****************************
    !****************************
    call Calc_TrX2(xmat,trx2,myrank)
    call MPI_Bcast(trx2,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
    if(trx2.GE.RCUT)then
        do isite=1,nsite_local
            do idim=1,ndim
                !$omp parallel
                !$omp do
                do imat=1,nmat_block
                    do jmat=1,nmat_block
                        delh_xmat(imat,jmat,idim,isite)=&
                            delh_xmat(imat,jmat,idim,isite)&
                            +2d0*g_R*xmat(imat,jmat,idim,isite)&
                            /dble(nsite_local*nsublat)
                    end do
                end do
            !$omp end do
            !$omp end parallel
            end do
        end do
    end if
    !******************************
    !******************************
    !*** Plane wave deformation ***
    !******************************
    !******************************
    if(nbmn.EQ.1)then
        !***************************
        !*** deriv. of mass term ***
        !***************************
        do isite=1,nsite_local
            !$omp parallel
            !$omp do
            do imat=1,nmat_block
                do jmat=1,nmat_block
                    temp=dcmplx(flux*flux*lattice_spacing*dble(nmat_block*nblock))
                    do idim=1,3
                        delh_xmat(imat,jmat,idim,isite)=&
                            delh_xmat(imat,jmat,idim,isite)&
                            +xmat(imat,jmat,idim,isite)*temp
                    end do
                    temp=dcmplx(0.25d0*flux*flux*lattice_spacing&
                        &*dble(nmat_block*nblock))
                    do idim=4,9
                        delh_xmat(imat,jmat,idim,isite)=&
                            delh_xmat(imat,jmat,idim,isite)&
                            +xmat(imat,jmat,idim,isite)*temp
                    end do
                end do
            end do
        !$omp end do
        !$omp end parallel
        end do
        !******************
        !*** cubic term ***
        !******************
        do isite=1,nsite_local
            com12=(0d0,0d0)
            com23=(0d0,0d0)
            com31=(0d0,0d0)
            !$omp parallel
            !$omp do
            do imat=1,nmat_block
                do jmat=1,nmat_block
                    do kmat=1,nmat_block*nblock
                        com12(imat,jmat)=com12(imat,jmat)&
                            &+xmat_row(imat,kmat,1,isite)&
                            &*xmat_column(kmat,jmat,2,isite)&
                            &-xmat_row(imat,kmat,2,isite)&
                            &*xmat_column(kmat,jmat,1,isite)
                        com23(imat,jmat)=com23(imat,jmat)&
                            &+xmat_row(imat,kmat,2,isite)&
                            &*xmat_column(kmat,jmat,3,isite)&
                            &-xmat_row(imat,kmat,3,isite)&
                            &*xmat_column(kmat,jmat,2,isite)
                        com31(imat,jmat)=com31(imat,jmat)&
                            &+xmat_row(imat,kmat,3,isite)&
                            &*xmat_column(kmat,jmat,1,isite)&
                            &-xmat_row(imat,kmat,1,isite)&
                            &*xmat_column(kmat,jmat,3,isite)
                    end do
                end do
            end do
            !$omp end do
            !$omp end parallel
            temp=(0d0,3d0)*dcmplx(flux*lattice_spacing*dble(nmat_block*nblock))
            !$omp parallel
            !$omp do
            do imat=1,nmat_block
                do jmat=1,nmat_block
                    delh_xmat(imat,jmat,1,isite)=&
                        delh_xmat(imat,jmat,1,isite)+com23(imat,jmat)*temp
                    delh_xmat(imat,jmat,2,isite)=&
                        delh_xmat(imat,jmat,2,isite)+com31(imat,jmat)*temp
                    delh_xmat(imat,jmat,3,isite)=&
                        delh_xmat(imat,jmat,3,isite)+com12(imat,jmat)*temp
                end do
            end do
        !$omp end do
        !$omp end parallel
        end do

    end if

    !****************************
    !*** fixing polyakov loop ***
    !****************************
    if((ngauge.EQ.0).AND.(g_pol.NE.0))then
        if(myrank.EQ.0)then
            !call Calc_Polyakov(nmat_block*nblock,alpha,Pol)
            re_pol=0d0
            im_pol=0d0
            do imat=1,nmat_block*nblock
                re_pol=re_pol+dcos(alpha(imat))
                im_pol=im_pol+dsin(alpha(imat))
            end do
            re_pol=re_pol/dble(nmat_block*nblock)
            im_pol=im_pol/dble(nmat_block*nblock)
            pol=dsqrt(re_pol**2d0+im_pol**2d0)
            if(pol.GT.(polfix+0.5d0*pol_fix_width))then
                !ham=ham+0.5d0*g_pol*(pol-(polfix+0.5d0*pol_fix_width))**2d0
                do imat=1,nmat_block*nblock
                    delh_alpha(imat)=delh_alpha(imat)&
                        &+g_pol*(pol-(polfix+0.5d0*pol_fix_width))&
                        &/dble(nmat_block*nblock)&
                        &*(-re_pol*dsin(alpha(imat))+im_pol*dcos(alpha(imat)))
                end do
            else if(pol.LT.(polfix-0.5d0*pol_fix_width))then
                !ham=ham+0.5d0*g_pol*(pol-(polfix-0.5d0*pol_fix_width))**2d0
                do imat=1,nmat_block*nblock
                    delh_alpha(imat)=delh_alpha(imat)&
                        &+g_pol*(pol-(polfix-0.5d0*pol_fix_width))&
                        &/dble(nmat_block*nblock)&
                        &*(-re_pol*dsin(alpha(imat))+im_pol*dcos(alpha(imat)))
                end do
            end if
        end if
    end if
          !****************************
      !*** constraint for Myers observable ***
      !****************************
    if(g_myers.NE.0)then
        call Calc_Myers(xmat,myers,myrank)
        ! factor 1/3 included in definition compared to action term
        ! myers_local=dble((0d0,-1d0)*(trx123-trx132))/dble(nmat_block*nblock)/dble(nsite_local*nsublat)
        temp=dcmplx(0d0,-1d0/dble(nmat_block*nblock)/dble(nsite_local*nsublat))
        if(myers.GT.(myersfix+0.5d0*myers_fix_width))then
            !ham=ham+0.5d0*g_myers*(myers-(myersfix+0.5d0*myers_fix_width))**2d0
           temp=temp*dcmplx(g_myers*(myers-(myersfix+0.5d0*myers_fix_width)),0d0)
        else if(myers.LT.(myersfix-0.5d0*myers_fix_width))then
            !ham=ham+0.5d0*g_myers*(myers-(myersfix-0.5d0*myers_fix_width))**2d0
           temp=temp*dcmplx(g_myers*(myers-(myersfix-0.5d0*myers_fix_width)),0d0)
        end if

        do isite=1,nsite_local
            ! this re-computation could be avoided.
            com12=(0d0,0d0)
            com23=(0d0,0d0)
            com31=(0d0,0d0)
            do imat=1,nmat_block
                do jmat=1,nmat_block
                    do kmat=1,nmat_block*nblock
                        com12(imat,jmat)=com12(imat,jmat)&
                            &+xmat_row(imat,kmat,1,isite)&
                            &*xmat_column(kmat,jmat,2,isite)&
                            &-xmat_row(imat,kmat,2,isite)&
                            &*xmat_column(kmat,jmat,1,isite)
                        com23(imat,jmat)=com23(imat,jmat)&
                            &+xmat_row(imat,kmat,2,isite)&
                            &*xmat_column(kmat,jmat,3,isite)&
                            &-xmat_row(imat,kmat,3,isite)&
                            &*xmat_column(kmat,jmat,2,isite)
                        com31(imat,jmat)=com31(imat,jmat)&
                            &+xmat_row(imat,kmat,3,isite)&
                            &*xmat_column(kmat,jmat,1,isite)&
                            &-xmat_row(imat,kmat,1,isite)&
                            &*xmat_column(kmat,jmat,3,isite)
                    end do
                end do
            end do

            do imat=1,nmat_block
                do jmat=1,nmat_block
                    delh_xmat(imat,jmat,1,isite)=&
                        delh_xmat(imat,jmat,1,isite)+com23(imat,jmat)*temp
                    delh_xmat(imat,jmat,2,isite)=&
                        delh_xmat(imat,jmat,2,isite)+com31(imat,jmat)*temp
                    delh_xmat(imat,jmat,3,isite)=&
                        delh_xmat(imat,jmat,3,isite)+com12(imat,jmat)*temp
                end do
            end do
        end do
    end if

    return

END SUBROUTINE Calc_Force
               
