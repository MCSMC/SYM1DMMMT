subroutine Calc_myers(xmat,myers,myrank)
    implicit none
  include 'mpif.h'
  include 'size_parallel.h'
    !***** input *****
    integer myrank
    double complex xmat(1:nmat_block,1:nmat_block,1:ndim,&
        -(nmargin-1):nsite_local+nmargin)
    !***** output *****
    double precision myers
    !******************
    double precision myers_local

    integer isite!,isite_p1
    integer idim,jdim
    integer imat,jmat,kmat

    double complex xmat_row(1:nmat_block,1:nmat_block*nblock,1:ndim,1:nsite_local)
    double complex xmat_column(1:nmat_block*nblock,1:nmat_block,&
        &1:ndim,1:nsite_local)
    double complex x23(1:nmat_block,1:nmat_block),&
        x32(1:nmat_block,1:nmat_block),&
        trx123,trx132
    integer IERR

    !call who_am_i(myrank,isublat,iblock,jblock)
    !move i-th row and j-th row of xmat to (i,j)-th node.
    call mpi_xmat_row(xmat,xmat_row,myrank)
    call mpi_xmat_column(xmat,xmat_column,myrank)

    trx123=(0d0,0d0)
    trx132=(0d0,0d0)
    do isite=1,nsite_local
        x23=(0d0,0d0)
        x32=(0d0,0d0)
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
        do imat=1,nmat_block
            do jmat=1,nmat_block
                trx123=trx123+x23(imat,jmat)*dconjg(xmat(imat,jmat,1,isite))
                trx132=trx132+x32(imat,jmat)*dconjg(xmat(imat,jmat,1,isite))
            end do
        end do
    end do
    myers_local=dble((0d0,-1d0)*(trx123-trx132))/dble(nmat_block*nblock)/dble(nsite_local*nsublat)

!    call MPI_Reduce(myers_local,myers,1,MPI_DOUBLE_PRECISION,&
!        MPI_SUM,0,MPI_COMM_WORLD,IERR)
    call MPI_Allreduce(myers_local,myers,1,MPI_DOUBLE_PRECISION,&
        MPI_SUM,MPI_COMM_WORLD,IERR)


    return

END subroutine Calc_myers
