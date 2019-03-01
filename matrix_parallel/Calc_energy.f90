subroutine Calc_energy(temperature,xmat,alpha,energy,myers,myrank,nbmn,flux,&
     &acoeff_md,bcoeff_md,acoeff_pf,bcoeff_pf,&
     &nbc,max_err,max_iteration,ngauge,purebosonic)

  implicit none

  include 'mpif.h'
  include 'size_parallel.h'
  !***** input *****
  integer myrank,nbmn
  integer ngauge,purebosonic
  double precision temperature,flux
  double complex xmat(1:nmat_block,1:nmat_block,1:ndim,&
       -(nmargin-1):nsite_local+nmargin)
  double precision alpha(1:nmat_block*nblock)
  double precision acoeff_md(0:nremez_md),bcoeff_md(1:nremez_md)
  double precision acoeff_pf(0:nremez_pf),bcoeff_pf(1:nremez_pf)
  integer nbc,max_iteration
  double precision max_err
  !***** output *****
  double precision energy,myers
  !******************
  double precision action,kinetic,potential,potential_BMN,energy_local,myers_local
  !double precision lattice_spacing

  double complex commutator(1:nmat_block,1:nmat_block)
  double complex uxumx(1:nmat_block,1:nmat_block)
  integer isite!,isite_p1
  integer idim,jdim
  integer imat,jmat,kmat
  double complex ei,ej

  double complex xmat_row(1:nmat_block,1:nmat_block*nblock,1:ndim,1:nsite_local)
  double complex xmat_column(1:nmat_block*nblock,1:nmat_block,&
       &1:ndim,1:nsite_local)
  double complex x23(1:nmat_block,1:nmat_block),&
       x32(1:nmat_block,1:nmat_block),&
       trx123,trx132
  double precision trx2_123,trx2_456789

  
  integer i,kblock_send,kblock_rcv,ishift,iblock,jblock,isublat,info
  double precision sum_pf
  
  
  double complex Gamma10d(1:ndim,1:nspin,1:nspin),gam123(1:nspin,1:nspin),gam12(1:nspin,1:nspin)
  integer ispin,jspin,kspin
  integer iremez
  integer iteration
  !***** For MPI *****                             
  double complex mat_send(1:nmat_block,1:nmat_block,1:nspin,1:nsite_local)
  double complex mat_rcv(1:nmat_block,1:nmat_block,1:nspin,1:nsite_local)
  integer IERR,IREQ,send_rank,receive_rank,tag
  integer STATUS(MPI_STATUS_SIZE)


  call who_am_i(myrank,isublat,iblock,jblock)
  !move i-th row and j-th row of xmat to (i,j)-th node.
  call mpi_xmat_row(xmat,xmat_row,myrank)
  call mpi_xmat_column(xmat,xmat_column,myrank)
  !nprocs=nsublat*nmat_block*nmat_block
  !**********************
  !*** potential term ***
  !**********************
  potential=0d0
  do isite=1,nsite_local
     !isite=0, nsite_local+1 are considered at neighboring nodes.
     do idim=1,ndim-1
        do jdim=idim+1,ndim
           commutator=(0d0,0d0)  
           do imat=1,nmat_block
              do jmat=1,nmat_block
                 do kmat=1,nmat_block*nblock
                    commutator(imat,jmat)=commutator(imat,jmat)&
                         &+xmat_row(imat,kmat,idim,isite)&
                         &*xmat_column(kmat,jmat,jdim,isite)&
                         &-xmat_row(imat,kmat,jdim,isite)&
                         &*xmat_column(kmat,jmat,idim,isite)
                 end do
              end do
           end do
           do imat=1,nmat_block
              do jmat=1,nmat_block
                 potential=potential&
                      +dble(commutator(imat,jmat)*dconjg(commutator(imat,jmat)))
              end do
           end do          
        end do
     end do
  end do
  !note that we took some over idim<jdim only. So the coefficient is 0.75*2=1.5.
  potential=potential*1.5d0
  !******************************  
  !*** Plane wave deformation ***
  !******************************
  potential_BMN=0d0
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
     potential_BMN=potential_BMN+flux*flux*(trx2_123+0.25d0*trx2_456789)
     !******************
     !*** cubic term ***
     !******************
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
     potential_BMN=potential_BMN+dble((0d0,7.5d0)*(trx123-trx132))*flux
     myers_local=dble((0d0,-1d0)*(trx123-trx132))/dble(nmat_block*nblock)/dble(nsite_local*nsublat)
  end if

  if(purebosonic.eq.0) then
    call calc_energy_fermion(temperature,xmat,xmat_row,xmat_column,&
       &alpha,GAMMA10d,nbmn,flux,myrank,sum_pf,max_err,max_iteration,nbc)
  end if
  
  energy_local=(potential+potential_BMN)/dble(nmat_block*nblock)/dble(nsite_local*nsublat)&
       &+sum_pf

  call MPI_Reduce(energy_local,energy,1,MPI_DOUBLE_PRECISION,&
       MPI_SUM,0,MPI_COMM_WORLD,IERR)

  call MPI_Reduce(myers_local,myers,1,MPI_DOUBLE_PRECISION,&
       MPI_SUM,0,MPI_COMM_WORLD,IERR)


  return

END subroutine Calc_energy
