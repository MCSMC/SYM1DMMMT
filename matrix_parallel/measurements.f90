SUBROUTINE measurements(xmat,alpha,nbc,nbmn,myrank,temperature,flux,&
     &GAMMA10d,neig_max,neig_min,ham_init,ham_fin,itraj,ntrial,iteration,&
     &max_err,max_iteration,ncv,n_bad_CG,nacceptance,nsmear,s_smear,&
     &acoeff_md,bcoeff_md,acoeff_pf,bcoeff_pf,ngauge,purebosonic)

  implicit none
  include 'size_parallel.h'
  include '../unit_number.inc'
  !input
  integer ngauge,purebosonic
  double complex xmat(1:nmat_block,1:nmat_block,1:ndim,&
       &-(nmargin-1):nsite_local+nmargin)
  double precision alpha(1:nmat_block*nblock)
  integer nbc,nbmn,myrank,ntrial,iteration,itraj,nacceptance,nsmear
  double precision s_smear
  integer ncv,n_bad_CG
  integer neig_max,neig_min
  double precision temperature,flux
  double complex Gamma10d(1:ndim,1:nspin,1:nspin)
  double precision max_err
  integer max_iteration

  double precision sum_trx2,trx2(1:ndim),com2,Pol,largest_eig,smallest_eig,&
       &energy,ham_fin,ham_init
  double complex xmat_smeared(1:nmat_block,1:nmat_block,1:ndim,&
       &-(nmargin-1):nsite_local+nmargin)
  double precision acoeff_md(0:nremez_md),bcoeff_md(1:nremez_md)
  double precision acoeff_pf(0:nremez_pf),bcoeff_pf(1:nremez_pf)

  double precision myers
  
  !measurements. 
  call Calc_TrX2_each(xmat,sum_trx2,myrank,trx2)
  call Calc_Com2(xmat,com2,myrank)
  call Calc_energy(temperature,xmat,alpha,energy,myers,myrank,nbmn,flux,&
       &acoeff_md,bcoeff_md,acoeff_pf,bcoeff_pf,&
       &nbc,max_err,max_iteration,ngauge,purebosonic)
    
  if(myrank.EQ.0)then
     call Calc_Polyakov(nmat_block*nblock,alpha,Pol)
  end if
  !largest eigenvalue of D^dag*D. 
  if(neig_max.GT.0)then
          !smear
     call smearing_xmat(xmat,xmat_smeared,s_smear,myrank,nsmear)
     call Largest_eigenvalue(temperature,&
          &xmat_smeared,alpha,GAMMA10d,neig_max,largest_eig,myrank,nbc,nbmn,flux)
  end if
  !smallest eigenvalue of D^dag*D. 
  if(neig_min.GT.0)then
     !smear
     call smearing_xmat(xmat,xmat_smeared,s_smear,myrank,nsmear)
     call Smallest_eigenvalue(temperature,&
          &xmat_smeared,alpha,GAMMA10d,neig_min,smallest_eig,myrank,nbc,&
          &max_err,max_iteration,nbmn,flux)
  end if

  !**************
  !*** output ***
  !**************
  if(myrank.EQ.0)then
     if((neig_max.EQ.0).AND.(neig_min.EQ.0))then
        
        write(unit_measurement,'(I8,1x,f15.9,1x,I4,1x,I8,1x,I8,15(1x,f15.7))')&
             &itraj,ham_fin-ham_init,ncv,n_bad_CG,iteration,energy,&
             &Pol,sum_trx2,&
             &trx2(1),trx2(2),trx2(3),trx2(4),trx2(5),trx2(6),trx2(7),trx2(8),trx2(9),&
             &com2,myers,dble(nacceptance)/dble(ntrial)
        
        write(*,'(I8,1x,f15.9,1x,I4,1x,I8,1x,I8,15(1x,f15.7))')&
             &itraj,ham_fin-ham_init,ncv,n_bad_CG,iteration,energy,&
             &Pol,sum_trx2,&
             &trx2(1),trx2(2),trx2(3),trx2(4),trx2(5),trx2(6),trx2(7),trx2(8),trx2(9),&
             &com2,myers,dble(nacceptance)/dble(ntrial)
        
     else  if((neig_max.GT.0).AND.(neig_min.EQ.0))then
        
        write(unit_measurement,'(I8,1x,f15.9,1x,I4,1x,I8,1x,I8,16(1x,f15.7))')&
             &itraj,ham_fin-ham_init,ncv,n_bad_CG,iteration,energy,&
             &Pol,sum_trx2,&
             &trx2(1),trx2(2),trx2(3),trx2(4),trx2(5),trx2(6),trx2(7),trx2(8),trx2(9),&             
             &com2,myers,largest_eig,&
             &dble(nacceptance)/dble(ntrial)
        write(*,'(I8,1x,f15.9,1x,I4,1x,I8,1x,I8,16(1x,f15.7))')&
             &itraj,ham_fin-ham_init,ncv,n_bad_CG,iteration,energy,&
             &Pol,sum_trx2,&
             &trx2(1),trx2(2),trx2(3),trx2(4),trx2(5),trx2(6),trx2(7),trx2(8),trx2(9),&
             &com2,myers,largest_eig,&
             &dble(nacceptance)/dble(ntrial)
        
     else  if((neig_max.EQ.0).AND.(neig_min.GT.0))then
        
        write(unit_measurement,'(I8,1x,f15.9,1x,I4,1x,I8,1x,I8,16(1x,f15.7))')&
             &itraj,ham_fin-ham_init,ncv,n_bad_CG,iteration,energy,&
             &Pol,sum_trx2,&
             &trx2(1),trx2(2),trx2(3),trx2(4),trx2(5),trx2(6),trx2(7),trx2(8),trx2(9),&
             &com2,myers,smallest_eig,&
             &dble(nacceptance)/dble(ntrial)
        write(*,'(I8,1x,f15.9,1x,I4,1x,I8,1x,I8,16(1x,f15.7))')&
             &itraj,ham_fin-ham_init,ncv,n_bad_CG,iteration,energy,&
             &Pol,sum_trx2,&
             &trx2(1),trx2(2),trx2(3),trx2(4),trx2(5),trx2(6),trx2(7),trx2(8),trx2(9),&
             &com2,myers,smallest_eig,&
             &dble(nacceptance)/dble(ntrial)
        
     else  if((neig_max.GT.0).AND.(neig_min.GT.0))then
        write(unit_measurement,'(I8,1x,f15.9,1x,I4,1x,I8,1x,I8,17(1x,f15.7))')&
             &itraj,ham_fin-ham_init,ncv,n_bad_CG,iteration,energy,&
             &Pol,sum_trx2,&
             &trx2(1),trx2(2),trx2(3),trx2(4),trx2(5),trx2(6),trx2(7),trx2(8),trx2(9),&
             &com2,myers,largest_eig,smallest_eig,&
             &dble(nacceptance)/dble(ntrial)
        write(*,'(I8,1x,f15.9,1x,I4,1x,I8,1x,I8,17(1x,f15.7))')&
             &itraj,ham_fin-ham_init,ncv,n_bad_CG,iteration,energy,&
             &Pol,sum_trx2,&
             &trx2(1),trx2(2),trx2(3),trx2(4),trx2(5),trx2(6),trx2(7),trx2(8),trx2(9),&
             &com2,myers,largest_eig,smallest_eig,&
             &dble(nacceptance)/dble(ntrial)
     end if

     write(unit_Polyakov_phase,*)alpha
     
  end if

  return
 
END SUBROUTINE measurements
