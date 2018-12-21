! pf2 = (-i)*D*(lattice spacing) * Pf1
!Please be careful about the overall factor (-i)*(lattice spacing). 
SUBROUTINE calc_energy_fermion(temperature,xmat,xmat_row_2,xmat_column_2,&
     &alpha,GAMMA10d,nbmn,flux,myrank,sum_pf,max_err,max_iteration,nbc)

  implicit none

  include 'mpif.h'
  include 'size_parallel.h'
  !***** input *****
  integer nbmn,myrank,nbc
  double precision temperature,flux
  double precision alpha(1:nmat_block*nblock)
  double complex GAMMA10d(1:ndim,1:nspin,1:nspin)
  double complex xmat(1:nmat_block,1:nmat_block,1:ndim,&
       -(nmargin-1):nsite_local+nmargin)
  double complex xmat_row_2(1:nmat_block,1:nmat_block*nblock,1:ndim,&
       &1:nsite_local)
  double complex xmat_column_2(1:nmat_block*nblock,1:nmat_block,1:ndim,&
       &1:nsite_local)
  double precision max_err
  integer max_iteration
  !***** output *****
  double precision sum_pf
  !******************
  double complex xmat_row(1:nmat_block,1:nmat_block*nblock,1:ndim,1:nsite_local)
  double complex xmat_column(1:nmat_block*nblock,1:nmat_block,1:ndim,&
       &1:nsite_local)
  double complex Gam123(1:nspin,1:nspin),Gam12(1:nspin,1:nspin)
  double complex phase,phase2,GAM
  double precision aij
  double precision lattice_spacing
  integer imat,jmat,kmat
  integer idim
  integer ispin,jspin,kspin
  integer isite
  integer isublat
  integer iblock,jblock,kblock_send,kblock_rcv
  integer ishift

  integer nremez
  double precision bcoeff(1:1)
  !Phi is random Gaussian used for the fermionic part
  double complex phi(1:nmat_block,1:nmat_block,&
       &1:nspin,-(nmargin-1):nsite_local+nmargin)
  !use multimass CG solver to obtain chi, with nremez=1
  double complex chi(1:1,1:nmat_block,1:nmat_block,&
       &1:nspin,-(nmargin-1):nsite_local+nmargin)

  double complex pf1(1:nmat_block,1:nmat_block,&
       1:nspin,-(nmargin-1):nsite_local+nmargin)
 double complex pf2(1:nmat_block,1:nmat_block,&
      1:nspin,-(nmargin-1):nsite_local+nmargin)

 double precision r1,r2
 
 double complex trace_local(1:nsublat,1:nspin,1:nsite_local)
 double complex trace(1:nsublat,1:nspin,1:nsite_local)
 
 integer i,info
 integer iteration
 
  !***** For MPI *****
  double complex mat_send(1:nmat_block,1:nmat_block,1:nspin,1:nsite_local)
  double complex mat_rcv(1:nmat_block,1:nmat_block,1:nspin,1:nsite_local)
  integer IERR,IREQ,send_rank,receive_rank,tag
  integer STATUS(MPI_STATUS_SIZE)

  do idim=1,ndim
     do imat=1,nmat_block
        do jmat=1,nmat_block*nblock
           do isite=1,nsite_local
              xmat_row(imat,jmat,idim,isite)=&
                   &dconjg(xmat_column_2(jmat,imat,idim,isite))
              xmat_column(jmat,imat,idim,isite)=&
                   &dconjg(xmat_row_2(imat,jmat,idim,isite))
           end do
        end do
     end do
  end do
  call who_am_i(myrank,isublat,iblock,jblock)
  lattice_spacing=1d0/temperature/dble(nsite_local*nsublat)
  !********************
  !********************
  !********************
  !*** fermion part ***
  !********************
  !********************
  !********************
!  sum_pf=0d0
  !****************************************************
  !****************************************************
  !*** Firstly we make Gaussian complex matrix Phi. ***
  !****************************************************
  !****************************************************
  
  !**************************************************************************
  !**** we must be careful about the normalization of the Gaussian term. ****
  !**************************************************************************
  if(myrank.GT.0)then
        !throw away some random numbers in order to avoid the communication
     do i=1,nmat_block*nmat_block*nspin*nsite_local*myrank
        call BoxMuller(r1,r2)
     end do
  end if
  
  do imat=1,nmat_block     
     do jmat=1,nmat_block
        do ispin=1,nspin
           do isite=1,nsite_local
              call BoxMuller(r1,r2)
              phi(imat,jmat,ispin,isite)=&
                   &(dcmplx(r1)+dcmplx(r2)*(0D0,1D0))/dcmplx(dsqrt(2d0))
           end do
        end do
     end do
  end do

  if(myrank.LT.nsublat*nblock*nblock-1)then
     !throw away some random numbers in order to avoid the communication
     do i=1,nmat_block*nmat_block*nspin*nsite_local*&
          (nsublat*nblock*nblock-1-myrank)
        call BoxMuller(r1,r2)
     end do
  end if
  !****************************
  !*** traceless projection ***
  !****************************
  call who_am_i(myrank,isublat,iblock,jblock)
  trace_local=(0d0,0d0)
  if(iblock.EQ.jblock)then
     do ispin=1,nspin
        do isite=1,nsite_local
           do imat=1,nmat_block
              trace_local(isublat,ispin,isite)=&
                   &trace_local(isublat,ispin,isite)&
                   &+phi(imat,imat,ispin,isite)
           end do
        end do
     end do
  end if
  
  call MPI_Allreduce(trace_local,trace,nsublat*nspin*nsite_local,&
       MPI_DOUBLE_COMPLEX,&
       &MPI_SUM,MPI_COMM_WORLD,IERR)
  
  trace=trace/dcmplx(nmat_block*nblock)
  
  call who_am_i(myrank,isublat,iblock,jblock)
  if(iblock.EQ.jblock)then
     do ispin=1,nspin
        do isite=1,nsite_local
           do imat=1,nmat_block
              phi(imat,imat,ispin,isite)=phi(imat,imat,ispin,isite)&
                   &-trace(isublat,ispin,isite)
           end do
        end do
     end do
  end if
  !******************************
  !*** adjust margin and b.c. ***
  !******************************
  call Adjust_margin_and_bc_pf(phi,myrank,nbc)
  
  !***********************************************************
  !***********************************************************
  !*** Now Gaussian complex matrix Phi has been generated. ***
  !***********************************************************
  !***********************************************************
  
  nremez=1
  bcoeff(1)=0d0
  !*************************************
  !*** phi -> chi=(M^dag*M)^{-1}*phi ***
  !*************************************
  call MakeGamma(Gamma10d)
  !note that Gamma10d is (-i) times actual gamma matrices. 
  call solver_biCGm(nbc,nremez,bcoeff,temperature,&
       &xmat,alpha,phi,chi,GAMMA10d,max_err,max_iteration,iteration,myrank,nbmn,flux,info)
  !Note that the normalization is still problematic. It will be adjusted later. 
  
  !***********************************
  !***********************************
  !*** Next we calculate pf1=K*Phi ***
  !***********************************
  !***********************************
  pf1=(0d0,0d0)
  !*****************************************************
  !*****************************************************
  !***  interaction part; MPI communication needed.  ***
  !*****************************************************
  !*****************************************************
  
  !###################################################################
  !### A part of xmat_column*pf does not involve MPI communication ###
  !###################################################################     
  do ispin=1,nspin
     do jspin=1,nspin
        do idim=1,ndim
           if(abs(GAMMA10d(idim,ispin,jspin)).GT.1d-2)then
              do isite=1,nsite_local
                 do imat=1,nmat_block
                    do jmat=1,nmat_block
                       do kmat=1,nmat_block
                          pf1(imat,jmat,ispin,isite)=&   
                               &pf1(imat,jmat,ispin,isite)+&
                               &GAMMA10d(idim,ispin,jspin)&
                               &*xmat_column((iblock-1)&
                               &*nmat_block+imat,kmat,idim,isite)&     
                               &*phi(kmat,jmat,jspin,isite)*(0d0,1.5d0)
                       end do
                    end do
                 end do
              end do
           end if
        end do
     end do
  end do
  !############################################     
  !### MPI-communication for xmat_column*pf ###
  !############################################
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
     send_rank=(isublat-1)*nblock*nblock&
          &+(kblock_send-1)*nblock+jblock-1
     receive_rank=(isublat-1)*nblock*nblock&
          &+(kblock_rcv-1)*nblock+jblock-1
     mat_send=(0d0,0d0)
     do ispin=1,nspin
        do jspin=1,nspin
           do idim=1,ndim
              if(abs(GAMMA10d(idim,ispin,jspin)).GT.1d-2)then
                 do isite=1,nsite_local
                    do imat=1,nmat_block
                       do jmat=1,nmat_block
                          do kmat=1,nmat_block
                             mat_send(imat,jmat,ispin,isite)=&
                                  &mat_send(imat,jmat,ispin,isite)+&
                                  &GAMMA10d(idim,ispin,jspin)&
                                  &*xmat_column((kblock_send-1)&
                                  &*nmat_block+imat,kmat,idim,isite)&
                                  &*phi(kmat,jmat,jspin,isite)*(0d0,1.5d0)
                          end do
                       end do
                    end do
                 end do
              end if
           end do
        end do
     end do
     
     tag=1
     call MPI_Isend(mat_send(1,1,1,1),&
          &nmat_block*nmat_block*nsite_local*nspin,&
          &MPI_DOUBLE_COMPLEX,&
          &send_rank,tag,MPI_COMM_WORLD,ireq,ierr)
     call MPI_Recv(mat_rcv(1,1,1,1),&
          &nmat_block*nmat_block*nsite_local*nspin,&
          &MPI_DOUBLE_COMPLEX,&
          &receive_rank,tag,MPI_COMM_WORLD,status,ierr)
     call MPI_Wait(ireq,status,ierr)
     do ispin=1,nspin
        do isite=1,nsite_local
           
           do imat=1,nmat_block
              do jmat=1,nmat_block
                 pf1(imat,jmat,ispin,isite)=&
                      &pf1(imat,jmat,ispin,isite)&
                      &+mat_rcv(imat,jmat,ispin,isite)
              end do
           end do
        end do
     end do
  end do
  !################################################################
  !### A part of pf*xmat_row does not involve MPI communication ###
  !################################################################
  do ispin=1,nspin
     do jspin=1,nspin
        do idim=1,ndim
           if(abs(GAMMA10d(idim,ispin,jspin)).GT.1d-2)then
              do isite=1,nsite_local
                 do imat=1,nmat_block
                    do jmat=1,nmat_block
                       do kmat=1,nmat_block
                          pf1(imat,jmat,ispin,isite)=&     
                               &pf1(imat,jmat,ispin,isite)-&
                               &GAMMA10d(idim,ispin,jspin)&
                               &*phi(imat,kmat,jspin,isite)&
                               &*xmat_row(kmat,(jblock-1)&
                               &*nmat_block+jmat,idim,isite)*(0d0,1.5d0)
                       end do
                    end do
                 end do
              end do
           end if
        end do
     end do
  end do
  !#########################################
  !### MPI-communication for pf*xmat_row ###
  !#########################################
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
     mat_send=(0d0,0d0)
     do ispin=1,nspin
        do jspin=1,nspin   
           do idim=1,ndim
              if(abs(GAMMA10d(idim,ispin,jspin)).GT.1d-2)then
                 do isite=1,nsite_local
                    do imat=1,nmat_block
                       do jmat=1,nmat_block
                          do kmat=1,nmat_block
                             mat_send(imat,jmat,ispin,isite)=&    
                                  &mat_send(imat,jmat,ispin,isite)-&
                                  &GAMMA10d(idim,ispin,jspin)&
                                  &*phi(imat,kmat,jspin,isite)&
                                  &*xmat_row(kmat,(kblock_send-1)&
                                  &*nmat_block+jmat,idim,isite)*(0d0,1.5d0)
                          end do
                       end do
                    end do
                 end do
              end if
           end do
        end do
     end do
     call MPI_Isend(mat_send(1,1,1,1),&
          &nmat_block*nmat_block*nsite_local*nspin,&
          &MPI_DOUBLE_COMPLEX,&
          &send_rank,tag,MPI_COMM_WORLD,ireq,ierr)
     call MPI_Recv(mat_rcv(1,1,1,1),&
          &nmat_block*nmat_block*nsite_local*nspin,&
          &MPI_DOUBLE_COMPLEX,&
          &receive_rank,tag,MPI_COMM_WORLD,status,ierr)
     call MPI_Wait(ireq,status,ierr)
     do ispin=1,nspin
        do isite=1,nsite_local
           do imat=1,nmat_block
              do jmat=1,nmat_block
                 pf1(imat,jmat,ispin,isite)=&
                      &pf1(imat,jmat,ispin,isite)&
                      &+mat_rcv(imat,jmat,ispin,isite)
              end do
           end do 
        end do
     end do
  end do
  !*************************************************************
  !*************************************************************
  !*** Plane wave deformation; no need for MPI communication ***
  !*************************************************************
  !*************************************************************
  if(nbmn.EQ.1)then
     Gam12=(0d0,0d0)
     Gam123=(0d0,0d0)
     do jspin=1,nspin
        do ispin=1,nspin
           do kspin=1,nspin
              Gam12(ispin,jspin)=Gam12(ispin,jspin)&
                   +Gamma10d(1,ispin,kspin)*Gamma10d(2,kspin,jspin)
           end do
        end do
     end do
     do jspin=1,nspin
        do ispin=1,nspin
           do kspin=1,nspin
              Gam123(ispin,jspin)=Gam123(ispin,jspin)&
                   +Gam12(ispin,kspin)*Gamma10d(3,kspin,jspin)
           end do
        end do
     end do
     !So far, Gamma123 is actually (+i) times Gamma123. 
     Gam123=Gam123*dcmplx(flux)*(-0.75d0,0d0)
     do isite=1,nsite_local
        do ispin=1,nspin
           do jspin=1,nspin                                      
              do imat=1,nmat_block
                 do jmat=1,nmat_block
                    pf1(imat,jmat,ispin,isite)=&
                         pf1(imat,jmat,ispin,isite)&
                         +Gam123(ispin,jspin)*phi(imat,jmat,jspin,isite)
                 end do
              end do                        
           end do
        end do
     end do
  end if
  !**********************************
  !*** Now we obtained pf1=K*Phi. ***
  !**********************************

  !*********************************
  !*** Next, pf2=M^dagger*K*Phi. ***
  !*********************************
  !note that we use xmat_row_2 and xmat_column_2 for stupid historical reason.
  call Adjust_margin_and_bc_pf(pf1,myrank,nbc)
  call Multiply_Dirac_dagger(temperature,xmat_row_2,xmat_column_2,&
     &alpha,pf1,pf2,GAMMA10d,nbmn,flux,myrank)
  !************************************
  !*** Inner product of chi and pf2 ***
  !************************************
  sum_pf=0d0
  do isite=1,nsite_local
     do ispin=1,nspin                 
        do imat=1,nmat_block
           do jmat=1,nmat_block
              sum_pf=sum_pf+dble(dconjg(chi(1,imat,jmat,ispin,isite))&
                   &*pf2(imat,jmat,ispin,isite)*(0d0,-1d0))
              !factor (-i) is needed because of the stupid normalization issue
           end do
        end do
     end do
  end do
  !********************************
  !*** Adjust the normalization ***
  !********************************
  sum_pf=sum_pf/dble(nmat_block*nblock*nmat_block*nblock)&
       &*temperature*lattice_spacing&
       &*0.5d0


  !****************************
  !****************************
  !****************************
  !*** fermion part is done ***
  !****************************
  !****************************
  !****************************
  return

END SUBROUTINE Calc_Energy_Fermion
