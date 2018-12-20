!##############################################################################
!######              BFSS matrix model on lattice                     #########
!######                                                               #########
!######                 written by Masanori Hanada                    #########
!######                                                               #########
!######            ver.1  bosonic, HMC, no parallelization.           #########
!######            ver.2  gauge fixed version is available.           #########
!######            ver.3  parallel code.                              #########
!######            ver.4  fermionic part is introduced.               #########
!######                   trace part of the fermion is left as it is. #########
!######            ver.5  U(1) part of pseudo-fermion is removed.     #########
!######                   potential term for alpha, which enforce     #########
!######                   them to satisfy the constraint, is introduced. ######
!######            ver.6  Constraint for alpha is improved.           #########
!######            ver.7  Constraint for TrX^2 is introduced.         #########
!######            ver.8  Fourier acceleration is introduced.         #########
!######                   Does not use Fast Fourier Transform yet.    #########
!######            ver.9  Flux introduced. nbmn=0 -> BFSS, nbmn=1 -> BMN ######
!######            ver.10 Mersenne twistor is introduced.             #########
!######                   (Taken from M.Matsumoto's webpage @ Hiroshima-U.) ###
!######            ver.11 "Fourier_acceleration_optimize" was at wrong place.##
!######                   Fast Fourier Transform is introduced. (A.Tsuchiya) ##
!######                   Debug (flux part of Dirac operator)         #########
!##############################################################################
!######                Parallelization of matrix products.            #########
!######            ver.1  Conservation of Hamiltonian is confirmed.   #########
!######            ver.2  Measurements added, + cleaning up.          #########
!######            ver.3  Dirac operator is modified;                 #########
!######                   Input/output format is modified.            #########
!######            ver.4  Dirac operator is imprroved;                #########
!######                   MPI communication is modified for that purpose. #####
!######            ver.5  Increase the readability. (with I.Kanamori) #########
!######            ver.6  Make it easier to read. (with I.Kanamori)   #########
!######            ver.7  Minor improvement in output format.         #########
!######                   The internal ebergy is now available for BMN.######## 
!######                   Intermediate configurations can be saved.   #########
!######            ver.8  imetropolis=1 -> no Metropolis test         #########
!######            ver.9  Now you can take a log of CG solver.        #########
!######            ver.10 Subtraction of U(1)-part of fermion is corrected. ###
!######            ver.10.1 Multiplication of Dirac has been improved #########
!######                     following LLNL IT-guy's advice.           #########
!######            ver.10.2 Smearing has been introduced. Probably useless. ###
!######            ver.10.3 Both gauged and ungauged can be studied.  #########
!######                     Measurements for BMN added.               #########
!######            ver.10.4 Some bugs in the BMN part removed,        #########
!######                      some measurements added.                 #########
!##############################################################################
!Mersenne twister
include 'mt19937.f90'

program BFSS
  
  use mtmod !Mersenne twistor

  implicit none
  
  include 'matrix_parallel/size_parallel.h'
  include 'mpif.h'
  include 'Fourier.inc'
  include 'unit_number.inc'
  include 'matrix_parallel/include_parallel.h' 
  include 'matrix_parallel/remez_md.dat'
  include 'matrix_parallel/remez_pf.dat'
  !---------------------------------
  call MPI_INIT(IERR)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,NPROCS, IERR)
  call MPI_COMM_RANK(MPI_COMM_WORLD,MYRANK, IERR)
  !*************************
  !**** read parameters ****
  !*************************
  open(unit=unit_input_para,status='OLD',file='input_matrix_v10.dat',&
       &action='READ')
  read(unit_input_para,*) input_config
  read(unit_input_para,*) output_config
  read(unit_input_para,*) data_output
  read(unit_input_para,*) Pol_phase
  read(unit_input_para,*) intermediate_config
  read(unit_input_para,*) acc_input
  read(unit_input_para,*) acc_output
  read(unit_input_para,*) CG_log
  read(unit_input_para,*) nbc
  read(unit_input_para,*) nbmn
  read(unit_input_para,*) ngauge
  read(unit_input_para,*) nsmear
  read(unit_input_para,*) s_smear
  read(unit_input_para,*) init
  read(unit_input_para,*) iaccelerate
  read(unit_input_para,*) isave
  read(unit_input_para,*) temperature
  read(unit_input_para,*) flux
  read(unit_input_para,*) ntraj
  read(unit_input_para,*) nskip
  read(unit_input_para,*) nsave
  read(unit_input_para,*) ntau
  read(unit_input_para,*) dtau_xmat
  read(unit_input_para,*) dtau_alpha
  read(unit_input_para,*) upper_approx
  read(unit_input_para,*) max_err
  read(unit_input_para,*) max_iteration
  read(unit_input_para,*) g_alpha
  read(unit_input_para,*) g_R
  read(unit_input_para,*) RCUT
  read(unit_input_para,*) neig_max
  read(unit_input_para,*) neig_min
  read(unit_input_para,*) nfuzzy
  read(unit_input_para,*) mersenne_seed
  read(unit_input_para,*) imetropolis
  read(unit_input_para,*) purebosonic


  close(unit_input_para)
  !Construc Gamma matrices. 
  call MakeGamma(Gamma10d)
  !initialize the variables for optimizeing Fourier acceleration. 
  imeasure=0
  fluctuation=0d0
  !***************************************
  !*** Rescaling of Remez coefficients ***
  !***************************************
  acoeff_md(0)= acoeff_md(0)*upper_approx**(-0.25d0)
  do iremez=1,nremez_md
     acoeff_md(iremez)=acoeff_md(iremez)*upper_approx**(0.75d0)
     bcoeff_md(iremez)=bcoeff_md(iremez)*upper_approx
  end do
  
  acoeff_pf(0)= acoeff_pf(0)*upper_approx**(0.125d0)
  do iremez=1,nremez_pf
     acoeff_pf(iremez)=acoeff_pf(iremez)*upper_approx**(1.125d0)
     bcoeff_pf(iremez)=bcoeff_pf(iremez)*upper_approx
  end do
  !*************************************
  !*** Set the initial configuration ***
  !*************************************
  call initial_configuration(xmat,alpha,acceleration,itraj,init,iaccelerate,&
       &nfuzzy,input_config,acc_input,flux,mersenne_seed)
  if(ngauge.EQ.1)then
     !ungauged
     alpha=0d0
     dtau_alpha=0d0
  end if
  !***********************************
  !******  Make the output file ******
  !***********************************
  call output_header(myrank,data_output,temperature,flux,&
       &ntau,dtau_xmat,dtaU_alpha,neig_max,neig_min,nbc,nbmn,&
       &init,input_config,output_config,iaccelerate,acc_input,acc_output,&
       &g_alpha,g_R,RCUT,upper_approx,max_err,max_iteration,CG_log,Pol_phase,&
       &isave,nsave,intermediate_config,imetropolis,ngauge,purebosonic)
  
  !*******************************************************
  !******  Make the intermediate configuration file ******
  !*******************************************************
  if(myrank.EQ.0)then
     if(isave.EQ.0)then
        open(unit=unit_intermediate_config,status='REPLACE',&
             &file=intermediate_config,action='WRITE')
     end if
  end if
  
  !**********************************
  !**** initialize counters *********
  !**********************************
  nacceptance=0 !number of acceptance
  ntrial=0 !number of trial 
  ncv=0 !number of the constraint violation.
  n_bad_CG=0 !count how many times CG solver failed to converge.
  
  !*************************************
  !*************************************
  !***  Loop of molecular evolutions ***
  !*************************************
  !*************************************  
  do while (itraj.LE.ntraj)
     !******************************
     !**** Molecular Evolution. ****
     !******************************
     if(myrank.EQ.0)then
        !Take CG_log
        write(unit_CG_log,*)"itraj=",itraj
     end if

     call RHMC_evolution(xmat,alpha,ncv,n_bad_CG,nacceptance,nbc,nbmn,&
          &temperature,flux,GAMMA10d,ntau,dtau_xmat,dtau_alpha,&
     &acceleration,g_alpha,g_R,RCUT,&
     &acoeff_md,bcoeff_md,acoeff_pf,bcoeff_pf,&
     &max_err,max_iteration,iteration,&
     &ham_init,ham_fin,ntrial,imetropolis,nsmear,s_smear,ngauge,purebosonic)

     !######################
     !#### measurements ####
     !######################
     if(MOD(itraj,nskip).EQ.0)then
        !keep the hermiticity with a good precision. 
        !Not really needed but just in case.
        call hermitian_projection(xmat,myrank)
        
        call measurements(xmat,alpha,nbc,nbmn,myrank,temperature,flux,&
             &GAMMA10d,neig_max,neig_min,ham_init,ham_fin,itraj,ntrial,iteration,&
             &max_err,max_iteration,ncv,n_bad_CG,nacceptance,nsmear,s_smear,&
             &acoeff_md,bcoeff_md,acoeff_pf,bcoeff_pf)
        
        
        !optimization of the Fourier acceleration parameters.
        call Fourier_transform_xmat(xmat,xmat_mom,myrank,x2p)
        call Fourier_acceleration_optimize(xmat_mom,fluctuation,myrank,imeasure,1)
        
     end if
     if(isave.EQ.0)then
        if(MOD(itraj,nsave).EQ.0)then
           call Save_Intermediate_Config(xmat,alpha,itraj)
        end if
     end if
     itraj=itraj+1
     
  end do
  !************************************************
  !************************************************
  !***  End of the loop of molecular evolutions ***
  !************************************************
  !************************************************    
  
  !close the output file for the measurement results. 
  if(myrank.eq.0)then
     close(unit_measurement)
     close(unit_CG_log)
     close(unit_Polyakov_phase)
     if(isave.EQ.0)then
        close(unit_intermediate_config)
     end if
  end if
  !*************************************
  !*** Save the final configuration ****
  !*************************************
 
  call Save_Final_Config(xmat,alpha,itraj,output_config)
  call Save_Acceleration_Parameters(fluctuation,imeasure,acc_output)

  call MPI_FINALIZE(IERR)
  !*******************************************
  !*** It's the end, have a good weekend! ****
  !*******************************************

end program BFSS
