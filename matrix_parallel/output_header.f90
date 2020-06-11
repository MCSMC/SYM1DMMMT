  !*************************************************************************
  !******  Open the output file (for measurements) and set the header ******
  !*************************************************************************
SUBROUTINE output_header(myrank,data_output,temperature,flux,&
     ntau,dtau_xmat,dtaU_alpha,neig_max,neig_min,nbc,nbmn,&
     &init,input_config,output_config,iaccelerate,acc_input,acc_output,&
     &g_alpha,g_R,RCUT,upper_approx,max_err,max_iteration,CG_log,Pol_phase,&
     &isave,nsave,intermediate_config,imetropolis,ngauge,purebosonic,&
     &polfix,g_pol,pol_fix_width,myersfix,g_myers,myers_fix_width)

  implicit none
  include 'size_parallel.h'
  include '../unit_number.inc'
  integer isave,nsave,ngauge,purebosonic
  character(1000) data_output,input_config,output_config,&
       &acc_input,acc_output,intermediate_config,CG_log,Pol_phase
  character(100)input_config_short,output_config_short,&
       &acc_input_short,acc_output_short,intermediate_config_short,CG_log_short
  double precision temperature
  double precision flux
  double precision g_alpha,g_R,RCUT,upper_approx,max_err
  integer max_iteration
  integer neig_max,neig_min,nbmn,nbc,ntau,myrank,init,iaccelerate,imetropolis
  double precision dtau_xmat,dtau_alpha
  double precision polfix,g_pol,pol_fix_width,myersfix,g_myers,myers_fix_width

  !************************************************************
  !*** They are too long to write in the output file... :' ****
  !************************************************************
  input_config_short=input_config
  output_config_short=output_config
  acc_input_short=acc_input
  acc_output_short=acc_output
  intermediate_config_short=intermediate_config
  CG_log_short=CG_log
  !*******************
  !*** data_output ***
  !*******************
  if(myrank.eq.0)then
     open(unit=unit_measurement,status='REPLACE',file=data_output,&
          &action='WRITE')
     if(nbc.EQ.0)then
        write(unit_measurement,*)"#pbc"
     else if(nbc.EQ.1)then
        write(unit_measurement,*)"#apbc"
     end if
     if(nbmn.EQ.0)then
        write(unit_measurement,*)"#BFSS matrix model"
     else if(nbmn.EQ.1)then
        write(unit_measurement,*)"#BMN matrix model"
     end if
     if(ngauge.eq.0) then
        write(unit_measurement,*)"#gauged theory"
     else
        write(unit_measurement,*)"#ungauged theory"
     end if
     if(purebosonic.eq.1) then
        write(unit_measurement,*)"#pure bosonic theory"
     else
        write(unit_measurement,*)"#full theory"
     end if
     if(nimprove.EQ.1)then
        write(unit_measurement,*) "#improved lattice action"
     else if(nimprove.EQ.0)then
        write(unit_measurement,*) "#naive lattice action"
     end if
     write(unit_measurement,*) "#size of the gauge group: nmat=",&
          nmat_block*nblock
     write(unit_measurement,*) "#size of the lattice: Nt=",&
          nsite_local*nsublat
     write(unit_measurement,*) "#number of nodes: nsublat*nblock*nblock=",&
          &nsublat*nblock*nblock
     write(unit_measurement,*) "#number of processes along time: nsublat=",&
          nsublat
     write(unit_measurement,*) "#number of processes for matrix components: &
          &nblock*nblock=",nblock*nblock
     write(unit_measurement,*) "#number of MPI processes in this run : &
          &nsublat*(nblock)^2=",nsublat*nblock*nblock
     write(unit_measurement,*) "#Temperature=",temperature
     if(nbmn.EQ.1)then
        write(unit_measurement,*)"#flux paramerer: mu=", flux
     end if
     write(unit_measurement,*) "#g_alpha=",g_alpha
     write(unit_measurement,*) "#g_R=",g_R
     write(unit_measurement,*) "#RCUT=",RCUT
     write(unit_measurement,*) "#(TrX^2 is constrained so that TrX^2 < RCUT)"
     write(unit_measurement,*) "#ntau=",ntau
     write(unit_measurement,*) "#dtau for xmat=",Dtau_xmat
     write(unit_measurement,*) "#dtau for alpha=",Dtau_alpha
     write(unit_measurement,*) "#CG stopping condition: relative error=",max_err
     write(unit_measurement,*) "#CG stopping condition: maximum CG iteration=",&
          max_iteration
     write(unit_measurement,*) "#upper limit of remez approximation &
          &(Largest eigenvalue of D^dag*D must be smaller than this):",&
          &upper_approx

     write(unit_measurement,*) "#init=",init
     write(unit_measurement,*) "#(init; 0 -> old config, 1 -> new config, &
          &2-> fuzzy sphere)"
     write(unit_measurement,*) "#imetropolis=",imetropolis
     write(unit_measurement,*) "#(imetropolis=1 -> no Metropolis test,&
          &just for thermalization)"
     if(init.EQ.0)then
          write(unit_measurement,*) "#input configuration:",input_config_short
     end if
     write(unit_measurement,*) "#output configuration:",output_config_short
     write(unit_measurement,*) "#iaccelerate=",iaccelerate
     write(unit_measurement,*) "#(iaccelerate; 0 -> read from acc_input, &
          &1-> naive)"
     if(iaccelerate.EQ.0)then
        write(unit_measurement,*) "#Fourier acceleration file &
             &used for this run:",&
             &acc_input_short
     end if
     write(unit_measurement,*) "#Fourier acceleration file obtained &
          &from this run:",acc_output_short
     if(isave.EQ.0)then
        write(unit_measurement,*) "#intermediate configurations are saved in:",&
             &intermediate_config_short
        write(unit_measurement,*) "#saved every nsave configurations, nsave=",&
             nsave
     end if
     write(unit_measurement,*) "#log of the CG solver is saved in:",&
          &CG_log_short
     write(unit_measurement,*) "# fixing Polyakov loop ",&
          & polfix," ",g_pol," ",pol_fix_width
     write(unit_measurement,*) "# fixing Myers term ",&
          & myersfix," ",g_myers," ",myers_fix_width

     !********************************************************
     !**** Explain the ordering of the measurement output ****
     !********************************************************
     if((neig_max.EQ.0).AND.(neig_min.EQ.0))then
        write(unit_measurement,*) "# traj,ham_fin-ham_init, &
             &number of constraint violation,&
             &number of CG-err,CG-iteration,energy, &
             &|Pol.|, trX^2(total; 1, 2, ..., 9), tr[X,Y]^2, cubic term, acceptance"    
     else  if((neig_max.GT.0).AND.(neig_min.EQ.0))then
        write(unit_measurement,*)"# iteration for evaluating largest eig=",&
             &neig_max
        write(unit_measurement,*) "# traj,ham_fin-ham_init, &
             &number of constraint violation,&
             &number of CG-err,CG-iteration,energy, &
             &|Pol.|, trX^2(total; 1, 2, ..., 9), tr[X,Y]^2, cubic term, largest eig, acceptance"    
     else  if((neig_max.EQ.0).AND.(neig_min.GT.0))then
        write(unit_measurement,*)"# iteration for evaluating smallest eig=",&
             &neig_min
        write(unit_measurement,*) "# traj,ham_fin-ham_init, &
             &number of constraint violation,&
             &number of CG-err,CG-iteration,energy, &
             &|Pol.|, trX^2(total; 1, 2, ..., 9), tr[X,Y]^2, cubic term, smallest eig, acceptance"  
        
     else  if((neig_max.GT.0).AND.(neig_min.GT.0))then
        write(unit_measurement,*)"# iteration for evaluating largest eig=",&
             &neig_max
        write(unit_measurement,*)"# iteration for evaluating smallest eig=",&
             &neig_min
        write(unit_measurement,*) "# traj,ham_fin-ham_init, &
             &number of constraint violation,&
             &number of CG-err,CG-iteration,energy, &
             &|Pol.|, trX^2(total; 1, 2, ..., 9), tr[X,Y]^2, cubic term, largest eig,&
             &smallest eig, acceptance"  
     end if
     write(unit_measurement,*)'#-----------------------------------------------'
  end if

  !**************
  !*** CG_log ***
  !**************
  if(myrank.eq.0)then
     open(unit=unit_CG_log,status='REPLACE',file=CG_log,&
          &action='WRITE')
  end if

  !***************************
  !*** Polyakov line phase ***
  !***************************
  if(myrank.eq.0)then
     open(unit=unit_Polyakov_phase,status='REPLACE',file=Pol_phase,&
          &action='WRITE')
  end if
  

  return
END SUBROUTINE output_header
