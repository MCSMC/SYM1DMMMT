
  integer nbc !boundary condition for fermions; 0 -> pbc, 1 -> apbc
  integer nbmn ! 0 -> BFSS, 1 -> BMN
  integer ngauge ! 0 -> gauged, 1 -> ungauged
  integer purebosonic ! 0 -> full theory, 1 -> pure bosonic part
  integer init !initial condition; 0 -> continue, 1 -> new
  integer isave !0 -> save intermediate config, 1 -> do not save
  integer nsave !saved every nsave trajectories 
!matrices
  double complex xmat(1:nmat_block,1:nmat_block,1:ndim,-nimprove:nsite_local+1+nimprove)
  double complex xmat_mom(1:nmat_block,1:nmat_block,1:ndim,1:nsite_local)

!remez coefficients
double precision acoeff_md(0:nremez_md),bcoeff_md(1:nremez_md)!molecular evolution
double precision acoeff_pf(0:nremez_pf),bcoeff_pf(1:nremez_pf)!pseudo fermion
double precision upper_approx!the largest eigenvalue of (M^¥dagger M)must be smaller than this.

!smearing
integer nsmear
double precision s_smear

!CG solver
integer max_iteration, iteration,n_bad_CG
double precision max_err
!For Mersenne Twister
integer mersenne_seed
!Fourier acceleration
double precision acceleration(1:nsite_local)
double precision fluctuation(1:nsite_local)
integer iaccelerate,imeasure
integer imetropolis

!number of CG iteration for calculating the largest and smallest eigenvalues of D=(M^dag*M)
integer neig_max,neig_min
!number of fuzzy sphere, when init=2.
integer nfuzzy
!Gamma matrices
double complex Gamma10d(1:ndim,1:nspin,1:nspin)




  double precision alpha(1:nmat_block*nblock)

  integer iremez
  integer ncv

double precision temperature
double precision flux
  !parameters for molecular evolution
  integer ntau
  double precision dtau_xmat,dtau_alpha

  !number of trajectories
  integer ntraj     !total number of trajectories at the end of the run
  integer itraj

  double precision ham_init,ham_fin


  !measurements
  integer nskip !measurement is performed every nskiptrajectories
  integer nacceptance !number of acceptance
  integer ntrial !number of trial
character(1000) input_config,data_output,output_config,acc_input,acc_output,intermediate_config,CG_log,Pol_phase


integer IERR,NPROCS,MYRANK



!coefficient for potential for alpha
double precision g_alpha
!coefficient for potential for Tr(x^2)
double precision g_R,RCUT
double precision polfix,g_pol,pol_fix_width,myersfix,g_myers,myers_fix_width

