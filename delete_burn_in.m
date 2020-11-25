%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% delete_burn_in
% Delete warm-up iterations from solutions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code last edited by CGP on 25 November 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MU(1:NN_burn_thin)=[];
NU(1:NN_burn_thin)=[];
PI_2(1:NN_burn_thin)=[];
DELTA_2(1:NN_burn_thin)=[];
SIGMA_2(1:NN_burn_thin)=[];
TAU_2(1:NN_burn_thin)=[];
GAMMA_2(1:NN_burn_thin)=[];
LAMBDA(1:NN_burn_thin)=[];
PHI(1:NN_burn_thin)=[];
A(1:NN_burn_thin,:)=[];
B(1:NN_burn_thin,:)=[];
L(1:NN_burn_thin,:)=[];
R(1:NN_burn_thin)=[];
Y_0(1:NN_burn_thin,:)=[];
Y(1:NN_burn_thin,:,:)=[];

TR_0(1:NN_burn_thin)=[];
ALPHA(1:NN_burn_thin)=[];
RHO(1:NN_burn_thin)=[];
OMEGA_2(1:NN_burn_thin)=[];
TR(1:NN_burn_thin,:)=[];
