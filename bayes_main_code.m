%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function bayes_main_code(save_tag,save_num)
%   DESCRIPTION: Main driver code. Input are:
%   * save_tag: prefix to be appended to file name (see bottom).
%   * save_num: iteration number to allow multiple chains (see bottom).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code last edited by CGP on 25 November 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function bayes_main_code(save_tag,save_num)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialize some variables (Matlab can give you grief otherwise)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mu=[]; alpha=[]; years=[];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% define save-out file name
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
save_tag=[save_tag,num2str(save_num)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define iteration parameters based on input
% (Values are defaults from the paper)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% parameters related to iterations
%multer = 1; % for testing purposes
multer = 100; % as in paper
NN_burn = 1000*multer;  % Number of burn-in iterations
NN_post = 1000*multer;  % Number of post-burn-in iterations
thin_period = 1*multer; % Number of iterations to thin by
NN_burn_thin=NN_burn/thin_period;    
NN_post_thin=NN_post/thin_period;    
NN=NN_burn+NN_post;                 
NN_thin=NN_burn_thin+NN_post_thin;   

% parameters related to space and time
la1=5;      % Southern latitudinal bounds of study region
la2=35;     % Northern latitudinal bounds "
lo1=-100;   % Western longitudinal bounds "
lo2=-60;    % Eastern longitudinal bounds "
minnum=10;  % Minimum number of data points to consider a record
coastcode=[902 904 906 908 912 916 930 931 ...
    932 934 936 938 939 940 950 941 960]; % PSMSL ID for coastlines
USLonCrit=-81.81;   % Westernmost latitude of US tide gauges
USLatCrit=34.25;    % Northernmost latitude of US tide gauges
t0=1908;            % start year
tf=2018;            % end year

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Prepare and load PSMSL annual tide gauge data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load PSMSL data
[DATA,LON,LAT,NAME,COAST,TIME]=prepare_data(la1,la2,lo1,lo2,minnum,...
    coastcode,t0,tf);
% For US gauges, only retain those records within the "southeast" region 
% of Thompson and Mitchum (2014) https://doi.org/10.1002/2014JC009999
ii=[]; ii=find(((LAT>USLatCrit)|(LON<USLonCrit))&((COAST==940)|(COAST==960)));
DATA(ii,:)=[];
LON(ii)=[];
LAT(ii)=[];
NAME(ii)=[];
COAST(ii)=[];

% In the PSMSL dataset, the Virginia Key tide gauge has a bogus datum
% offset of ~6m, whereas all other gauges have offsets of ~7m. Since this
% values is arbitrary and has no physical significance in this context,
% adjust manually so that the datum offset at Virginia Key is ~7m as all
% other gauges
for n=1:numel(NAME);
    II(n)=strcmp(NAME(n).name,'VIRGINIA KEY, FL');
end
II=find(II);
if ~isempty(II)
    DATA(II,:)=DATA(II,:)+1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define space and time parameters related to data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[N,K]=size(DATA);
D=EarthDistances([LON' LAT']); 
T=1:K; T=T-mean(T);
T0=T(1)-1;
M=sum(sum(~isnan(DATA'))~=0); % number of locations with at least one datum

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load Florida cable data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('20190228_annual_fs_transport_sv_1982_2018.mat')
gs_data=annuTranVal;
gs_erro=annuTranErr;
gs_year=years;
clear annuTran* years

gs_year(find(isnan(gs_data)))=[];
gs_erro(find(isnan(gs_data)))=[];
gs_data(find(isnan(gs_data)))=[];
for nn=1:numel(gs_data)
    gs_indx(nn)=find(TIME==gs_year(nn));    
end
J=numel(gs_data);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Determine which clusters the various target locations fall within
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ClusMult = determine_clusters(COAST);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set the seeds of the random number generators
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
randn('state', save_num*sum((1000+600)*clock))
rand('state', save_num*sum((1000+800)*clock))
%rng(1e6*now);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Allocate space for the sample arrays
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
initialize_output

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define the hyperparameter values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%set_hyperparameters
HP = set_hyperparameters(N,K,M,DATA,gs_data,gs_indx);

%%%%%%%%%%%%%%%%%%%%
% Set initial values
%%%%%%%%%%%%%%%%%%%%
set_initial_values

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set up selection matrices
%%%%%%%%%%%%%%%%%%%%%%%%%%%
H_master=double(~isnan(DATA));
M_k=sum(H_master);
for k=1:K
    gauges_with_data(k).indices=find(H_master(:,k)~=0);
    selection_matrix(k).H=zeros(M_k(k),N);
    selection_matrix(k).F=zeros(M_k(k),M);
    for m_k=1:M_k(k)
       selection_matrix(k).H(m_k,gauges_with_data(k).indices(m_k))=1;
       selection_matrix(k).F(m_k,gauges_with_data(k).indices(m_k))=1;
    end
    Z(k).z=squeeze(DATA(gauges_with_data(k).indices,k));
end
DEL=zeros(1,N);
for n=1:numel(NAME)
   if strcmp(NAME(n).name,'SETTLEMENT POINT') % eastern end point
       DEL(n)=1;
   elseif strcmp(NAME(n).name,'LAKE WORTH PIER') % western end point
       DEL(n)=-1;
   else
       DEL(n)=0;
   end
end
GMAT=zeros(J,K); % j is data k is process
for jj=1:J
    GMAT(jj,gs_indx(jj))=1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set up identity matrices and vectors of zeros or ones
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
I_N=eye(N);
I_M=eye(M);
I_K=eye(K);
I_J=eye(J);
ONE_N=ones(N,1);
ONE_M=ones(M,1);
ONE_K=ones(K,1);
ONE_J=ones(J,1);
ZERO_N=zeros(N,1);
ZERO_M=zeros(M,1);
ZERO_K=zeros(K,1);
ZERO_J=zeros(J,1);
for k=1:K
   I_MK(k).I=eye(M_k(k));
   ONE_MK(k).ONE=ones(M_k(k),1);
   ZERO_MK(k).ZERO=zeros(M_k(k),1);
end

E=ZERO_J;
for jj=1:J
   E(jj,jj)=gs_erro(jj)^2;
end
invE=inv(E);

pause(1)
%plot(TIME,tr), hold on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Loop through the Gibbs sampler with Metropolis step
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
for nn=1:NN, if mod(nn,25)==0, toc, disp([num2str(nn),' of ',num2str(NN),' iterations done.']), tic, end
    nn_thin=[]; nn_thin=ceil(nn/thin_period);
     
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Define matrices to save time
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    PiMat=pi_2*exp(-lambda*D); invPiMat=inv(PiMat);
    SigmaMat=sigma_2*(ClusMult.*exp(-phi*D)); invSigmaMat=inv(SigmaMat);
    LMat=exp(-lambda*D); invLMat=inv(LMat);
    SMat=ClusMult.*exp(-phi*D); invSMat=inv(SMat);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Sample from p(y_K|.)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    V_Y_K=[]; PSI_Y_K=[];
    V_Y_K=rho/omega_2*DEL'*(tr(K)-tr_0-alpha*T(K))+delta_2^(-1)*(selection_matrix(K).H'*(Z(K).z-selection_matrix(K).F*(l+a*T(K))))+...
    	invSigmaMat*(r*y(:,K-1)+(T(K)-r*T(K-1))*b);
    PSI_Y_K=((rho^2)/omega_2*DEL'*DEL+1/delta_2*selection_matrix(K).H'*selection_matrix(K).H+invSigmaMat)^(-1);
    y(:,K)=mvnrnd(PSI_Y_K*V_Y_K,PSI_Y_K)';
    clear V_Y_K PSI_Y_K   

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Sample from p(y_k|.)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for kk=(K-1):-1:1
    	V_Y_k=[]; PSI_Y_k=[];
      	if kk==1
        	V_Y_k=rho/omega_2*DEL'*(tr(kk)-tr_0-alpha*T(kk))+1/delta_2*(selection_matrix(1).H'*(Z(1).z-selection_matrix(1).F*(l+a*T(1))))+...
                invSigmaMat*(r*(y_0+y(:,2))+(1+r^2)*T(1)*b-r*(T0+T(2))*b);
        else
         	V_Y_k=rho/omega_2*DEL'*(tr(kk)-tr_0-alpha*T(kk))+1/delta_2*(selection_matrix(kk).H'*(Z(kk).z-selection_matrix(kk).F*(l+a*T(kk))))+...
            	invSigmaMat*(r*(y(:,kk-1)+y(:,kk+1))+(1+r^2)*T(kk)*b-r*(T(kk-1)+T(kk+1))*b);
        end
       	PSI_Y_k=inv((rho^2)/omega_2*DEL'*DEL+1/delta_2*selection_matrix(kk).H'*selection_matrix(kk).H+(1+r^2)*invSigmaMat);
       	y(:,kk)=mvnrnd(PSI_Y_k*V_Y_k,PSI_Y_k)';
      	clear V_Y_k PSI_Y_k 
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Sample from p(y_0|.)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    V_Y_0=[]; PSI_Y_0=[];
    V_Y_0=(HP.eta_tilde_y_0/HP.delta_tilde_y_0_2)*ONE_N+invSigmaMat*(r*y(:,1)-r*(T(1)-r*T0)*b);
    PSI_Y_0=inv(1/HP.delta_tilde_y_0_2*I_N+r^2*invSigmaMat);
    y_0=mvnrnd(PSI_Y_0*V_Y_0,PSI_Y_0)';
    clear V_Y_0 PSI_Y_0
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Sample from p(b|.)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   	V_B=[]; PSI_B=[]; SUM_K=ZERO_N;
    for kk=1:K
     	if kk==1
         	SUM_K=SUM_K+(T(1)-r*T0)*(y(:,1)-r*y_0);
        else
          	SUM_K=SUM_K+(T(kk)-r*T(kk-1))*(y(:,kk)-r*y(:,kk-1));
        end
   	end
    V_B=mu*invPiMat*ONE_N+invSigmaMat*SUM_K;
    PSI_B=inv(invPiMat+invSigmaMat*sum((T-r*[T0 T(1:K-1)]).^2));
    b=mvnrnd(PSI_B*V_B,PSI_B)';
    clear V_B PSI_B SUM_K

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Sample from p(a|.)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   	V_A=[]; PSI_A=[]; SUM_K1=ZERO_M; SUM_K2=zeros(M);
    for kk=1:K
    	SUM_K1=SUM_K1+T(kk)*selection_matrix(kk).F'*(Z(kk).z-...
         	selection_matrix(kk).H*y(:,kk)-selection_matrix(kk).F*l);
        SUM_K2=SUM_K2+T(kk)^2*(selection_matrix(kk).F'*selection_matrix(kk).F);
   	end
    V_A=(1/delta_2)*SUM_K1;
    PSI_A=inv(1/gamma_2*eye(M)+1/delta_2*SUM_K2);
    a=mvnrnd(PSI_A*V_A,PSI_A)';
    clear V_A PSI_A SUM_K*

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Sample from p(mu|.)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    V_MU=[]; PSI_MU=[];
    V_MU=HP.eta_tilde_mu/HP.delta_tilde_mu_2+ONE_N'*invPiMat*b;
   	PSI_MU=inv(1/HP.delta_tilde_mu_2+ONE_N'*invPiMat*ONE_N);
    mu=normrnd(PSI_MU*V_MU,sqrt(PSI_MU));
    clear V_MU PSI_MU  
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Sample from p(pi_2|.)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    inside1=[]; inside2=[];
    inside1=N/2;
    inside2=1/2*((b-mu*ONE_N)'*invLMat)*(b-mu*ONE_N);
    pi_2=1/randraw('gamma', [0,1/(HP.nu_tilde_pi_2+inside2),...
      	(HP.lambda_tilde_pi_2+inside1)], [1,1]);
   	clear inside*
      
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Sample from p(delta_2|.)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    SUM_K=0;
    for kk=1:K
      	xxx=[]; xxx=(Z(kk).z-selection_matrix(kk).H*y(:,kk)-selection_matrix(kk).F*(l+a*T(kk)));
      	SUM_K=SUM_K+xxx'*xxx;
   	end
    delta_2=1/randraw('gamma', [0,1/(HP.nu_tilde_delta_2+1/2*SUM_K),...
     	(HP.lambda_tilde_delta_2+1/2*sum(M_k))], [1,1]);    
    clear SUM_K
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Sample from p(r|.)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    V_R=0; PSI_R=0;
    for kk=1:K
     	if kk==1
         	V_R=V_R+((y_0-b*T0)')*invSigmaMat*(y(:,1)-b*T(1));
         	PSI_R=PSI_R+((y_0-b*T0)')*invSigmaMat*(y_0-b*T0);
        else
         	V_R=V_R+((y(:,kk-1)-b*T(kk-1))')*invSigmaMat*(y(:,kk)-b*T(kk));
          	PSI_R=PSI_R+((y(:,kk-1)-b*T(kk-1))')*invSigmaMat*(y(:,kk-1)-b*T(kk-1));
        end        
   	end
    PSI_R=inv(PSI_R);
    dummy=1;
    while dummy
      	sample=normrnd(PSI_R*V_R,sqrt(PSI_R));
      	if sample>HP.u_tilde_r&&sample<HP.v_tilde_r
         	r=sample;
          	dummy=0;
        end
    end
    clear V_R PSI_R dummy ctr

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Sample from p(sigma_2|.)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    SUM_K=0;
    for kk=1:K
     	if kk==1
         	DYKK=[];
          	DYKK=y(:,1)-r*y_0-(T(1)-r*T0)*b;
         	SUM_K=SUM_K+(DYKK')*invSMat*DYKK;           
        else
         	DYKK=[];
           	DYKK=y(:,kk)-r*y(:,kk-1)-(T(kk)-r*T(kk-1))*b;
         	SUM_K=SUM_K+(DYKK')*invSMat*DYKK;           
        end
    end
   	sigma_2=1/randraw('gamma', [0,1/(HP.nu_tilde_sigma_2+1/2*SUM_K),...
     	(HP.lambda_tilde_sigma_2+N*K/2)], [1,1]);
   	clear SUM_K DYKK
            
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Sample from p(phi|.)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Phi_now=log(phi);
    Phi_std=0.05;
    Phi_prp=normrnd(Phi_now,Phi_std);
    R_now=ClusMult.*exp(-exp(Phi_now)*D);
    R_prp=ClusMult.*exp(-exp(Phi_prp)*D);
    invR_now=inv(R_now);
    invR_prp=inv(R_prp);
    sumk_now=0;
    sumk_prp=0;
        
   	for kk=1:K
      	if kk==1
         	DYYK=y(:,1)-r*y_0-(T(1)-r*T0)*b;
          	sumk_now=sumk_now+(DYYK')*invR_now*DYYK;
          	sumk_prp=sumk_prp+(DYYK')*invR_prp*DYYK;
        else
         	DYYK=y(:,kk)-r*y(:,kk-1)-(T(kk)-r*T(kk-1))*b;
         	sumk_now=sumk_now+(DYYK')*invR_now*DYYK;
         	sumk_prp=sumk_prp+(DYYK')*invR_prp*DYYK;
        end
    end
        
 	ins_now=-1/(2*HP.delta_tilde_phi_2)*(Phi_now-HP.eta_tilde_phi)^2-1/(2*sigma_2)*sumk_now;
   	ins_prp=-1/(2*HP.delta_tilde_phi_2)*(Phi_prp-HP.eta_tilde_phi)^2-1/(2*sigma_2)*sumk_prp;
  	MetFrac=det(R_prp*invR_now)^(-K/2)*exp(ins_prp-ins_now);
   	success_rate=min(1,MetFrac);
   	if rand(1)<=success_rate
     	Phi_now=Phi_prp; 
    end
  	phi=exp(Phi_now);
  	clear Phi_now Phi_std Phi_prp mat_now mat_prp ins_* sumk MetFrac success_rate R_*

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Sample from p(lambda|.)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Lambda_now=log(lambda);
    Lambda_std=0.05;
    Lambda_prp=normrnd(Lambda_now,Lambda_std);
    R_now=exp(-exp(Lambda_now)*D);
    R_prp=exp(-exp(Lambda_prp)*D);
    invR_now=inv(R_now);
    invR_prp=inv(R_prp);
    sumk_now=0;
    sumk_prp=0;
        
    sumk_now=(b-mu*ONE_N)'*invR_now*(b-mu*ONE_N);
    sumk_prp=(b-mu*ONE_N)'*invR_prp*(b-mu*ONE_N);
        
 	ins_now=-1/(2*HP.delta_tilde_lambda_2)*(Lambda_now-HP.eta_tilde_lambda)^2-1/(2*pi_2)*sumk_now;
   	ins_prp=-1/(2*HP.delta_tilde_lambda_2)*(Lambda_prp-HP.eta_tilde_lambda)^2-1/(2*pi_2)*sumk_prp;
  	MetFrac=det(R_prp*invR_now)^(-1/2)*exp(ins_prp-ins_now);
   	success_rate=min(1,MetFrac);
   	if rand(1)<=success_rate
     	Lambda_now=Lambda_prp; 
    end
  	lambda=exp(Lambda_now);
  	clear Lambda_now Lambda_std Lambda_prp mat_now mat_prp ins_* sumk MetFrac success_rate R_*

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Sample from p(l|.)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    V_L=[]; PSI_L=[]; SUM_K1=ZERO_M; SUM_K2=zeros(M,M);
    for kk=1:K
     	SUM_K1=SUM_K1+(selection_matrix(kk).F')*(Z(kk).z-selection_matrix(kk).H*y(:,kk)-...
            selection_matrix(kk).F*a*T(kk));
        SUM_K2=SUM_K2+(selection_matrix(kk).F')*selection_matrix(kk).F;
    end
    V_L=nu/tau_2*ONE_M+1/delta_2*SUM_K1;
    PSI_L=inv(1/tau_2*I_M+1/delta_2*SUM_K2);
    l=mvnrnd(PSI_L*V_L,PSI_L)';

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Sample from p(nu|.)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    V_NU=[]; PSI_NU=[];
   	V_NU=HP.eta_tilde_nu/HP.delta_tilde_nu_2+1/tau_2*(ONE_M'*l);
    PSI_NU=inv(1/HP.delta_tilde_nu_2+M/tau_2);
    nu=normrnd(PSI_NU*V_NU,sqrt(PSI_NU));
    clear V_NU PSI_NU
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Sample from p(tau_2|.)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    tau_2=1/randraw('gamma', [0,1/(HP.nu_tilde_tau_2+...
     	1/2*(((l-nu*ONE_M)')*(l-nu*ONE_M))),(HP.lambda_tilde_tau_2+M/2)], [1,1]);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Sample from p(gamma_2|.)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    gamma_2=1/randraw('gamma', [0,1/(HP.nu_tilde_gamma_2+...
     	1/2*a'*a),(HP.lambda_tilde_gamma_2+M/2)], [1,1]);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Sample from p(tr|.)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    V_TR=[]; PSI_TR=[];
    V_TR=GMAT'*invE*gs_data'+1/omega_2*(rho*y'*DEL'+alpha*T'+tr_0*ONE_K);
    %(HP.eta_tilde_y_0/HP.delta_tilde_y_0_2)*ONE_N+invSigmaMat*(r*y(:,1)-r*(T(1)-r*T0)*b);
    PSI_TR=inv(GMAT'*invE*GMAT+1/omega_2*I_K);
    tr=mvnrnd(PSI_TR*V_TR,PSI_TR)';
    clear V_TR PSI_TR
   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Sample from p(tr_0|.)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    V_TR_0=[]; PSI_TR_0=[];
   	V_TR_0=HP.eta_tilde_tr_0/HP.delta_tilde_tr_0_2+1/omega_2*ONE_K'*(tr-rho*y'*DEL'-alpha*T');
    PSI_TR_0=inv(1/HP.delta_tilde_tr_0_2+K/omega_2);
    tr_0=normrnd(PSI_TR_0*V_TR_0,sqrt(PSI_TR_0));
    clear V_TR_0 PSI_TR_0

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Sample from p(alha|.)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    V_ALPHA=[]; PSI_ALPHA=[];
   	V_ALPHA=HP.eta_tilde_alpha/HP.delta_tilde_alpha_2+1/omega_2*T*(tr-tr_0*ONE_K-rho*y'*DEL');
    PSI_ALPHA=inv(1/HP.delta_tilde_alpha_2+1/omega_2*(T*T'));
    alpha=normrnd(PSI_ALPHA*V_ALPHA,sqrt(PSI_ALPHA));
    clear V_ALPHA PSI_ALPHA

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Sample from p(rho|.)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    V_RHO=[]; PSI_RHO=[];
   	V_RHO=HP.eta_tilde_rho/HP.delta_tilde_rho_2+1/omega_2*DEL*y*(tr-tr_0*ONE_K-alpha*T');
    PSI_RHO=inv(1/HP.delta_tilde_rho_2+1/omega_2*DEL*y*y'*DEL');
    rho=normrnd(PSI_RHO*V_RHO,sqrt(PSI_RHO));
    clear V_RHO PSI_RHO

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Sample from p(omega_2|.)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    pr=[]; pr=tr-tr_0*ONE_K-alpha*T'-rho*y'*DEL';
    omega_2=1/randraw('gamma', [0,1/(HP.nu_tilde_omega_2+...
     	1/2*pr'*pr),(HP.lambda_tilde_omega_2+K/2)], [1,1]);
    clear pr  
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Now update arrays
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    update_all_arrays

end
toc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% delete the burn-in period values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
delete_burn_in

%plot(TIME,tr), hold on


%%%%%%%%%%%%%
% save output
%%%%%%%%%%%%%
save(['bayes_model_solutions/experiment_',save_tag,'_',num2str(t0),'_',num2str(tf),'.mat'],...
    'MU','NU','PI_2','DELTA_2','SIGMA_2','TAU_2','GAMMA_2',...
    'PHI','LAMBDA','A','B','L','R','Y_0','Y','HP','DATA','LON','LAT','COAST',...
    'NAME','N','K','D','TR','TR_0','ALPHA','RHO','OMEGA_2','TIME')