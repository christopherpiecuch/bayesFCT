%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%v
% set_initial_values
% Define initial process and parameter values.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code last edited by CGP on 25 November 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% mean parameters
mu=normrnd(HP.eta_tilde_mu,sqrt(HP.delta_tilde_mu_2));
nu=normrnd(HP.eta_tilde_nu,sqrt(HP.delta_tilde_nu_2));
alpha=normrnd(HP.eta_tilde_alpha,sqrt(HP.delta_tilde_alpha_2));
rho=normrnd(HP.eta_tilde_rho,sqrt(HP.delta_tilde_rho_2));
tr_0=normrnd(HP.eta_tilde_tr_0,sqrt(HP.delta_tilde_tr_0_2));

% variance parameters
pi_2=min([1 1/randraw('gamma', [0,1/HP.nu_tilde_pi_2,HP.lambda_tilde_pi_2], [1,1])]); % use min to prevent needlessly large values
delta_2=min([1 1/randraw('gamma', [0,1/HP.nu_tilde_delta_2,HP.lambda_tilde_delta_2], [1,1])]); % use min to prevent needlessly large values
sigma_2=min([1 1/randraw('gamma', [0,1/HP.nu_tilde_sigma_2,HP.lambda_tilde_sigma_2], [1,1])]); % use min to prevent needlessly large values
tau_2=min([1 1/randraw('gamma', [0,1/HP.nu_tilde_tau_2,HP.lambda_tilde_tau_2], [1,1])]); % use min to prevent needlessly large values
gamma_2=min([1 1/randraw('gamma', [0,1/HP.nu_tilde_gamma_2,HP.lambda_tilde_gamma_2], [1,1])]); % use min to prevent needlessly large values
omega_2=min([1 1/randraw('gamma', [0,1/HP.nu_tilde_omega_2,HP.lambda_tilde_omega_2], [1,1])]); % use min to prevent needlessly large values


% inverse length scale parameters
phi=exp(normrnd(HP.eta_tilde_phi,sqrt(HP.delta_tilde_phi_2)));
lambda=exp(normrnd(HP.eta_tilde_lambda,sqrt(HP.delta_tilde_lambda_2)));

% spatial fields
a=(mvnrnd(zeros(M,1),gamma_2*eye(M)))';        
b=(mvnrnd(mu*ones(N,1),pi_2*exp(-lambda*D)))';        
l=(mvnrnd(nu*ones(M,1),tau_2*eye(M)))';        
        
% AR(1) parameter
r=HP.u_tilde_r+(HP.v_tilde_r-HP.u_tilde_r )*rand(1);

% process
y_0=zeros(N,1);
y=zeros(N,K);
tr=tr_0*ones(K,1);
