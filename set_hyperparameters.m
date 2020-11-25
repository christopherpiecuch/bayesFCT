%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set_hyperparameters
% This function assigns hyperparameter values in the prior distributions.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code last edited by CGP on 25 November 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function HP = set_hyperparameters(N,K,M,DATA,TRANS,INDEX)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% As preparation for setting hyperparameter values, perform ordinary least
% squares line fits to available data at each tide gauge, and then estimate
% the noise properties (autocorrelation, white noise variance) of the
% residuals. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t=1:K;
m=nan(M,1);
s=nan(M,1);
r=nan(M,1);
e=nan(M,1);
y0=nan(M,1);
l=nan(M,1);
for n=1:M
    y=[]; 
    y=squeeze(DATA(n,:)); 
    i=[]; 
    i=find(~isnan(y));
        p=[]; q=[]; [p,q]=polyfit(t(i),y(i),1);
        w=(inv(q.R)*inv(q.R)')*q.normr^2/q.df; 
        m(n)=p(1);
        s(n)=w(1,1);
        [a b]=aryule(y(i)-p(1)*t(i)-p(2),1);
        d(n)=std(y(i)-p(1)*t(i)-p(2));
        r(n)=-a(2);
        e(n)=sqrt(b);
        l(n)=p(1)*mean(t)+p(2);
        y0(n)=p(2)-l(n);
    clear y i p w a b
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set variance inflation parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
var_infl=5^2;
var_infl2=10^2;
var0_infl=1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define the hyperparameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% y_0

HP.eta_tilde_y_0         	= nanmean(y0); % Mean of y_0 prior
HP.delta_tilde_y_0_2        = 0.25*var0_infl*(nanvar(y0)); % Variance of y_0 prior

% r
HP.u_tilde_r                = 0;% Lower bound of r prior
HP.v_tilde_r                = 1; % Upper bound of r prior

% mu
HP.eta_tilde_mu             = nanmean(m);% Mean of mu prior
HP.delta_tilde_mu_2         = var_infl2*nanvar(m); % Variance of mu prior

% nu
HP.eta_tilde_nu             = nanmean(l); % Mean of nu prior
HP.delta_tilde_nu_2         = var_infl*nanvar(l); % Variance of nu prior

% pi_2
HP.lambda_tilde_pi_2        = 1/2; % Shape of pi_2 prior
HP.nu_tilde_pi_2            = 1/2*nanvar(m); % Inverse scale of pi_2 prior

% delta_2
HP.lambda_tilde_delta_2     = 1/2; % Shape of delta_2 prior
HP.nu_tilde_delta_2         = 1/2*1e-4; % Guess (1 cm)^2 error variance

% sigma_2
HP.lambda_tilde_sigma_2     = 1/2; % Shape of sigma_2 prior
HP.nu_tilde_sigma_2         = 1/2*nanmean(e.^2); % Inverse scale of sigma_2 prior

% tau_2
HP.lambda_tilde_tau_2       = 1/2; % Shape of tau_2 prior
HP.nu_tilde_tau_2           = 1/2*nanvar(l); % Inverse scale of tau_2 prior

% gamma_2
HP.lambda_tilde_gamma_2     = 1/2; % Shape of tau_2 prior
HP.nu_tilde_gamma_2         = 1/2*(1e-3)^2; % Guess (1 mm/yr)^2 error variance

% phi
HP.eta_tilde_phi            = -7; % "Mean" of phi prior
HP.delta_tilde_phi_2        = 5; % "Variance" of phi prior

% lambda (this one's strongly constrained; 95% within 500,2000 km)
HP.eta_tilde_lambda         = -6.9; % "Mean" of phi prior
HP.delta_tilde_lambda_2     = 0.1225; % "Variance" of phi prior


% tr_0
HP.eta_tilde_tr_0           = nanmean(TRANS);% Mean of tr_0 prior (time-mean transport)
HP.delta_tilde_tr_0_2       = var_infl*nanvar(TRANS); % Variance of tr_0 prior

% rho
HP.eta_tilde_rho            = 0; %nanstd(TRANS)/nanmean(d);% Mean of rho prior (transfer function Sv/m between SL and GS)
HP.delta_tilde_rho_2        = var_infl*(nanstd(TRANS)/nanmean(d))^2; % Variance of rho prior

% alpha
HP.eta_tilde_alpha          = 0;% Mean of alpha prior (error trend in GS)
HP.delta_tilde_alpha_2      = (0.25)^2; % Variance of alpha prior % Sv/yr

% omega_2
HP.lambda_tilde_omega_2     = 1/2; % Shape of omega_2 prior
HP.nu_tilde_omega_2         = 1/2*nanvar(TRANS); % Guess (1 mm/yr)^2 error variance
