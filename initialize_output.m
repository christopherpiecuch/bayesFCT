%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialize_output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code last edited by CGP on 25 November 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A = nan(NN_thin,M);         % Spatial field of error trend
B = nan(NN_thin,N);         % Spatial field of process linear trend
L = nan(NN_thin,M);         % Spatial field of observational biases
R = nan(NN_thin,1);         % AR(1) coefficient of the process
Y = nan(NN_thin,K,N);       % Process values
Y_0 = nan(NN_thin,N);       % Process initial conditions
MU = nan(NN_thin,1);        % Mean value of process linear trend
NU = nan(NN_thin,1);        % Mean value of observational biases
PHI = nan(NN_thin,1);       % Inverse range of process innovations
LAMBDA = nan(NN_thin,1);    % Inverse range of process trend field
PI_2 = nan(NN_thin,1);      % Spatial variance of process linear trend
SIGMA_2 = nan(NN_thin,1);   % Sill of the process innovations
DELTA_2 = nan(NN_thin,1);   % Instrumental error variance 
TAU_2 = nan(NN_thin,1);     % Spatial variance in observational biases
GAMMA_2 = nan(NN_thin,1);   % Spatial variance in error trends

RHO = nan(NN_thin,1);       % Inverse range of process trend field
ALPHA = nan(NN_thin,1);     % Inverse range of process trend field
OMEGA_2 = nan(NN_thin,1);   % Inverse range of process trend field
TR_0 = nan(NN_thin,1);      % Inverse range of process trend field
TR = nan(NN_thin,K);        % Process values
