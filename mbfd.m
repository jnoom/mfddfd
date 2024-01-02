clear all
addpath('functions')
rng(2123)                   % seed for random variables 
seed0 = randi([0 1],1,32);  % random seed for PRBS

%% Variables for model-based fault diagnosis algorithm
s2        = 8;              % VARX order for model-based fault diagnosis
lambda_mb = 0.3;            % tuning parameter for 1-norm for model-based fault diagnosis
m_error   = 0.4;            % model error for erroneous model
obs_poles = [0.2, 0.21];    % observer poles

%% Input/output dimensions
nu  = 1;            % no. inputs
ny  = 1;            % no. outputs
nd  = 1;            % no. disturbances
tf  = 0.06;         % final simulation time [s]
fsr = 20*10^3;      % sampling rate [Hz] (fsr*tf should result in an integer)
N   = tf*fsr-1;     % number of samples

%% Dictionary
nth = 50;                        % defined number of possible faults
fd  = linspace(1000,1900,nth)';  % frequencies for variable load resistance [Hz]
th  = reshape(square(2*pi*(fd*(1:N+1))/fsr),[nd,nth,N+1]);  % dictionary

%% Model
VB    = 10;             % input voltage [V]
L     = 220*10^-6;      % inductance [H]
C     = 330*10^-6;      % capacitance [C]
Rn    = 5;              % nominal load resistance [Ohm]
Rv    = 1;              % magnitude variable load resistance [Ohm]
dutyn = 0.45;           % nominal duty ratio
dutyv = 0.05;           % magnitude of variable duty ratio
fs    = 200*10^3;       % switching frequency [Hz] (fs/fsr should result in an integer)
R     = 10^-3*eye(ny);  % Measurement noise > 0

[PRBS1,seed1] = prbs(32,N+1,seed0);
[PRBS2,~]     = prbs(32,N+1,seed1);

u  = dutyn + dutyv*(2*PRBS1-1);     % input for identification set
u2 = dutyn + dutyv*(2*PRBS2-1);     % input for validation set
x0 = [dutyn*VB/Rn,dutyn*VB];        % initial condition [inductor current, resistor voltage]

%% Fault signals
nz     = 3;                         % number of nonzero elements in z
zmag   = 1;                         % magnitude of z
z      = zeros(nth,1);              % faults
nzi    = randperm(nth,nz)';         % indices of nonzeros in z
z(nzi) = zmag;                      % magnitudes of nonzeros in z

for k=1:N+1
    d(:,k) = th(:,:,k)*z;           % fault signals
end

%% Output sequences
[t,yr]   = buck_converter(VB,L,C,Rn,fs,fsr,tf,d,u,x0);      % noise-free identification set
y        = yr + sqrt(R)*randn(ny,N+1);                      % measured output for identification set

[t2,yr2] = buck_converter(VB,L,C,Rn,fs,fsr,tf,d,u2,x0);     % noise-free validation set
y2       = yr2 + sqrt(R)*randn(ny,N+1);                     % measured output for validation set

%% Model-based fault diagnosis
% Averaged models of the buck converter
Ac    = [0, -1/L; 1/C, -1/(Rn*C)];
Bc    = [VB/L; 0];
Cc    = [0, 1];
Dc    = 0;
Fc    = [0; dutyn*VB/(C*Rn^2)];
sysd  = c2d(ss(Ac,Bc,Cc,Dc),1/fsr);
sysd2 = c2d(ss(Ac,Fc,Cc,Dc),1/fsr);
Ad=sysd.A; Bd=sysd.B; Cd=sysd.C; Dd=sysd.D; Fd=sysd2.B;
Ade   = Ad + [0, 0; m_error, 0];        % state transition matrix for erroneous model

% Design observers
K  = place(Ad',Cd',obs_poles)';   % observer gain for error-free model
Ke = place(Ade',Cd',obs_poles)';  % observer gain for erroneous model

% Construct VARX parameters
for j = 1:s2
    Bv((j-1)*nu+1:j*nu,1:ny)  = (Cd * (Ad-K*Cd)^(j-1) * Bd)';
    Fv((j-1)*nd+1:j*nd,1:ny)  = (Cd * (Ad-K*Cd)^(j-1) * Fd)';
    Kv((j-1)*ny+1:j*ny,1:ny)  = (Cd * (Ad-K*Cd)^(j-1) * K)';
    Bve((j-1)*nu+1:j*nu,1:ny) = (Cd * (Ade-Ke*Cd)^(j-1) * Bd)';
    Fve((j-1)*nd+1:j*nd,1:ny) = (Cd * (Ade-Ke*Cd)^(j-1) * Fd)';
    Kve((j-1)*ny+1:j*ny,1:ny) = (Cd * (Ade-Ke*Cd)^(j-1) * Ke)';
end

% Construct data matrices
for i=1:N+1
    uv2(i,:)  = {reshape(u(:,i),[1,length(u(:,1))])};
    yv2(i,:)  = {reshape(y(:,i),[1,length(y(:,1))])};
    thv2(i,:) = {reshape(th(:,:,i)',[1,numel(th(:,:,1))])};
end
Tu2  = cell2mat(uv2(toeplitz(s2:N,s2:-1:1)));
Ty2  = cell2mat(yv2(toeplitz(s2:N,s2:-1:1)));
Tth2 = cell2mat(thv2(toeplitz(s2:N,s2:-1:1)));

% Perform optimization
zh  = lasso_cvx(y(:,s2+1:N+1)',Tu2,Ty2,Tth2,Bv,Fv,Kv,lambda_mb,nth,ny);     % Estimated faults using error-free model
zhe = lasso_cvx(y(:,s2+1:N+1)',Tu2,Ty2,Tth2,Bve,Fve,Kve,lambda_mb,nth,ny);  % Estimated faults using erroneous model

%% Figure
figure;plot(fd,zh,'o','Color','#D95319');hold on;
plot(fd,zhe,'d','Color','#EDB120')
plot([fd(nzi(1)) fd(nzi(1))],[-0.2, 1.2],'--','Color','#0072BD')
plot([fd(nzi(2)) fd(nzi(2))],[-0.2, 1.2],'--','Color','#0072BD')
plot([fd(nzi(3)) fd(nzi(3))],[-0.2, 1.2],'--','Color','#0072BD');
xlim([min(fd) max(fd)]);ylim([-0.2, 1.2])
xlabel('Square wave frequency [Hz]','Interpreter','latex')
ylabel('Normalized value [-]','Interpreter','latex')
legend('$\hat{z}$ using error-free model','$\hat{z}$ using erroneous model','Nonzero entries of $z$','Interpreter','latex')
