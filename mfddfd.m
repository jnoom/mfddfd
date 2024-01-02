clear all
addpath('functions')
rng(2123)                   % seed for random variables 
seed0 = randi([0 1],1,32);  % random seed for PRBS

%% Variables for model-free data-driven fault diagnosis algorithm
s         = 3;      % VARX order for model-free data-driven fault diagnosis
lambda1   = 10;     % tuning parameter for nuclear norm 
lambda2   = 10;     % tuning parameter for 1-norm
gamma     = 10^-3;  % step size of PPXA algorithm
Nmax      = 10^4;   % maximum number of iterations of PPXA algorithm
eps       = 10^-1;  % threshold for fault diagnosis

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

%% Model-free data-driven fault diagnosis
% Construct data matrices
for i=1:N+1
    uv(i,:)  = {reshape(u(:,i),[1,length(u(:,1))])};
    yv(i,:)  = {reshape(y(:,i),[1,length(y(:,1))])};
    thv(i,:) = {reshape(th(:,:,i)',[1,numel(th(:,:,1))])};
end
Tu  = cell2mat(uv(toeplitz(s:N,s:-1:1)));
Ty  = cell2mat(yv(toeplitz(s:N,s:-1:1)));
Tth = cell2mat(thv(toeplitz(s:N,s:-1:1)));
H   = kron(eye(ny),[Tu,Tth,Ty]);
Y   = reshape(y(:,s+1:N+1)',[],1);

% Perform optimization
X0     = zeros(s*(nu+nd*nth)+s*ny,ny); 
X_ppxa = ppxa(Y,H,vec(X0),0.5*lambda1,0.5*lambda2,gamma,s,nu,ny,nd,nth,Nmax);

% Extract model parameters B and K
Bh = X_ppxa(1:s*nu,1:ny);
Kh = X_ppxa(s*(nu+nd*nth)+1:s*(nu+nd*nth)+s*ny,1:ny);

% Extract model parameter F and faults z
Fz = X_ppxa(s*nu+1:s*(nu+nd*nth),1:ny);
for j=1:nd
    for i=1:s
        F2z(ny*(j-1)+(i-1)*ny*nd+1:ny*j+(i-1)*ny*nd,:) = ...
            Fz(nth*(j-1)+(i-1)*nth*nd+1:nth*j+(i-1)*nth*nd,:)';
    end
end
[U,S,V] = svd(F2z);

for k=1:N+1
    dh(:,k) = th(:,:,k)*V(:,1);
end
Fh = reshape(U(:,1)*S(1,1),ny,nd*s)';

ind = abs(V(:,1))>max(abs(V(:,1)))*eps;
zh1 = ind.*V(:,1);

%% Second step of MFDDFD disregarding the 1-norm
thr  = th(:,ind,:);         % reduced size of theta
nth2 = length(thr(1,:,1));  % no. diagnosed faults

% Update data matrix
for i=1:N+1
    th2v(i,:) = {reshape(thr(:,:,i)',[1,numel(thr(:,:,1))])};
end
Tth2 = cell2mat(th2v(toeplitz(s:N,s:-1:1)));
H2   = kron(eye(ny),[Tu,Tth2,Ty]); 

% Perform optimization without regularization on the 1-norm
X02 = zeros(s*(nu+nd*nth2)+s*ny,ny);  % Initial condition for optimization
X2  = ppxa(Y,H2,vec(X02),0.5*lambda1,0,gamma,s,nu,ny,nd,nth2,Nmax);

% Extract updated model parameters B and K
Bh2 = X2(1:s*nu,1:ny);
Kh2 = X2(s*(nu+nd*nth2)+1:s*(nu+nd*nth2)+s*ny,1:ny);

% Extract updated model parameter F and faults z
Fz2 = X2(s*nu+1:s*(nu+nd*nth2),1:ny);
for j=1:nd
    for i=1:s
        F2z2(ny*(j-1)+(i-1)*ny*nd+1:ny*j+(i-1)*ny*nd,:) = ...
          Fz2(nth2*(j-1)+(i-1)*nth2*nd+1:nth2*j+(i-1)*nth2*nd,:)';
    end
end
[U2,S2,V2] = svd(F2z2);

for k=1:N+1
    dhr(:,k)  = thr(:,:,k)*V2(:,1);
end
Fh2 = reshape(U2(:,1)*S2(1,1),ny,nd*s)';


%% Validation
for k = s+1:N+1
    yh2(:,k) = Bh2'*vec(u2(:,k-1:-1:k-s)) + Kh2'*vec(y2(:,k-1:-1:k-s)) + Fh2'*vec(dhr(:,k-1:-1:k-s));
end

VAF = vaf(yr2(:,s+1:N+1)',yh2(:,s+1:N+1)');  % Variance accounted for

true  = sort(nzi);  % real faults
diagn = find(ind);  % diagnosed faults

if length(find(ind)) == nz
    acc = min(sort(nzi) ==find(ind));  % Diagnosis performance
else
    acc = 0;
end

%% Figure
figure;
subplot(4,1,1);plot(t(s+1:N+1),yr2(:,s+1:N+1)');hold on;plot(t(s+1:N+1),yh2(:,s+1:N+1)');
xlim([0 0.02]);ylabel('$V_R(k)$ [V]','Interpreter','latex');
legend('noise-free system output','identified model output','Interpreter','latex')
subplot(4,1,2);plot(t(s+1:N+1),u2(s+1:N+1));ylim([0.3,0.6]);xlim([0 0.02])
ylabel('Duty ratio of $S_B(k)$ [-]','Interpreter','latex');
subplot(4,1,3);plot(t(s+1:N+1),d(s+1:N+1)+Rn);ylim([0,10]);xlim([0 0.02])
ylabel('Resistance $R_B(k)$ [$\Omega$]','Interpreter','latex');
xlabel('Time $t(k)$ [s]','Interpreter','latex')

subplot(4,1,4);plot(fd,V(:,1),'o','Color','#D95319');hold on;
plot([fd(nzi(1)) fd(nzi(1))],[-0.1, 0.7],'--','Color','#0072BD')
plot([fd(nzi(2)) fd(nzi(2))],[-0.1, 0.7],'--','Color','#0072BD')
plot([fd(nzi(3)) fd(nzi(3))],[-0.1, 0.7],'--','Color','#0072BD');xlim([min(fd) max(fd)])
xlabel('Square wave frequency [Hz]','Interpreter','latex')
ylabel('Normalized value [-]','Interpreter','latex')
legend('$\hat{z}$ from SVD','Nonzero entries of $z$','Interpreter','latex')
