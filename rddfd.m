clear all
addpath('functions')
rng(2024)                   % seed for random variables 
seed0 = randi([0 1],1,32);  % random seed for PRBS

%% Variables for recursive model-free data-driven fault diagnosis algorithm
s         = 3;      % VARX order for model-free data-driven fault diagnosis
lambda1   = 10;     % tuning parameter for nuclear norm 
lambda2   = 10;     % tuning parameter for 1-norm
fw        = 10^3;   % finite window length
gamma     = 10^-3;  % step size of PPXA algorithm
Nmax      = 10;     % maximum number of iterations of PPXA algorithm
eps       = 10^-1;  % threshold for fault diagnosis

%% Input/output dimensions
nu  = 1;            % no. inputs
ny  = 1;            % no. outputs
nd  = 1;            % no. disturbances
tf  = 0.3;          % final simulation time [s]
tzc = 0.2;          % fraction of simulation time when model/faults change
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
Rn2   = 3.5;            % nominal load resistance after model change [Ohm]
Rv    = 1;              % magnitude variable load resistance [Ohm]
dutyn = 0.45;           % nominal duty ratio
dutyv = 0.05;           % magnitude of variable duty ratio
fs    = 200*10^3;       % switching frequency [Hz] (fs/fsr should result in an integer)
R     = 10^-3*eye(ny);  % Measurement noise > 0

[PRBS1,seed1] = prbs(32,N+1,seed0);
[PRBS2,~]     = prbs(32,N+1,seed1);

u   = dutyn + dutyv*(2*PRBS1-1);     % input for identification set
u2  = dutyn + dutyv*(2*PRBS2-1);     % input for validation set
x0  = [dutyn*VB/Rn,dutyn*VB];        % initial condition [inductor current, resistor voltage]
Rnv = [Rn*ones(1,2*tzc*(N+1)), Rn2*ones(1,tzc*(N+1)), Rn*ones(1,round((1-3*tzc)*(N+1)))];  % Resulting nominal load resistance

%% Fault signals
nz       = 3;                       % number of nonzero elements in z
zmag     = 1;                       % magnitude of z
z        = zeros(nth,1);            % faults (first fault set)
z2       = zeros(nth,1);            % faults (second fault set)
nzi      = randperm(nth,nz)';       % indices of nonzeros in z (first fault set)
nzi2     = randperm(nth,nz)';       % indices of nonzeros in z2 (second fault set)
z(nzi)   = zmag;                    % magnitudes of nonzeros in z
z2(nzi2) = zmag;                    % magnitudes of nonzeros in z2

for k=1:N+1
    if k <= round(tzc*(N+1))
        d(:,k) = th(:,:,k)*z;       % fault signals (first fault set)
    elseif k <= round(3*tzc*(N+1))
        d(:,k) = th(:,:,k)*z2;      % fault signals (second fault set) 
    else
        d(:,k) = th(:,:,k)*z;       % fault signals (first fault set)
    end
end

%% Output sequences
[t,yr]   = buck_converter(VB,L,C,Rnv,fs,fsr,tf,d,u,x0);      % noise-free identification set
y        = yr + sqrt(R)*randn(ny,N+1);                       % measured output for identification set

[t2,yr2] = buck_converter(VB,L,C,Rnv,fs,fsr,tf,d,u2,x0);     % noise-free validation set
y2       = yr2 + sqrt(R)*randn(ny,N+1);                      % measured output for validation set

%% Online execution
Pk     = eye(s*(nu+nd*nth+ny));         % initial condition of covariance matrix
rk     = zeros(s*(nu+nd*nth+ny)*ny,1);  % initial condition of rk
w(:,s) = zeros(s*(nu+nd*nth+ny)*ny,1);  % initial condition of PPXA algorithm

for k = s+1:N+1
    tic
    
    % Update data matrices
    thnk = [];
    for i=k-1:-1:k-s
        thnk = [thnk;reshape(th(:,:,i)',[1,numel(th(:,:,1))])'];
    end
    Hk   = [reshape(u(:,k-1:-1:k-s),[],1);thnk;reshape(y(:,k-1:-1:k-s),[],1)]';   % new measurements
    HkPk = Hk*Pk;
    if k < s+1+fw
        Pk = Pk - (3*gamma*HkPk'/(eye(ny)+3*gamma*(HkPk*Hk'))) * HkPk;            % update Pk for adding new measurements
        for iy = 1:ny 
            rk((iy-1)*s*(nu+nd*nth+ny)+1:iy*s*(nu+nd*nth+ny)) = ...
                rk((iy-1)*s*(nu+nd*nth+ny)+1:iy*s*(nu+nd*nth+ny)) + Hk'*y(iy,k);  % update rk for adding new measurements
        end    
    else
        thnl = [];
        for i=k-fw-1:-1:k-fw-s
            thnl = [thnl;reshape(th(:,:,i)',[1,numel(th(:,:,1))])'];
        end
        Hl    = [reshape(u(:,k-fw-1:-1:k-fw-s),[],1);thnl;reshape(y(:,k-fw-1:-1:k-fw-s),[],1)]';   % old measurements
        Pk1   = Pk - (3*gamma*HkPk'/(eye(ny)+3*gamma*(HkPk*Hk'))) * HkPk;         % update Pk for adding new measurements
        HlPk1 = Hl*Pk1;
        Pk    = Pk1 - (3*gamma*HlPk1'/(-eye(ny)+3*gamma*(HlPk1*Hl'))) * HlPk1;    % update Pk for removing old measurements
        for iy = 1:ny
            rk((iy-1)*s*(nu+nd*nth+ny)+1:iy*s*(nu+nd*nth+ny)) = ...
                rk((iy-1)*s*(nu+nd*nth+ny)+1:iy*s*(nu+nd*nth+ny)) + Hk'*y(iy,k) - Hl'*y(iy,k-fw);  % update rk for adding new measurements and removing old measurements
        end    
    end
    
    % Run PPXA algorithm
    w(:,k) = ppxa_fw(Pk,rk,w(:,k-1),0.5*lambda1,0.5*lambda2,gamma,s,nu,ny,nd,nth,Nmax);
    t_calc(k) = toc;
    
    % Extract model parameters B and K
    w2 = reshape(w(:,k),[],ny);
    Bh1 = w2(1:s*nu,1:ny);
    Kh1 = w2(s*(nu+nd*nth)+1:s*(nu+nd*nth)+s*ny,1:ny);
    
    % Extract model parameter F and faults z
    Fz = w2(s*nu+1:s*(nu+nd*nth),1:ny);
    for j=1:nd
        for i=1:s
            F2z(ny*(j-1)+(i-1)*ny*nd+1:ny*j+(i-1)*ny*nd,:) = Fz(nth*(j-1)+(i-1)*nth*nd+1:nth*j+(i-1)*nth*nd,:)';
        end
    end
    [U,S,V] = svd(F2z);
    zh(:,k) = V(:,1);
    for k2=1:N+1
        dh(:,k2)  = th(:,:,k2)*zh(:,k);
    end
    Fh1 = reshape(U(:,1)*S(1,1),ny,nd*s)';
    
    % Validation
    if k <=round(tzc*(N+1))
        yrv = yr2(1:round(tzc*(N+1)));
        yv = y2(1:round(tzc*(N+1)));
        uv = u2(1:round(tzc*(N+1)));
        dhv  = dh(:,1:round(tzc*(N+1)));
    elseif k <=round(2*tzc*(N+1))
        yrv = yr2(round(tzc*(N+1))+1:round(2*tzc*(N+1)));
        yv = y2(round(tzc*(N+1))+1:round(2*tzc*(N+1)));
        uv = u2(round(tzc*(N+1))+1:round(2*tzc*(N+1)));
        dhv  = dh(:,round(tzc*(N+1))+1:round(2*tzc*(N+1)));
    elseif k <=round(3*tzc*(N+1))
        yrv = yr2(round(2*tzc*(N+1))+1:round(3*tzc*(N+1)));
        yv = y2(round(2*tzc*(N+1))+1:round(3*tzc*(N+1)));
        uv = u2(round(2*tzc*(N+1))+1:round(3*tzc*(N+1)));
        dhv  = dh(:,round(2*tzc*(N+1))+1:round(3*tzc*(N+1)));
    else 
        yrv = yr2(round(3*tzc*(N+1))+1:round(4*tzc*(N+1)));
        yv = y2(round(3*tzc*(N+1))+1:round(4*tzc*(N+1)));
        uv = u2(round(3*tzc*(N+1))+1:round(4*tzc*(N+1)));
        dhv  = dh(:,round(3*tzc*(N+1))+1:round(4*tzc*(N+1)));
    end
    
    for k2 = s+1:round(tzc*(N+1))
        yh1(:,k2) = Bh1'*vec(uv(:,k2-1:-1:k2-s)) + Kh1'*vec(yv(:,k2-1:-1:k2-s)) + Fh1'*vec(dhv(:,k2-1:-1:k2-s));
    end
    VAF_y1(:,k) = vaf(yrv(:,s+1:k2)',yh1(:,s+1:k2)');  % variance accounted for
end


%% Figure
figure
subplot(3,1,1);imagesc([(s+1)/(N+1)*tf, tf],[min(fd) max(fd)],abs(zh(:,s+1:end)),[0, 0.9]);
c = colorbar('location','eastoutside');
c.Label.Interpreter = 'latex';
c.Label.String = 'Normalized magnitude of $\hat{z}$';
colormap('turbo');
title('Regularized','Interpreter','latex')

hold on;
line([tzc*tf,tzc*tf], [min(fd) max(fd)]', 'Color', 'm', 'LineStyle', '--','LineWidth',1);
line([2*tzc*tf,2*tzc*tf], [min(fd) max(fd)]', 'Color', 'm', 'LineStyle', '--','LineWidth',1);
line([3*tzc*tf,3*tzc*tf], [min(fd) max(fd)]', 'Color', 'm', 'LineStyle', '--','LineWidth',1);
ylabel('Square wave freq. [Hz]','Interpreter','latex')
xlim([0 tf])

subplot(3,1,2);hold on;plot(t(s+1:end),[sqrt(mean(zh(setdiff(1:end,nzi),s+1:round(tzc*(N+1))).^2))'; sqrt(mean(zh(setdiff(1:end,nzi2),round(tzc*(N+1))+1:3*round(tzc*(N+1))).^2))'; sqrt(mean(zh(setdiff(1:end,nzi),3*round(tzc*(N+1))+1:end).^2))'],'Color','#D95319')
line([tzc*tf,tzc*tf], [0 0.2]', 'Color', 'k', 'LineStyle', '--','LineWidth',1);
line([2*tzc*tf,2*tzc*tf], [0 0.2]', 'Color', 'k', 'LineStyle', '--','LineWidth',1);
line([3*tzc*tf,3*tzc*tf], [0 0.2]', 'Color', 'k', 'LineStyle', '--','LineWidth',1);
ylim([0,0.15]);grid minor
ylabel('RMSE inactive freq. [-]','Interpreter','latex')
legend('Regularized','Interpreter','latex')

subplot(3,1,3);hold on;plot(t(s+1:end),VAF_y1(s+1:end)','Color','#D95319');ylim([70,100])
ylabel('VAF (\%)','Interpreter','latex');grid minor
line([tzc*tf,tzc*tf], [0 100]', 'Color', 'k', 'LineStyle', '--','LineWidth',1);
line([2*tzc*tf,2*tzc*tf], [0 100]', 'Color', 'k', 'LineStyle', '--','LineWidth',1);
line([3*tzc*tf,3*tzc*tf], [0 100]', 'Color', 'k', 'LineStyle', '--','LineWidth',1);
xlabel('Time $t(k)$ [s]','Interpreter','latex')
legend('Regularized','Interpreter','latex')