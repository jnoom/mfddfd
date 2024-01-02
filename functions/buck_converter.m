function [T,Y] = buck_converter(V,L,C,Rn,fs,fsr,tf,Rv,duty,x0)
% state-space matrices 
B = [1/L;0];
c = [0,1];
D = 0;

Rv2   = imresize(Rv, [1 fs*tf], 'nearest');     % Upsampling of variable load resistance
Rn2   = imresize(Rn, [1 fs*tf], 'nearest');     % Upsampling of nominal load resistance
R     = Rn2+Rv2;                                % Total load resistance
duty2 = imresize(duty, [1 fs*tf], 'nearest');   % Upsampling of variable duty ratio
y     = zeros(length(duty2),2);                 % Output voltage (initialization)

warning off
t = 0;
for i=1:fs*tf
    A         = [0, -1/L; 1/C, -1/(R(i)*C)];     % Variable state transition matrix 
    sysc      = ss(A,B,c,D);
    t         = [t, ((i-1)+duty2(i))/fs, i/fs];  % Calculate time instances based on duty ratio
    [y2,~,x]  = lsim(sysc,[V V],t(end-2:end-1),x0(end,:));  % ON-mode
    y(i,1)    = y2(end);                         % Final output after ON-mode
    [y2,~,x0] = lsim(sysc,[0 0],t(end-1:end),x(end,:));     % OFF-mode
    y(i,2)    = y2(end);                         % Final output after OFF-mode
end
warning on

T = t(3:2*fs/fsr:end);
Y = y(1:fs/fsr:end,2)';
end

