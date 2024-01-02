function out = proxnn(x,gamma,tau,s,nu,ny,nd,nth)
out = x;  % keep original values for part of the vector at which the nuclear norm does not apply

for i=1:ny
    Fz(:,i) = x((i-1)*s*(nu+nth*nd+ny) + s*nu+1 : (i-1)*s*(nu+nth*nd+ny) + s*(nu+nth*nd));  % forward transormation 1
end
Fsz = reshape(Fz,nth,[]);                   % forward transormation 2

[U,S,V] = svd(Fsz,'econ');
S2      = diag(max(0,diag(S)-gamma*tau));   % soft-thresholding

Fszr = reshape(U*S2*V',s*nth*nd,[]);        % backward transformation 1
for i=1:ny
    out((i-1)*s*(nu+nth*nd+ny) + s*nu+1 : (i-1)*s*(nu+nth*nd+ny) + s*(nu+nth*nd))...
        = Fszr(:,i);                        % backward transformation 2
end
end