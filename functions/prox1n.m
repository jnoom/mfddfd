function out = prox1n(x,gamma,lambda,s,nu,ny,nd,nth)
out = x;  % keep original values for part of the vector at which the 1-norm does not apply

for i=1:ny
    Fz(:,i) = x((i-1)*s*(nu+nth*nd+ny) + s*nu+1 : (i-1)*s*(nu+nth*nd+ny) + s*(nu+nth*nd));  % forward transormation 1
end
Fzv = reshape(Fz,[],1);  % forward transormation 2
    
Fzvr = reshape(sign(Fzv).*(max(0,abs(Fzv) - gamma*lambda)),s*nth*nd,[]);  % soft-thresholding & backward transformation 1

for i=1:ny
    out((i-1)*s*(nu+nth*nd+ny) + s*nu+1 : (i-1)*s*(nu+nth*nd+ny) + s*(nu+nth*nd))...
        = Fzvr(:,i);     % backward transformation 2
end
end