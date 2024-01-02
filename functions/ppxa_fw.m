function out = ppxa_fw(Pk,rk,x0,lambda1,lambda2,gamma,s,nu,ny,nd,nth,Nmax)
v(:,1:3) = repmat(x0,1,3);
x        = mean(v(:,:,1),2);

j = 2;
while j < Nmax
    p(:,1) = proxls_fw(v(:,1),3*gamma,rk,Pk,ny);
    p(:,2) = proxnn(v(:,2),3*gamma,lambda1,s,nu,ny,nd,nth);
    p(:,3) = prox1n(v(:,3),3*gamma,lambda2,s,nu,ny,nd,nth);
    pm(:)  = mean(p,2);
    for i=1:3
        v(:,i) = v(:,i) + 2*pm(:) - x - p(:,i);
    end
    x = x + pm(:) - x;
    j = j+1;
end
out = x;
end