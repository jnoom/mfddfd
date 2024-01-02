function out = ppxa(Y,H,x0,lambda1,lambda2,gamma,s,nu,ny,nd,nth,Nmax)
v(:,1:3) = repmat(x0,1,3);
x        = mean(v(:,:,1),2);
Hinv     = inv(3*gamma*(H'*H) + eye(length(H(1,:))));

j = 2;
while j < Nmax
    p(:,1) = proxls(v(:,1),3*gamma,Y,H,Hinv);
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