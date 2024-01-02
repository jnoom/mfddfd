function out = proxls_fw(x,gamma,rk,Pk,ny)
x2  = reshape(x,[],ny);
rk2 = reshape(rk,[],ny);
n   = length(x2(:,1));
for i=1:ny
    out((i-1)*n+1:i*n,1) = Pk*(x2(:,i) + gamma*rk2(:,i));
end
end