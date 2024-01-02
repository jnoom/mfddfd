function out = proxls(x,gamma,Y,H,Hinv)
if nargin <5
    out = (gamma*(H'*H) + eye(length(H(1,:))))\(x + gamma*H'*Y);
else
    out = Hinv*(x + gamma*H'*Y);
end
end