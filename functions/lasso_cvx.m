function z = lasso_cvx(Y,Tu,Ty,Tth,Bv,Fv,Kv,lambda,nth,ny)
    cvx_begin quiet
        variables z(nth,1)
        ls_term = 0;
        for i=1:ny
            ls_term = ls_term + sum_square_abs(Y(:,i)-[Tu,Ty]*[Bv(:,i);Kv(:,i)] - Tth*kron(Fv(:,i),z));
        end
        minimize( ls_term + lambda*norm(z,1) )
    cvx_end
end