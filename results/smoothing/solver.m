function [z,error] = solver(y,x,x0, Params,A,tol)
%% Initialization
if isempty(x0)
    z  = init(y,Params.L,A);
else
    z = x0;
end
error   = norm(y-abs(A(z)),'fro')/norm(y,'fro');
fprintf('Iter = 0  Error = %f \n', error);

u = Params.u0;

%% main loop
sgd    = Params.m;
batch  = Params.B;
gamma  = Params.y;
gamma1 = Params.y1;
mu     = Params.mu;

for t = 1: Params.T
    for i=1:batch:sgd-batch+1
        
        grad = compute_grad(z,u,y,Params,i);
        
        z = z - mu*grad;
        
        if norm(grad,'fro') < gamma*u
            u = gamma1*u;
        end
    end
    error(t + 1)   = norm(y-abs(A(z)),'fro')/norm(y,'fro');
    
    if error(end)<tol
        break;
    end
    
    fprintf('Iter = %d  Error = %f normgrad = %f mu = %f\n',t, error(end),norm(grad,'fro'),u);
end
return