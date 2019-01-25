function [z,error,mec] = solver(y,x, Params,A,tol)
%% Initialization

z       = init(y,Params.L,A,x);
error   = norm(y-abs(A(z)),'fro')/norm(y,'fro');
mec     = metric(x,z);
fprintf('Iter = 0  Error = %f Metric = %f\n', error,mec);

u = Params.u0;

%% main loop
sgd    = Params.m;
batch  = Params.B;
gamma  = Params.y;
gamma1 = Params.y1;
mu     = Params.mu;

for t = 1: Params.T
    for i=1:batch:sgd-batch+1
        
        grad = compute_grad(z,u,y,Params,A,i);
        
        z = z - mu*grad;
        
        if norm(grad,'fro') < gamma*u
            u = gamma1*u;
        end
    end
    error(t + 1)   = norm(y-abs(A(z)),'fro')/norm(y,'fro');
    mec(t+1)       = metric(x,z);
    
    if norm(grad,'fro')<tol
        break;
    end
    
    fprintf('Iter = %d  Error = %f Metric = %f normgrad = %f mu = %f\n',t, error(end),mec(end),norm(grad,'fro'),u);
end
return