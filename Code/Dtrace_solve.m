function  [Delta, obj_fun] = Dtrace_solve(S1, S2, lambda, varargin)

argin = inputParser;
argin.addRequired('S1', @(x)  isnumeric(x));
argin.addRequired('S2', @(x)  isnumeric(x));
argin.addRequired('lambda', @(x) isnumeric(x) && x>=0);

argin.addParamValue('alpha', 1, @(x) isnumeric(x) && x>0);
argin.addParamValue('MAX_ITER', 5e2, @(x) isnumeric(x) && x>0);
argin.addParamValue('TOL', 1e-5, @(x) isnumeric(x) && x>0);
argin.parse(S1, S2, lambda, varargin{:});

%% Copy from params object
alpha = argin.Results.alpha;
MAX_ITER = argin.Results.MAX_ITER;
TOL = argin.Results.TOL;

%% Initialization
p = size(S1,1);
Delta = zeros(p);
Delta_old = Delta;
obj_fun = 0;

%% Iteration
for iter = 1:MAX_ITER
    
    Y = Delta + (iter/(iter+3))*(Delta - Delta_old);
    
    while 1
        [f_Y, grad_Y] = Dtrace_grad(S1, S2, Y);
        Z = soft_thresh( Y-alpha*grad_Y, alpha*lambda, true);
        [f_Z, ~] = Dtrace_grad(S1, S2, Z);
        if f_Z <= (f_Y + sum( diag(grad_Y*(Z-Y)) ) + (norm(Z-Y,'fro'))^2 /2 /alpha)
            break
        end
        alpha = alpha * 0.5;
    end
    
    Delta_old = Delta;
    Delta = Z;
    [f_Delta, ~] = Dtrace_grad(S1, S2, Delta);
    obj_fun_old = obj_fun;
    obj_fun = f_Delta + sum(sum( abs(lambda * Delta) ));
    
    if iter > 1
        if (abs(obj_fun_old - obj_fun)) < (TOL * abs(obj_fun))
            break
        end
    end
    
end


