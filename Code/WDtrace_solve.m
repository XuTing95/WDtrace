function  [Delta, obj_fun] = WDtrace_solve(S1, S2, lambda, Prior, wpara, varargin)

% Complete algorithm of the WD-trace model which is provided in Algorithm 1 in the reference. 
% For details of the input and output, please refer to the reference.
%
% Reference: 
% T. Xu and X. F. Zhang (2017)
% Identifying gene network rewiring by integrating gene expression and gene network data
%
% COPYRIGHT  Central China Normal University
% Ting Xu <tingxu@mails.ccnu.edu.cn> 


%% Parse inputs
argin = inputParser;
argin.addRequired('S1', @(x)  isnumeric(x));
argin.addRequired('S2', @(x)  isnumeric(x));
argin.addRequired('lambda', @(x) isnumeric(x) && x>=0);
argin.addRequired('Prior', @(x) isnumeric(x));
argin.addRequired('wpara', @(x) isnumeric(x));

argin.addParamValue('alpha', 1, @(x) isnumeric(x) && x>0);
argin.addParamValue('MAX_ITER', 5e2, @(x) isnumeric(x) && x>0);
argin.addParamValue('TOL', 1e-5, @(x) isnumeric(x) && x>0);
argin.parse(S1, S2, lambda, Prior, wpara, varargin{:});

%% Copy from params object
alpha = argin.Results.alpha;
MAX_ITER = argin.Results.MAX_ITER;
TOL = argin.Results.TOL;

%% Initialization
p = size(S1,1);
Delta = zeros(p);
Delta_old = Delta;
obj_fun = 0;

Weight = Prior;
Weight( Weight==1 ) = wpara;
Weight( Weight==0 ) = 1;

%% Iteration
for iter = 1:MAX_ITER
    
    Y = Delta + (iter/(iter+3))*(Delta - Delta_old);
    
    while 1
        [f_Y, grad_Y] = Dtrace_grad(S1, S2, Y);
        Z = soft_thresh_w( Y-alpha*grad_Y, alpha*lambda, Weight, true);
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
    obj_fun = f_Delta + sum(sum( abs(lambda * Weight .* Delta) ));
    
    if iter > 1
        if (abs(obj_fun_old - obj_fun)) < (TOL * abs(obj_fun))
            break
        end
    end
    
end

Delta = Delta-diag(diag(Delta));
