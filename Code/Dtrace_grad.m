function [f, grad] = Dtrace_grad(S1, S2, Delta)

% Compute the value and gradient of D-trace loss function. 
% For details of the input and output, please refer to the reference.
%
% Reference: 
% T. Xu and X. F. Zhang (2017)
% Identifying gene network rewiring by integrating gene expression and gene network data
%
% COPYRIGHT  Central China Normal University
% Ting Xu <tingxu@mails.ccnu.edu.cn> 


    b = S1 - S2;
    A = S1 * Delta * S2;
    B = S2 * Delta * S1;
    grad = 0.5 * (A+B) - b;
    f = 0.5 * sum(diag(Delta * A)) - sum(diag(Delta * b));
    
end