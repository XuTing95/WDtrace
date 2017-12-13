% --------------------- Weighted SOFT-THRESHOLDING OPERATOR ------------------- %
%
% -------------- minimize 1/2*||X - Y||_2^2 + lambda*Weight.*||X||_1 ---------- %
%
% -------------------------- LAST UPDATE: 12/13/2016 -------------------------- %
%
%
% Reference: 
% T. Xu and X. F. Zhang (2017)
% Identifying gene network rewiring by integrating gene expression and gene network data
%
% COPYRIGHT  Central China Normal University
% Ting Xu <tingxu@mails.ccnu.edu.cn> 


function X = soft_thresh_w(Y, lambda, Weight, penalize_diagonal)

if nargin == 1
    disp('Not enough inputs');
    disp('Enter both Y and lambda');
    return
end

X = sign(Y).*(max(abs(Y) - lambda.*Weight, 0));

if ~penalize_diagonal
    X = X - diag(diag(X)) + diag(diag(Y));
end
