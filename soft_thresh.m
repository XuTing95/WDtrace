%% -------------------- SOFT-THRESHOLDING OPERATOR ------------------- %%


% ---------------- minimize 1/2*||X - Y||_2^2 + lambda*||X||_1 -------- %



% ----------------- LAST UPDATE: 4/19/2012 ---------------------------- %

function X = soft_thresh(Y,lambda, penalize_diagonal)

if nargin == 1
    disp('Not enough inputs');
    disp('Enter both Y and lambda');
    return
end

X = sign(Y).*(max(abs(Y) - lambda,0));

if ~penalize_diagonal
    X = X - diag(diag(X)) + diag(diag(Y));
end
