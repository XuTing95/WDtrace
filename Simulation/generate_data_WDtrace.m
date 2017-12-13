function [S_1, S_2, Delta_true, Prior] = generate_data_WDtrace(p, n, diff_rate, pri_rate, umin_sparse,umax_sparse)

% Generate the simulation data using the procedure provised in Section III-A of the reference.
% For details of the input and output, please refer to the reference.
%
% Reference: 
% T. Xu and X. F. Zhang (2017)
% Identifying gene network rewiring by integrating gene expression and gene network data
%
% COPYRIGHT  Central China Normal University
% Ting Xu <tingxu@mails.ccnu.edu.cn> 


%% STEP 1
A = SFNG(p, 2, 1);


%% STEP 2
dense = (ones(p,p) * umin_sparse) + (rand(p,p) * (umax_sparse - umin_sparse));
dense = triu(dense.*A,1);
dense = dense + dense';
Net_1 = dense;
Net_2 = dense;


%% STEP 3
A = triu(A, 1);
B = A;
edgen = p * (p-1) / 2;
[i,j] = find( triu(ones(p),1) );
ind_edge = [i,j];

diffn = floor( edgen * diff_rate );
diff_edge = sort( randsample(edgen, diffn) );
Net_2 = triu(Net_2, 1);
for k = 1:diffn
    ind = ind_edge(diff_edge(k),:);
    if A(ind(1),ind(2))==0
        B(ind(1),ind(2)) = 1;
        Net_2(ind(1),ind(2)) = umin_sparse + rand(1) * (umax_sparse-umin_sparse);
    else
        B(ind(1),ind(2)) = unidrnd(2) - 1;
        Net_2(ind(1),ind(2)) = B(ind(1),ind(2)) * ( umin_sparse + rand(1) * (umax_sparse-umin_sparse) );
    end
end
Net_2 = Net_2 + Net_2';


%% STEP 4
unknown = floor( edgen * (1-pri_rate) );

unknow1_edge = sort( randsample(edgen, unknown) );
for k = 1:unknown
    ind = ind_edge(unknow1_edge(k),:);
    A(ind(1),ind(2)) = 0;
end

unknow2_edge = sort( randsample(edgen, unknown) );
for k = 1:unknown
    ind = ind_edge(unknow2_edge(k),:);
    B(ind(1),ind(2)) = 0;
end

A = A + A';
B = B + B';
Prior = ( A | B ) * 1;


%% check
flipper = unidrnd(2,p,p);
flipper(flipper == 2) = -1;
Net_1 = Net_1 .* flipper;
Net_2 = Net_2 .* flipper;

Net_1_sym = triu(Net_1) + triu(Net_1)';
Net_2_sym = triu(Net_2) + triu(Net_2)';
issym=@(x) all(all(tril(x)==triu(x).'));
if ~issym(Net_1_sym)
    error('input Theta_1 is not symmetric.');
end
if ~issym(Net_2_sym)
    error('input Theta_2 is not symmetric.');
end


%% STEP 5
min_eigval_1 = min(eig(Net_1_sym));
min_eigval_2 = min(eig(Net_2_sym));
min_eigval = min(min_eigval_1, min_eigval_2);

Theta_true_1 = Net_1_sym + (eye(p)*(abs(min_eigval) + 0.1));
Theta_true_2 = Net_2_sym + (eye(p)*(abs(min_eigval) + 0.1));

if (min(eig(Theta_true_1)) <= 0)
    error('input Theta_1 is not PSD.');
end
if (min(eig(Theta_true_2)) <= 0)
    error('input Theta_2 is not PSD.');
end


%% STEP 6
X_1 = mvnrnd(zeros(n,p),inv(Theta_true_1) ); 
S_1 = cov(X_1);
X_2 = mvnrnd(zeros(n,p),inv(Theta_true_2) ); 
S_2 = cov(X_2); 
Delta_true = Theta_true_2 - Theta_true_1;
Delta_true = Delta_true - diag(diag(Delta_true));

end
