function [S_1, S_2,Theta_true_1, Theta_true_2, Prior, X_1, X_2] = generate_data_WDtrace(p, n, diff_ratio, pri_ratio, umin_sparse,umax_sparse, model_type, mlinks, distribution_type,  df)

switch model_type
    case 'scale_free'
        A = SFNG(p, mlinks, 1);
    case 'ER'
        A = ER_network(p);
    case 'community'
        A = community_network(p);
end

%%
%¡¡STEP 1 : Set the "dense matrix" and make it sparse

dense = (ones(p,p) * umin_sparse) + (rand(p,p) * (umax_sparse - umin_sparse));

dense = triu(dense.*A,1);
dense = dense + dense';
Net_1 = dense;
Net_2 = dense;


%% 
%¡¡STEP 2 :¡¡Generate difference between Net_1 or Net_2

A = triu(A, 1);
B = A;
edgen = p * (p-1) / 2;
[i,j] = find( triu(ones(p),1) );
ind_edge = [i,j];

diffn = floor( edgen * diff_ratio );
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


%%
% STEP 3 : Get the prior information matrix

unknown = floor( edgen * (1-pri_ratio) );

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


%%
% STEP 4 : Make all distributions "U[-umin, -umax] union U[umin, umax]" include negative parts
% (No reason to do this sooner.)

% Randomly select, 50% of values, to be negative.

flipper = unidrnd(2,p,p);
flipper(flipper == 2) = -1;
Net_1 = Net_1 .* flipper;
Net_2 = Net_2 .* flipper;


%% 
% STEP 5 : "Conditioning" the data; it must have the following Characteristics:

% 1. Symmetric
Net_1_sym = triu(Net_1) + triu(Net_1)';
Net_2_sym = triu(Net_2) + triu(Net_2)';

% check
issym=@(x) all(all(tril(x)==triu(x).'));
if ~issym(Net_1_sym)
    error('input Theta_1 is not symmetric.');
end
if ~issym(Net_2_sym)
    error('input Theta_2 is not symmetric.');
end


% 2. Positive Semi Definite
min_eigval_1 = min(eig(Net_1_sym));
min_eigval_2 = min(eig(Net_2_sym));

min_eigval = min(min_eigval_1, min_eigval_2);

% adjust the diagonals
Theta_true_1 = Net_1_sym + (eye(p)*(abs(min_eigval) + 0.1));
Theta_true_2 = Net_2_sym + (eye(p)*(abs(min_eigval) + 0.1));

% check
if (min(eig(Theta_true_1)) <= 0)
    error('input Theta_1 is not PSD.');
end
if (min(eig(Theta_true_2)) <= 0)
    error('input Theta_2 is not PSD.');
end


%%
% STEP 6 : Generate samples and get the sample covariance matrix

switch distribution_type
    case 'Gaussian'
        X_1 = mvnrnd(zeros(n,p),inv(Theta_true_1) ); 
        % Generates ssize samples from True Covariance Matrix i
        S_1 = cov(X_1); % Sample Covariance matrix i
        
        X_2 = mvnrnd(zeros(n,p),inv(Theta_true_2) ); 
        % Generates ssize samples from True Covariance Matrix i
        S_2 = cov(X_2); % Sample Covariance matrix i
        
    case 'student'
        s = diag(inv(Theta_true_1));
        if (any(s~=1))
            C = inv(Theta_true_1) ./ sqrt(s * s');
        end
        sigma = repmat(sqrt(s'),n,1);
        X_1 = mvtrnd(C,df, n);
        X_1 = sigma.*X_1;
        S_1 = cov(X_1);
        
        s = diag(inv(Theta_true_2));
        if (any(s~=1))
            C = inv(Theta_true_2) ./ sqrt(s * s');
        end
        sigma = repmat(sqrt(s'),n,1);
        X_2 = mvtrnd(C,df, n);
        X_2 = sigma.*X_2;
        S_2 = cov(X_2);
end

end
