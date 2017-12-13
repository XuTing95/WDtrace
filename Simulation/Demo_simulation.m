%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 		                     	                                                                     %
% 		                   TEST FILE FOR WD-trace using simulated data                               %
% 		                     	                                                                     %
% 		                     	                                                                     %
%  Refer to the paper:  T. Xu and X. F. Zhang (2017)                                                 %
%     Identifying gene network rewiring by integrating gene expression and gene network data         %
%                                                                                                    %
%    CONTACT   Ting Xu (tingxu@mails.ccnu.edu.cn) for any questions or comments on the code.	     %
% 		                     	                                                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% experiment settings
clear 
close all 

p = 50;
n = 100;
diff_rate = 0.1;
pri_rate = 0.8;
umin_sparse = 0.3;
umax_sparse = 0.6;
lambda = 0.5;
w = 0.1;

%% generate simulation data
[S_1, S_2, Delta_true, Prior] = generate_data_WDtrace(p, n, diff_rate, pri_rate, umin_sparse, umax_sparse);


%% estimate the differential netwoks using WD-trace
[Delta, ~] = WDtrace_solve(S_1, S_2, lambda, Prior, w);


%% Compare estimated networks to true networks
subplot(1,2,1);
imagesc(Delta_true);
colorbar
title('True \Delta');

subplot(1,2,2);
imagesc(Delta)
colorbar
title('Estimated \Delta');


%% Calculate the Pre, Rec, TPR and FPR of estimated result
Diff_true = abs(Delta_true) > 10^(-10);
Diff = abs(Delta) >  10^(-10);
TP = sum(sum(triu(Diff_true,1).* triu(Diff,1)));
FP = sum(sum(triu(~Diff_true,1).* triu(Diff,1)));
TN = sum(sum(triu(~Diff_true,1).* triu(~Diff,1)));
FN = sum(sum(triu(Diff_true,1).* triu(~Diff,1)));

Pre = TP / (TP + FP+eps);
Rec = TP / (TP + FN+eps);
TPR = TP / (TP + FN+eps);
FPR = FP / (FP + TN+eps);
