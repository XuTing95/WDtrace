%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 		                     	                                                                     %
% 		          TEST FILE FOR WD-trace using Breast Cancer gene expression data                    %
% 		                     	                                                                     %
% 		                     	                                                                     %
%  Refer to the paper:  T. Xu and X. F. Zhang (2017)                                                 %
%     Identifying gene network rewiring by integrating gene expression and gene network data         %
%                                                                                                    %
%    CONTACT   Ting Xu (tingxu@mails.ccnu.edu.cn) for any questions or comments on the code.	     %
% 		                     	                                                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
close all


%% load data and set parameter
load BRCA
load BRCA_Prior
Gene = BRCA.Gene;
X1 = zscore(BRCA.data{1,1});
X2 = zscore(BRCA.data{2,1});
S1 = cov(X1);
S2 = cov(X2);

w = 0.40;
lambda  = 0.49;


%% estimate the differential networks using WD-trace with lambda = 0.34, w = 0.35
[Delta_WD, ~] = WDtrace_solve(S1, S2, lambda, Prior, w);


%% get the differential edges estimated by WD-trace
[i,j] = find( Delta_WD~=0 );
kk = i>j;
Edge_WD = [Gene(i(kk))',Gene(j(kk))'];

%% get the top 10 hub genes of WD-trace estimated differential network
result_WD = tabulate(Gene(i));
[~,id_count] = sort(cell2mat(result_WD(:,2)),'descend');
result_WD = result_WD(id_count,:);
Hub_Genes_WD = result_WD(id_count(1:10),1);

