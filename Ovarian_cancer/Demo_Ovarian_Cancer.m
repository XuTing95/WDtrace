%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 		                     	                                                                     %
% 		         TEST FILE FOR WD-trace using Ovarian Cancer gene expression data                    %
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
load OV
load OV_Prior

Gene = OV.Gene;
S1 = OV.cov{1,1};
S2 = OV.cov{2,1};

w = 0.35;
lambda  = 0.34;


%% estimate the differential networks using WD-trace with lambda = 0.34, w = 0.35
[Delta_WD, ~] = WDtrace_solve(S1, S2, lambda, Prior, w);


%% get the differential edges estimated by WD-trace
[i,j] = find( Delta_WD~=0 );
kk = i>j;
Edge_WD = [Gene(i(kk)),Gene(j(kk))];

%% get the top 10 hub genes of WD-trace estimated differential network
result_WD = tabulate(Gene(i));
[~,id_count] = sort(cell2mat(result_WD(:,2)),'descend');
result_WD = result_WD(id_count,:);
Hub_Genes_WD = result_WD(id_count(1:10),1);
