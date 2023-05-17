clear all
close all

% % LOAD AIRPORT NETWORK ----------------------------------------------------
 T = readtable('T_T100D_MARKET_ALL_CARRIER.csv');
% % load the air traffic data form US database of 2017
% % col 1: passengers carried (weight)
% % col 2: distance 
% % col 5: origin airport ID (numeric)
% % col 6: origin airport ID (letters)
% % col 8: destination airport ID (numeric) 
% % col 9: destination airport ID (letters) 
% 
A = sparse(T{:,"ORIGIN_AIRPORT_ID"}, T{:,"DEST_AIRPORT_ID"}, T{:,"PASSENGERS"});
save('saveA.mat', 'A')
% -------------------------------------------------------------------------


alpha = 0.05;   % univariate significance level
a = 1;          % Polya parameter
apr_lvl = 10;   % approximation level for the polya filter



% k = full(sum(A>0)); % Degree sequence (in-degree) (sum over column)
% s = full(sum(A)); % Stregth degree (in-stregth) 
% 
%%% PLOT THE 10 HEAVIEST LINK IN THE NETWORK %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [ind1,ind2] = find(A > 0);
% 
% [sk,ind_sort] = sort(k); % Sorted degree sequence
% 
% % Printing to screen the IDs of the top 10 networks in terms of in-degree
% for h = 0:9
%    
%     fprintf('Degree = %d, ID = %d\n',sk(end-h),ind_sort(end-h))
%     
% end
% 
% fprintf('\n')


L = nnz(A); % Computing the number of links in A

alpha = alpha/L; % Bonferroni correction


% % load the BA weighted directed network %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BA = load("saveweightedBA.mat");
% BA = BA.W;
% 
% % GET THE NUMBER OF NODES CONNECTED IN THE NETWORK
% 
% N = length(BA(:,1)); % total number of nodes
% 
% % compute the real number of nodes connected in the network
% aux = find(BA>0);
% [row,col] = ind2sub( [N, N], aux);
% new = [row col];
% N_A = length(unique(new));
% 
% % GET THE NUMBER OF LINK IN THE NETWORK
% L = nnz(BA); 
% 
% alpha = 0.1;   % univariate significance level
% a = 0.01;          % Polya parameter
% apr_lvl = 10;   % approximation level for the polya filter
% 
% alpha = alpha/L; % Bonferroni correction


%%% BACKBONE EXTRACTION (LIST OF EDGES) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%backbone = disp_filter(A,alpha);  
%backbone = hypergeom_filter(A,alpha);  
backbone = polya_filter(A, a, alpha, apr_lvl);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%backbone = disp_filter(BA,alpha);       
%backbone = hypergeom_filter(BA,alpha); 
%backbone = polya_filter(BA, a, alpha, apr_lvl);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% list of indeces of the validated link in the backbone, in the form of a
% single index
ind = sub2ind(size(A), backbone(:,1), backbone(:,2));

%BACKBONE 
BB = sparse(backbone(:,1), backbone(:,2), A(ind));

% pad the matrix with zeros to make it the same dimension of the original
% network
BB = sparse([full(BB) zeros(size(BB, 1), size(A,1)-size(BB,2))]);
BB = sparse([full(BB); zeros(size(A,2)-size(BB,1), size(BB,2))]);

save('saveBB.mat', 'BB')

% G = digraph(BB);
% 
% P = shortestpath(G, 15167, 15000)
% 
% plot(G, 'Layout','circle')

