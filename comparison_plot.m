clear all 
close all

%%% This file aims to plot measurement of the different performance of the
%%% 3 filters implemented and, in function of the multivariate
%%% significance level

% % load the original network (airports)
% A = matfile("saveA.mat"); 
% A = A.A;
% 
% % GET THE NUMBER OF NODES CONNECTED IN THE NETWORK
% 
% N = length(A(:,1)); % total number of nodes
% 
% % compute the real number of nodes connected in the network
% [row,col] = find(A>0);
% new = [row col];
% N_A = length(unique(new));
% 
% % GET THE NUMBER OF LINK IN THE NETWORK
% L = nnz(A); 


% % load the ER weighted network
% ER = load("saveweightedER.mat")
% ER = ER.W;
% 
% load the BA weighted directed network 
A = load("saveweightedBA.mat");
A = A.W;

% GET THE NUMBER OF NODES CONNECTED IN THE NETWORK

N = length(A(:,1)); % total number of nodes

% compute the real number of nodes connected in the network
aux = find(A>0);
[row,col] = ind2sub( [N, N], aux);
new = [row col];
N_A = length(unique(new));

% GET THE NUMBER OF LINK IN THE NETWORK
L = nnz(A); 


%%% FRACTION OF NODES RETAINED %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initialize vectors of nodes retained for 
fr_node_DF = [];       % disparity filter
fr_node_HF = [];       % hypergeometric filter
fr_node_PF_min = [];   % polya filter min 
fr_node_PF_max = [];   % polya filter max 
fr_node_PF_1 = [];     % polya filter a=1 

%%% FRACTION OF LINKS RETAINED %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initialize vectors of links retained for 
fr_link_DF = [];       % disparity filter
fr_link_HF = [];       % hypergeometric filter
fr_link_PF_min = [];   % polya filter min 
fr_link_PF_max = [];   % polya filter max
fr_link_PF_1 = [];     % polya filter a=1

%%% JACCARD SIMILARITY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initialize vectors of jaccard similarity for 
jacc_DF = [];       % disparity filter
jacc_HF = [];       % hypergeometric filter
jacc_PF_min = [];   % polya filter min 
jacc_PF_max = [];   % polya filter max
jacc_PF_1 = [];     % polya filter a=1



% generate 30 values of multivariate significance level equally spaced 
% between 10^-8 and 10^-1
alpha = logspace(-8, 0, 10);


apr_lvl = 10;   % approximation level for the polya filter
a1 = 0.01;      % min val for a in polya filter 
a2 = 6;         % max val for a in polya filter 


for i = 1:length(alpha)
    tic
    disp(i)
    disp("Disparity filter")
    % Computing backbone with disparity filter function (edge list)
    backbone = disp_filter(A,alpha(i)); 

    fr_node_DF(i) = length(unique(backbone)) / N_A;
    fr_link_DF(i) = length(backbone) / L;

    % converting edge list to matrix
    if backbone > 0
        ind = sub2ind(size(A), backbone(:,1), backbone(:,2));
        BB = sparse(backbone(:,1), backbone(:,2), A(ind));
        % pad the matrix with zeros to make it square
        BB = sparse([full(BB) zeros(size(BB, 1), size(A,1)-size(BB,2))]);
        BB = sparse([full(BB); zeros(size(A,2)-size(BB,1), size(BB,2))]);
    
        jacc_DF(i) = Jaccard_weighted(A, BB);
    end

    disp("Hypergeometric filter")
    % Computing backbone with hypergeometric filter function (edge list)
    backbone = hypergeom_filter(A,alpha(i)); 

    fr_node_HF(i) = length(unique(backbone)) / N_A;
    fr_link_HF(i) = length(backbone) / L;

    % converting edge list to matrix
    if backbone > 0
        ind = sub2ind(size(A), backbone(:,1), backbone(:,2));
        BB = sparse(backbone(:,1), backbone(:,2), A(ind));
        % pad the matrix with zeros to make it square
        BB = sparse([full(BB) zeros(size(BB, 1), size(A,1)-size(BB,2))]);
        BB = sparse([full(BB); zeros(size(A,2)-size(BB,1), size(BB,2))]);
    
        jacc_HF(i) = Jaccard_weighted(A, BB); 
    end

    disp("Polya filter, a min")
    % Computing backbone with polya filter function (edge list)
    backbone = polya_filter(A, a1, alpha(i), apr_lvl);

    fr_node_PF_min(i) = length(unique(backbone)) / N_A;
    fr_link_PF_min(i) = length(backbone) / L;

    if backbone > 0
        % converting edge list to matrix
        ind = sub2ind(size(A), backbone(:,1), backbone(:,2));
        BB = sparse(backbone(:,1), backbone(:,2), A(ind));
        % pad the matrix with zeros to make it square
        BB = sparse([full(BB) zeros(size(BB, 1), size(A,1)-size(BB,2))]);
        BB = sparse([full(BB); zeros(size(A,2)-size(BB,1), size(BB,2))]);
    
        jacc_PF_min(i) = Jaccard_weighted(A, BB);
    end

    disp("Polya filter, a max")
    backbone = polya_filter(A, a2, alpha(i), apr_lvl);

    fr_node_PF_max(i) = length(unique(backbone)) / N_A;
    fr_link_PF_max(i) = length(backbone) / L;

    if backbone > 0
        % converting edge list to matrix
        ind = sub2ind(size(A), backbone(:,1), backbone(:,2));
        BB = sparse(backbone(:,1), backbone(:,2), A(ind));
    
        % pad the matrix with zeros to make it square
        BB = sparse([full(BB) zeros(size(BB, 1), size(A,1)-size(BB,2))]);
        BB = sparse([full(BB); zeros(size(A,2)-size(BB,1), size(BB,2))]);
    
        jacc_PF_max(i) = Jaccard_weighted(A, BB);
    end 


    disp("Polya filter, a=1")
    backbone = polya_filter(A, 1, alpha(i), apr_lvl);

    fr_node_PF_1(i) = length(unique(backbone)) / N_A;
    fr_link_PF_1(i) = length(backbone) / L;

    if backbone > 0
        % converting edge list to matrix
        ind = sub2ind(size(A), backbone(:,1), backbone(:,2));
        BB = sparse(backbone(:,1), backbone(:,2), A(ind));
    
        % pad the matrix with zeros to make it square
        BB = sparse([full(BB) zeros(size(BB, 1), size(A,1)-size(BB,2))]);
        BB = sparse([full(BB); zeros(size(A,2)-size(BB,1), size(BB,2))]);
    
        jacc_PF_1(i) = Jaccard_weighted(A, BB);
    end
    toc
end

%%% PLOT GRAPH WITH FRACTION OF RETAINED NODES VS ALPHA %%%%%%%%%%%%%%%%%%%
semilogx(alpha, fr_node_DF, alpha, fr_node_HF, alpha, ...
    fr_node_PF_max, alpha, fr_node_PF_min, alpha, fr_node_PF_1, 'LineWidth',1.5);

hold on 
patch([alpha fliplr(alpha)], [fr_node_PF_min fliplr(fr_node_PF_max)], 'b', 'FaceAlpha',.3)

xlabel("alpha")
ylabel("Fraction of nodes retained")
legend('DF', 'HF', 'PF_{MAX}', 'PF_{MIN}','PF (a=1)', 'Location','northwest')

savefig(gcf, 'Fr_nodes_BA.fig')

hold off 

%%% PLOT GRAPH WITH FRACTION OF RETAINED LINKS VS ALPHA %%%%%%%%%%%%%%%%%%%

semilogx(alpha, fr_link_DF, alpha, fr_link_HF, ...
    alpha, fr_link_PF_max, alpha, fr_link_PF_min, alpha, fr_link_PF_1, 'LineWidth',1.5);

hold on 
patch([alpha fliplr(alpha)], [fr_link_PF_min fliplr(fr_link_PF_max)], 'b', 'FaceAlpha',.3)

xlabel("alpha")
ylabel("Fraction of links retained")
legend('DF', 'HF', 'PF_{MAX}', 'PF_{MIN}', 'PF (a=1)', 'Location','northwest')
savefig(gcf, 'Fr_links_BA.fig')

hold off 


%%% PLOT GRAPH WITH JACCARD SIMILARITY VS ALPHA %%%%%%%%%%%%%%%%%%%%%%%%%%%

semilogx(alpha, jacc_DF, alpha, jacc_HF, ...
    alpha, jacc_PF_max, alpha, jacc_PF_min, alpha, jacc_PF_1, 'LineWidth',1.5);

hold on 
patch([alpha fliplr(alpha)], [jacc_PF_min fliplr(jacc_PF_max)], 'b', 'FaceAlpha',.3)

xlabel("alpha")
ylabel("Jaccard similarity")
legend('DF', 'HF', 'PF_{MAX}', 'PF_{MIN}', 'PF (a=1)', 'Location','northwest')
savefig(gcf, 'Jacc_BA.fig')

hold off 


%%% FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [J] = Jaccard_weighted(W,B)
% this function compute the weighted Jaccard similarity between two network 
% whose adiacency matrices are W and B (ideally W is the original network 
% and B is the backbone filtered. 
% For each node, it compares the weights of outgoing links of the node in W
% and in B, and then averages over all the nodes that has at least one 
% outgoing link. 
% W and B are sparse objects

    N = length(W); % total number of nodes

    % compute the real number of nodes connected in the network
    [row,col] = find(W>0);
    %new = [row col];

    % list of indeces of nodes which has at least an outgoing link
    ind = unique(row);  
    N_A = length(ind);

    % initialize vector of jaccard similarities for each node 
    jacc = [];

    % loop over connected nodes
    for i = 1:N_A

        jacc(i) = sum(min( W(ind(i),:), B(ind(i),:) )) / sum( max(W(ind(i),:), B(ind(i),:)) );
    end
    % compute the average
    J = mean(jacc);

end


