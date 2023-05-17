clear all
close all

% load the original network (airports)
A = matfile("saveA.mat"); 
A = A.A;

% load the backbone (check in main_network_analysis.m which one)
B = matfile("saveBB.mat");
B = B.BB;

%%% ---- TEST AREA --------------------------------------------------------

J = Jaccard_weighted(A,B)

% G = digraph(A);
% 
% D = 1./A;
% D(D==inf) = 0;
% 
% [P, d] = shortestpath(G, 1, 3001 )
% 
% T = get_SP_tree(15167, D);
% 
% sigma = indicator(A);

% tic
% S = salience(B);
% toc


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


function [J] = Jaccard(W, B)
% this function compute the Jaccard similarity between two network whose 
% adiacency matrices are W and B (ideally W is the original network and B 
% is the backbone filtered. 
% For each node, it compares the links of the node in W and in B, and then
% averages over all the nodes. 
% W and B are sparse objects
    N = length(W);
    %loop over nodes
    jacc = [];
    for i=1:N
        %disp(i);

        x = find(W(i,:)~=0); % links in W for node i, 
        % this is an array of all the indeces of the nodes i is connected to 
        y = find(B(i,:)~=0); % links in B for node i

        xx = zeros(N,1);   % array with length equal to the total number of nodes
        xx(x) = 1;       % 1 means there is a link 
        yy = zeros(N,1);
        yy(y) = 1;

        a = sum( (xx+yy)>=2 ); % links for node i in both W and B
        b = sum( (xx-yy)>0 );  % links exclusively in W
        c = sum( (yy-xx)>0 );  % links exclusively in B

        if (a+b+c) > 0
            jacc(i) = a / (a + b + c); 
        else
            jacc(i) = 0;
        end
    end

    J =  1/N * sum(jacc);
end



function [S] = salience(W)
% this function computes the salience of each link in the W network

    N = length(W(:,1)); % number of nodes
    S = sparse(N,N); % matrix of saliences (s_ij = salience of link ij), initialized to 0 
    
    % compute the inverse weights matrix
    D = 1./W;
    D(D==inf) = 0;

    G = digraph(D);

    % loop over the beginning node of the path (the reference node x) 
    for x=1:N
        disp(x);

        %check if the node is connected to something
        if max(D(x,:))~=0 

            %create the shortest path matrix for the reference node x
            T = sparse(N,N); 

            % loop over the ending node of the path, among all the N nodes
            for y=1:N

                %check if the node is at the end of some link 
                if max(D(:,y))~=0
    
                    % compute the shortest path
                    [P,d] = shortestpath(G,x,y);
                    
                    % add the links in the shortestpath to the T matrix
                    for k=1:(length(P)-1)

                        T(P(k),P(k+1)) = 1;
                    end   
                end
            end

            % add the component of the salience related to the reference node x
            % to the salience vector
            S = S + T;

        end
    end

    % average the sum for each link to get the salience 
    S = S ./ N;
end




function [sigma] = indicator(W)
% this function creates a cell array containing the shortest path between
% any two nodes in the network. The entrance sigma_ij = path array that
% goes from node i to node j,
% it's still a little slow

    N = length(W(:,1)); % number of nodes
    sigma = {}; % cell matrix
    
    % compute the inverse weights matrix
    D = 1./W;
    D(D==inf) = 0;

    G = digraph(D);

    % loop over the beginning node of the path
    for x=1:N
        disp(x);

        %check if the node is connected to something
        if max(D(x,:))~=0 

            % loop over the ending node of the path, among all the N nodes
            for y=1:N

                %check if the node is at the end of some link 
                if max(D(:,y))~=0
    
                    % compute the shortest path
                    [P,d] = shortestpath(G,x,y);
        
                    % add each link in the SP to the sigma matrix
                    sigma{x,y}= P;
                end
            end
        end
    end

end



function [S] = salience_old(i, j, W)
% compute the salience of the link between node i and j 

    D = 1./W; % compute the distances matrix 
    D(D==inf) = 0;

    N = length(W(:,1)); % number of nodes

    for r = 1:N   % loop on all nodes 
        T = get_SP_tree(r, D);  % compute the shortest path tree matrix 
        S = S + T; % add the new terms to salience
    end

    S = S / N; % compute the mean                               
end

function [T] = get_SP_tree(r, D)
% compute the shortest path tree matrix of the reference node r, 
% D is the distance matrix of the network, where d_ij = 1 / w_ij
    
    N = length(D(:,1)); % number of nodes

    G = digraph(D); % put it in a graph form  

    P = zeros(N, 15 );     % initialise shortest paths matix
    d = [];                % initialise distances vector

    for i = 1:N % compute shortest path between node r and all the other nodes in the network
        
        if mod( i, N/10 ) == 0
            disp(i);
        end
        
        [p, d(i)] = shortestpath( G, r, i );
        P( i, 1:length(p) ) = p;
        
    end

    T = zeros( N, N ); % initialize T(r) matrix

    for i = 1:N
        disp(i)
        aux = find( P==i );
        
        if aux>0

            for j = 1:N
                aux2 = ( P(aux+1)==j );
                if aux2 > 0
                    T(i,j) = 1;
                end
            end
        end
    end 
end