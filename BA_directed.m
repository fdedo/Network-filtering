%%% THIS CODE GENERATE A DIRECTED BARABASI-ALBERT NETWORK
%%% It uses the copying mechanisms to generate new links directed toward
%%% already existing nodes, and for each link it choose randomly if
%%% directing it from the new node or towards it.

%%% THINGS YET TO IMPLEMENT:
    % A fitness model
    % The possibility to generate links between two already existing nodes
    % A boost for new links
    % Some sort of ageing 
    % ...

clear all
close all

N = 3000;   % Number of final nodes
m0 = 6;     % Starting fully connected nodes
m = 3;      % Number of new links added at each step
p = 1;   % probability to connect to the randomly chosen link (1-p is the prob to copy one of the links of the chosen node)
alfa = 0.5; % probability of generating a link outgoing from the new node

%A = BA_network(N,m);

A = BA_Dnet(N,m0, m, p, alfa);

save('saveBA.mat', 'A');


%%% PLOT THE DEGREE DISTRIBUTION
k_in = full(sum(A));
k_out = full(sum(A'));

% incoming degree distribution
h_in = histogram(k_in,100,'Normalization','pdf');
x_in = h_in.BinEdges;
x_in = x_in(2:end);
y_in = h_in.Values;

f = find(y_in == 0);
x_in(f) = []; y_in(f) = [];

% outgoing degree distribution
h_out = histogram(k_out, 100, 'Normalization','pdf');
x_out = h_out.BinEdges;
x_out = x_out(2:end);
y_out = h_out.Values;

f = find(y_out == 0);
x_out(f) = []; y_out(f) = [];

subplot(2,1,1)
loglog(x_in,y_in,'ob','MarkerSize',6,'MarkerFaceColor','b')
hold on
% fit con una power law
% OCCHIO CHE I COEFFICIENTI NON SONO QUELLI CHE TROVA LUI 
c_in = polyfit(log(x_in),log(y_in),1);
loglog(x_in, 50*exp(c_in(2))*x_in.^-3, 'b-','LineWidth',1.5)

subplot(2,1,2)
loglog(x_out,y_out,'ob','MarkerSize',6,'MarkerFaceColor','r')
hold on 

%fit con una power law
% OCCHIO CHE I COEFFICIENTI NON SONO QUELLI CHE TROVA LUI 
c_out = polyfit(log(x_out),log(y_out),1);
loglog(x_out, 50*exp(c_out(2))*x_out.^-3, 'r-','LineWidth',1.5)

xlabel('$k$','Interpreter','latex')
ylabel('$p(k)$','Interpreter','latex')
set(gca,'FontSize',22)

%close


G = digraph(A);
plot(G, 'Layout','force')


%%% function to generate a undirected BA network 
function A = BA_network(N,m)    


        % Initialising adjacency matrix of initial core (fully connected)
        A = ones(m) - eye(m);
        A = sparse(A);
        
        % Loop on number of nodes to be added 
        for i = 1:N-m
        
            k = full(sum(A)); % Degree sequence
            p = k/sum(k); % Probability of linking to a node
        
            %%% Our strategy to simulate preferential attachement will be
            %%% to convert the vector p of probabilities into a vector of
            %%% cumulative probabilities (i.e., a vector that sums up to 1,
            %%% which is akin to the partition of the unit interval into k
            %%% sub-intervals, whose length is proportional to the
            %%% probability of connecting to nodes).
            
            p = cumsum(p); 
            r = rand(m,1); % m random numbers in [0,1]
        
            %%% Adding new nodes and links
            ind = [];
        
            %%% Loop on new links formed by new node: the links are formed
            %%% with the preexisting nodes whose corresponding
            %%% sub-interval in the vector p contains the random numbers in
            %%% r
            for j = 1:m
                        
                aux = p - r(j);
                aux = find(aux > 0);
                ind = [ind; aux(1)];
        
            end
        
            ind = unique(ind)
        
            %%% Creating new rows and columns in adjacency matrix
            A = [A; zeros(1,size(A,2))];
            A = [A zeros(size(A,1),1)];
        
            %%% Adding new links
            A(end,ind) = 1;
            A(ind,end) = 1;
        
        end     

end

%%% QUESTIONS:
%%% WHEN INITIALISING THE NETWORK, IS IT OK TO USE THE SAME FULL CONNECTED
%%% STARTING POINT AS UNIDIRECTED?
%%% WHEN ADDING NEW NODES, IS IT OK TO INTRUDUCE THEM ONLY WITH OUTGOING
%%% LINK TO PRE-EXISTING NODES?


function A = BA_Dnet(N, m0, m, p, alfa)
%%% function to generate a DIRECTED BA network 
%%% N total number of nodes at the end
%%% m0 the initial number of nodes
%%% m links for each new node
%%% p the probability that drives the copying attachement 
%%% alfa the probability to generate a link from the new node from the
%%% existing one 

    % Initialising adjacency matrix of initial core (fully connected)
    A = ones(m0) - eye(m0);
    A = sparse(A);

    % Loop on number of nodes to be added 
    for i = 1:N-m0
    
        k_in = full(sum(A)); % Incoming degree sequence (sum over a column)
        p_in = k_in/sum(k_in); % Probability of linking to a node

        k_out = full(sum(A')); % Outgoing degree (sum over a line) 
        p_out = k_out / sum(k_out); % Probability of starting a link
    
        %%% Our strategy to simulate preferential attachement will be
        %%% to convert the vector p of probabilities into a vector of
        %%% cumulative probabilities (i.e., a vector that sums up to 1,
        %%% which is akin to the partition of the unit interval into k
        %%% sub-intervals, whose length is proportional to the
        %%% probability of connecting to nodes).
        
        p_in = cumsum(p_in);  % cumulative probability of linking to a node 
        p_out = cumsum(p_out); % cumulative prob of starting a link 

        r = rand(m,1); % m random numbers in [0,1]
    
        %%% Adding new nodes and links
        ind1 = [];
        ind2 = [];
    
        %%% Loop on new links formed by new node: the links are formed
        %%% with the preexisting nodes whose corresponding
        %%% sub-interval in the vector p_in contains the random numbers in
        %%% r
        %%% we use p_in since we want to introduce links
        for j = 1:m
            % chose the direction of the link to add
            if rand < alfa
                % randomly select a node, weighted on his incoming degree        
                aux1 = p_in - r(j);
                aux1 = find(aux1 > 0);
                ind1 = [ind1; aux1(1)];
            else 
                %randomly select a node, weigthed on his outgoing degree
                aux2 = p_out - r(j);
                aux2 = find(aux2 > 0);
                ind2 = [ind2; aux2(1)];
            end
    
        end
    
        ind1 = unique(ind1);
        ind2 = unique(ind2);
        
        %%% Creating new rows and columns in adjacency matrix
        A = [A; zeros(1,size(A,2))];
        A = [A zeros(size(A,1),1)];

        %%% FOR THE LINK FROM EXISTING TO NEW
        A(ind2, end) = 1;

        %%% FOR THE LINK FROM NEW TO EXISTING 
        if length(ind1) > 0 
            n_rand = rand(length(ind1),1); % generate a rand num for each new link
        
            for k = 1:length(n_rand) % loop sui 3 link da aggiungere 
        
                if n_rand(k) < p % collego il link al nodo scelto con prob p
                    A(end, ind1(k)) = 1;
        
                else
                    % find the nodes tho whom the selected one is connected
                    f = find(A(ind1(k), :) == 1);
                    % se il nodo scelto non ha link in uscita, connetti al
                    % nodo stesso
                    if length(f) == 0
                        A(end, ind1(k)) = 1;
                    else
                         % select a random node among these
                        index = randsample(f,1);
                        % copy the link
                        A(end, index) = 1;
                    end
                end
            end  
        end
    end     
end