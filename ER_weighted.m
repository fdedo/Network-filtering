clear all
close all


%%% GENERATE THE NETWORK %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N = 10000; % Number of nodes
p = 0.01; % Connectivity parameter

ER = ER_network(N,p); % Generating an ER network's adjacency matrix

disp("end network generation")

%%% RANDOMLY ASSIGN THE WEIGHTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

W = sparse(N,N);  %initialize weights matrix

mu = 10;            % power law exponent
sigma = 1;           % power law coefficient 

L = nnz(ER);  % count the number of links

% randomly extract the weights out of a gaussian distribution 
rand = normrnd(mu, sigma, L,1); 

aux = find(ER>0);

W(aux) = rand;

disp("end weight generation")

save('saveweightedER.mat', "W")

disp("end saving network")


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

k = full(sum(ER)); % Computing the degree sequence



%%% VISUALIZING THE NETWORK

G = graph(ER); %first create the graph object from the adiacience matric #F5D9A4 #AC291E

plot(G, 'Layout','auto', 'EdgeColor','white', 'NodeColor','#F5D9A4') 
%the layout specifics how to dispose the nodes
set(gca,'Color','k')




