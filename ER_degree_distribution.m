clear all
close all

%%% In the following we check that for large N the degree distribution of
%%% the ER model is approximately Gaussian

N = 2e3; % Number of nodes
p = 0.2; % Connectivity parameter
Nsims = 100; % Number of adjacency matrices of the ER network to generate

deg = []; % Array to store all degree sequences

for ns = 1:Nsims

    A = ER_network(N,p); % Generating one instance of an ER network's adjacency matrix

    k = full(sum(A)); % Computing the degree sequence : sum() sums over the columns, while full create an array of dimension N filling with 0 the elements missing 

    deg = [deg; k']; % Storing the degree sequence (it concatenates the simulations in adiacent columns)

end

%%% Calculating the mean degree for each node trough the simualtion: 
%%% (I ADDED THIS BE CAREFUL)
kmin = mean(deg, 2);


%%% Plotting histogram vs Gaussian distribution

histogram(k,50,'Normalization','pdf')

hold on

m = N*p;
s = sqrt(N*p);
x = linspace(m-5*s,m+5*s,1000);
y = exp(-(x-m).^2/(2*s^2))/sqrt(2*pi*s^2);

plot(x,y,'r','LineWidth',2)

xlabel('$k$','Interpreter','latex')
ylabel('$p(k)$','Interpreter','latex')
set(gca,'FontSize',20)
