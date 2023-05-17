# Network-filtering
Repository containing all the codes implemented in my thesis work
------------------------------------------------------------------
The codes starting with BA are related to Barabasi-Albert network: 
  BA_degree_distribution.m
    plot the degree distribution of a BA net
    
  BA_directed.m
    it containes the function to generate a directed BA net
  
  BA_network.m
     it containes the function to generate an undirected BA net
  
  BA_weighted.m
   it containes the function to generate a directed BA net with randomly assigned weights to each link
  
------------------------------------------------------------------  
The codes starting with ER are related to Erdos-Renyi network:
  ER_degree_distribution.m
     it plots the defree distribution for an ER net, checking that it tends to be gaussian for large nets
     
  ER_network.m
     it containes the function to generate an undirected ER net
  
  ER_weighted.m 
    it containes the function to generate an undirected ER net with randomly assigned weights to each link

------------------------------------------------------------------
There are three files related to filters:
  disp_filter.m
    containes the function that implement the disparity filter 
  
  hypergeometric_filter.m
    containes the function that implement the hypergeometric filter
  
  polya_filter.m
    containes the function that implement the polya filter
  
------------------------------------------------------------------
An in the end there are files related to the comparison of filters:
  comparison_function.m
    containes the implementation of Jaccard similarity and salience 
  
  comparison_plot.m
    This file aims to plot measurement of the different performance of the
    3 filters implemented in function of the multivariate significance level
  
  main_network_analysis.m
    This file can be used to quickly extract a backbone out of the desired network

------------------------------------------------------------------
randraw.m 
  is an imported function that extract numbers out of particular distribution 
  (used to extract from a power law distribution the weights of the synthetic networks)



