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
The other files are related to filters:
  comparison_function.m
    containes the implementation of Jaccard similarity and salience 
  
  comparison_plot.m
  
  disp_filter.m
  
  hypergeometric_filter.m
  
  polya_filter.m
  
  main_network_analysis.m
  
  randraw.m
  
  write_map.m 


