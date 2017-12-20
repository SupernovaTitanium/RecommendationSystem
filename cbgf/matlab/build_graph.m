function G = build_graph(Z)
% Create an undirected graph using ZZ' as adjacency matrix. 
% Node names are just the string literals of node id's. 
% E.g., name of node 42 is '42'.
% Since matlab is 1-indexed, node id's and node names are all 1-indexed. 
% 
% Input: n x k binary sparse matrix
% Output: undirected graph

ids = (1:size(Z,1))';
names = num2str(ids);
pretty_names = strtrim(cellstr(names));
G = graph(triu(Z * Z') > 0, pretty_names, 'upper');