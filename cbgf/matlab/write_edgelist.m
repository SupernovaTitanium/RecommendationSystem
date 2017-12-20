function graph_filename = write_edgelist(G)
% Write edge list of G into a file.
%
% Input:
% G: an undirected graph
% Output:
% filename: the target edge list file.

% Caution! Using a fixed filename for the graph.
graph_filename = 'mygraph.txt';
writetable(G.Edges, graph_filename, 'WriteVariableNames', false, 'Delimiter', '\t');