function G = load_graph(filename)
% Load undirected graph from SNAP-style edge list file.
% Treat lines starting with # as comments.
% Remove duplicated edges, e.g. (i,j) and (j,i).
% Some of the nodes may have self-loops, but removing them from the edge 
% list will remove some isolated nodes too, so remove self-loops after the 
% graph is constructed. 
% 
% Input: edge list file
% Output: undirected graph, no self loops.

% Read edgelist file from SNAP
fid = fopen(filename);
cell_content = textscan(fid, '%s %s', 'CommentStyle', '#', 'CollectOutput', true);
str_content = string(cell_content{:});
edgelist = cellstr(unique(sort(str_content,2), 'rows'));

% Create graph
G = graph(edgelist(:,1), edgelist(:,2));

% Remove self loops
n = G.numnodes;
selfloops = nonzeros(G.findedge(1:n, 1:n));
G = G.rmedge(selfloops);

fclose(fid);