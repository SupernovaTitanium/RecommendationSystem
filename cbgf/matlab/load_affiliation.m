function Z = load_affiliation(filename, graph)
% Load node affiliation from file. Skip empty line.
% 
% Input: 
% Affiliation file: contains k lines representing k communities. 
% Each line contains *name* of the nodes that are affiliated to the 
% corresponding community.
% Graph: graph object for translating node names to node id's.
% Output: 
% Z: a n x k sparse binary matrix, where each column is equivalent to each
% line in the affiliation file. 

fid = fopen(filename);
I = [];
J = [];
S = [];
k = 0;
while true
  line = fgetl(fid); % get line without '\n'
  if ~ischar(line) % test EOF
    break
  end
  line = strtrim(line);
  if isempty(line) % skip emtpy line
    continue
  end
  k = k + 1;
  names = strsplit(line)';
  I = [I; names];
  J = [J; k * ones(size(names))];
  S = [S; ones(size(names))];
end
fclose(fid);
Z = sparse(graph.findnode(I), J, S, graph.numnodes, k);