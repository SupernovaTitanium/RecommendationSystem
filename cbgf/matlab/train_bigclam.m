function cols_filename = train_bigclam(config)
% Run bigclam on graph with k communities. 
% 
% Input: config containing graph_filename and k.
% Output: sparse binary matrix filename.

% Path for the binary executable
exe = '/Users/xunzheng/workspace/snap/examples/bigclam/bigclam';

% Caution! This is the default in bigclam.  
cols_filename = 'cmtyvv.txt';

% cmd = sprintf('taskset -c 1,3,5,7,9,11,13,15 %s -i:%s -c:%d -nt:1', bin_filename, graph_filename, k);
cmd = sprintf('%s -i:%s -c:%d -nt:%d', exe, config.graph_filename, config.k, config.cpu);
[status, output] = system(cmd);
% status = system(cmd);
if status ~= 0
  error('status is not 0!');
end