function cols_filename = train_metropolis(config)
% Run metropolis algortihm on graph, using provided config struct.
% 
% Input: structure for configs
% Output: sparse binary matrix filename.

% Path for the binary executable
exe = '/Users/xunzheng/workspace/cbgf/cmake-build-release/cbgf';

% Caution! Using the same filename for all experiments. 
config_filename = 'myconfig.txt';
cols_filename = 'mycols.txt';

% Write config into a file
fid = fopen(config_filename, 'w');
fprintf(fid, 'cpu %d\n', config.cpu);
fprintf(fid, 'K %d\n', config.k);
fprintf(fid, 'maxEpoch %d\n', config.max_epoch);
fprintf(fid, 'weight %f\n', config.weight);
fprintf(fid, 'algorithm %s\n', config.alg);
fprintf(fid, 'edgeListFile %s\n', config.graph_filename);
fprintf(fid, 'outputByColsFile %s\n', cols_filename);
fclose(fid);

% cmd = sprintf('taskset -c 0,2,4,6,8,10,12,14 %s %s', exe, config_filename);
cmd = sprintf('%s %s', exe, config_filename);
[status, output] = system(cmd);
% status = system(cmd); % allow ouptut on screen
if status ~= 0
  error('status is not 0!');
end