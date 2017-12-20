% NMI on real datasets
clear;

% Common config for algorithms
config.graph_filename = '/Users/xunzheng/workspace/data/com-dblp.ungraph.txt';
config.ground_truth_filename = '/Users/xunzheng/workspace/data/com-dblp.all.cmty.txt';
config.cpu = 4;
config.max_epoch = 500;
config.weight = 0.0;

% Ground-truth community
G = load_graph(config.graph_filename);
Z = load_affiliation(config.ground_truth_filename, G);
fprintf('Z: %d x %d matrix with nnz/row = %g\n', size(Z,1), size(Z,2), nnz(Z)/size(Z,1));

% Repeated experiments for a fixed dataset
candidate_k = [13000 14000];
num_repeat = 1;
nmi_bernoulli = zeros(num_repeat, length(candidate_k));
nmi_uniform = zeros(num_repeat, length(candidate_k));
nmi_bigclam = zeros(num_repeat, length(candidate_k));
for r = 1:num_repeat
  for c = 1:length(candidate_k)
    config.k = candidate_k(c);
    fprintf('Repeat %d, candidate_k = %d\n', r, config.k);
    
    % Train bernoulli
    config.alg = 'bulk';
    cols_filename = train_metropolis(config);
    Z_bernoulli = load_affiliation(cols_filename, G);
    nmi_bernoulli(r,c) = normalized_mutual_information(Z, Z_bernoulli);
    
    % Train uniform
    config.alg = 'uniform';
    cols_filename = train_metropolis(config);
    Z_uniform = load_affiliation(cols_filename, G);
    nmi_uniform(r,c) = normalized_mutual_information(Z, Z_uniform);
    
    % Train bigclam
    cols_filename = train_bigclam(config);
    Z_bigclam = load_affiliation(cols_filename, G);
    nmi_bigclam(r,c) = normalized_mutual_information(Z, Z_bigclam);
  end
end

nmi_avg = [mean(nmi_bernoulli);
           mean(nmi_uniform);
           mean(nmi_bigclam)];
nmi_std = [std(nmi_bernoulli);
           std(nmi_uniform);
           std(nmi_bigclam)];

nmi_avg, nmi_std


