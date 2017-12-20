% NMI on synthetic datasets
clear;

% Size of random matrix Z
n = 1000;
k = 10;
row_nnz = 2;

% Common config for algorithms
config.cpu = 1;
config.max_epoch = 100;
config.weight = 1.0; % no balancing for synthetic data

% Prepare synthetic graph
Z = random_sparse_binary_matrix(n, k, row_nnz);
Z_dense = full(Z);
save -ascii 'Z.txt' Z_dense
G = build_graph(Z);
config.graph_filename = write_edgelist(G);
fprintf('Z: %d x %d matrix with nnz/row = %d\n', n, k, row_nnz);

% Repeated experiments for a fixed dataset
candidate_k = [40 60 80 100 120 140];
num_repeat = 2;
nmi_bernoulli = zeros(num_repeat, length(candidate_k));
nmi_uniform = zeros(num_repeat, length(candidate_k));
nmi_bigclam = zeros(num_repeat, length(candidate_k));
for r = 1:num_repeat
  for c = 1:length(candidate_k)
    config.k = candidate_k(c);
    fprintf('Repeat %d, candidate_k = %d\n', r, config.k);
    
    % Train bernoulli
    config.alg = 'bernoulli';
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



