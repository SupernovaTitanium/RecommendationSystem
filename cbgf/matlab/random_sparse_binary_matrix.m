function Z = random_sparse_binary_matrix(n, k, row_nnz)
% Generate a random n x k sparse binary matrix with fixed row nnz.

I = [];
J = [];
V = [];
for i = 1:n
  I = [I; i * ones(row_nnz,1)];
  J = [J; randperm(k, row_nnz)']; % use randperm to get unique samples
  V = [V; ones(row_nnz,1)];
end
Z = sparse(I, J, V, n, k);