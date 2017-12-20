% Simulated annealing for symmetric binary matrix factorization:
%
% max_{Z} f(Z) = \sum_{i \le j, (i,j) \in E} I(z_i^T z_j > 0) / (2*pos)
%              + \sum_{i \le j, (i,j) \notin E} I(z_i^T z_j = 0) / (2*neg)
%
% Z: n x k binary matrix
% z_i: i-th row of Z
% n: number of nodes in the graph
% E: edge set of the graph, including all self-loops.
% pos: number of 1's in upper triangle of adjacency matrix
% neg: number of 0's in upper triangle of adjacency matrix
%
% The objective is weighted to address label imbalance. 
% A perfect recovery will produce objective value of 1.
%
% Algorithm: simulated annealing via Metropolis-Gibbs.
% SA: at epoch m, draw samples from exp(f(Z)/T_m), T_m is the temperature. 
% Gibbs: direct sampling is hard, so draw from conditionals exp(h(z_i)/T_m).
% Metropolis: even conditionals are hard, so use Bernoulli(mu_i) as proposal,
% where mu_i is the average of node i's neighbors.

clear;

% Params
filename = '~/workspace/data/ca-AstroPh.txt';
k = 1000;
max_epoch = 20;

% Create graph
G = load_graph(filename);
n = G.numnodes;
A = (G.adjacency + speye(n,n)) > 0; % add self-loops
pos = nnz(triu(A)); % num of 1's in upper triangle
neg = n*(n+1)/2 - pos; % num of 0's in upper triangle
obj = @(Ahat) nnz(Ahat(triu(A==1)) == 1)/pos/2 ...
            + nnz(Ahat(triu(A==0)) == 0)/neg/2;

% Randomly init one assignment for each row of Z
Z = sparse(1:n, randi(k,n,1), ones(n,1), n, k);

% Gibbs chain
fprintf('%6s%12s%12s%12s%12s%12s%12s%12s%12s\n', ...
  'epoch', 'obj', 'T', 'delta', 'nnz(Z)/n', 'nnz(Z)/k', 'mh_rate', 'sec/epoch', 'elapsed');
fprintf('%6d%12.4f%12.4f%12.4f%12.4f%12.4f%12.4f%12.4f%12.1f\n', ...
  0, obj(Z * Z' > 0), 0, 0, nnz(Z)/n, nnz(Z)/k, 0, 0, 0);
elapsed = 0;
for epoch = 1:max_epoch
  tic;
  mh_count = 0;
  delta = 0.001 / (1 + epoch); % mu might contain 0 or 1, so add a little bit of slack
  T = 1 / log(1 + epoch); % temperature
  for i = randperm(n)
    % Set up
    z_old = Z(i,:)';
    nbr = G.neighbors(i);
    d = G.degree(i);
    y = sparse([nbr; i], ones(d+1,1), ones(d+1,1), n, 1); % ground-truth for node i
    % Proposal
    mu = sum(Z(nbr,:),1)' / d; % average of neighbors
    mu(mu == 1) = 1 - delta; % avoid numerical problems
    mu(mu == 0) = delta;
    z_new = rand(k, 1) < mu; % draw from Bernoulli(mu)
    theta = log(mu ./ (1-mu)); % natural parameter of Bernoulli(mu)
    c = theta' * (z_old - z_new); % correction term for proposal
    % Accept
    h = @(pred) nnz(pred(y==1)==1)*neg/pos + nnz(pred(y==0)==0);
    a = (h(Z * z_new > 0) - h(Z * z_old > 0)) / T + c;
    if log(rand()) < a
      Z(i,:) = z_new';
      mh_count = mh_count + 1;
    end
  end
  t = toc;
  elapsed = elapsed + t;
  fprintf('%6d%12.4f%12.4g%12.4g%12.4f%12.4f%12.4f%12.4f%12.1f\n', ...
    epoch, obj(Z * Z' > 0), T, delta, ...
    nnz(Z)/n, nnz(Z)/k, mh_count/n, t, elapsed);
end
