%% load data
load ../A_samples

%% run kmeans
tic; 
idx = kmeans(full(A_100), 10, 'display', 'iter', 'Replicates', 1);
toc

%% ari
disp(['ARI=',num2str(rand_index(idx, y_100, 'adjusted'))]);
