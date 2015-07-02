% Visualises a dat/txt file in sparse matrix format, so we can get a feel
% for the matrix

load ./Benchmark_result/network.dat
maxnodes = max(network(:, 2));
AM=zeros(maxnodes,maxnodes);

for i = 1:size(network)
    AM(network(i, 1), network(i, 2)) = network(i, 3);
end
[sortedAM, I] = sort(sum(AM>0), 'descend');

subplot(1, 2, 1);
imagesc(AM(I, I)>0); title('Benchmark Network');

load ./Actual_data/NNMthr_mean.mat

subplot(122)
[sortedNNMthr, I] = sort(sum(NNMthr_mean_not5>0), 'descend');

imagesc(NNMthr_mean_not5(I, I)>0); title('Brain Network');