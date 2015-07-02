%% Plots the difference and the similarities with a benchmark and a real brain network
% Brandon Lam, 11-07-2014

clear all; close all; clc;

%% Obtaining info - Benchmark
Ntwk2AdjMat
BenchDegree = sum(Adj_Mat~=0);

%% Obtaining info - Actual
load NNMthr_mean.mat
NNMthr = NNMthr_mean_not5;
ActualDegree = sum(NNMthr~=0);

%% Plotting the degrees
subplot(1, 3, 1);
A = 1:length(BenchDegree);
plot(A, sort(BenchDegree), 'b', A, sort(ActualDegree), 'm');
title('Degree distribution'); xlabel('Sorted Nodes'); ylabel('Degree');
legend('Benchmark', 'Actual');

%% Plotting the strengths
subplot(132);
plot(A, sort(sum(Adj_Mat)/max(sum(Adj_Mat))), 'b', A, sort(sum(NNMthr)/max(sum(NNMthr))), 'm')
title('Strength distribution'); xlabel('Sorted Nodes'); ylabel('Normalised strength');
legend('Benchmark', 'Actual');

%% Plotting the degree vs strength
subplot(133);
plot(BenchDegree, sum(Adj_Mat)/max(sum(Adj_Mat)), 'bx', ActualDegree,...
    sum(NNMthr)/max(sum(NNMthr)), 'mx')
title('Scatter plot'); xlabel('Degree'); ylabel('Normalised strength');
legend('Benchmark', 'Actual');