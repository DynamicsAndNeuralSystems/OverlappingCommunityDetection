%load benchmark network data
load('Computation_Module/Conversions/networks/network56.dat')

%Run OCDAs on benchmark
Computation(network56, {'OSLOM','Jerry','Shen','NNMF','Infomap'}, 1, 'networks/community56')

%Visualise performance
% cd ../
% cd Visualisation_Module
Visualisation({'Benchmark', 'OSLOM', 'Shen', 'NNMF', 'Infomap','SLPA'})
