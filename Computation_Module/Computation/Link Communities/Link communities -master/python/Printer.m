load network.dat

fid = fopen('unweighted.txt', 'wt');
fprintf(fid, '%d\t%d\n', network(:, 1:2)');
fclose(fid);
