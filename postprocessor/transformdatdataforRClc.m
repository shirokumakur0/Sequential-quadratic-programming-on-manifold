clear;
datset = [
    2,3;
    2,4;
    3,5;
    3,6;
    4,6;
    4,7;
    4,8;
];

tolKKTrespowerset = [2, 3, 4, 5, 6, 7, 8]; % 1e-* tolerance


[N,~] = size(datset);
for i = 1:N
    rdim = datset(i,1);
    cdim = datset(i,2);
    filename = sprintf('with_SQP_zz_RC_lc_RDim%dCDim%d.dat', rdim, cdim);
    M = readmatrix(filename);
    [NM,~] = size(M);
    
    for j = 1:NM
       outputdata = M(j,:);
       tolKKTres = M(j,17);
       transformedfilename = sprintf('with_SQP_zz_RC_lc_RDim%dCDim%dTol%d.dat', rdim, cdim, tolKKTres);
       dlmwrite(transformedfilename, outputdata, 'delimiter', ',', 'precision', 16, '-append');
    end 
end