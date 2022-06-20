function X=Hutch_AMM(A, B, s, h, block_sz)
% Block Hutch AMM with simple Hutchinson estimated sampling probability.
% Inputs:
%   A,B: matrix factor, their multiplication is to be estimated
%   s: sampling number for AMM
%   h: mat-vec number for simple Hutchinson
%   block_sz: average block size for each group, must be divisor


    % group splitting step
    k = size(A,2);
    num_group = ceil(k/block_sz);
    index_p = cell(1,num_group);
    W = cell(1,num_group);
    for i = 1:1:num_group
        index_p{i} =  ((i-1)*block_sz+1):(i*block_sz) ; %natural sequential
        W{i} = A(:,index_p{i})*B(index_p{i},:);
    end
    % AMM step
    p = zeros(1,num_group); 
    for j = 1:1:num_group
        p(j) = sqrt(simple_hutchinson(transpose(W{j})*W{j},h));
    end
    p = p/sum(p);
    hutchIndSample_list = randsample(num_group,s,true,p);
    Tr = zeros(size(A,1),size(B,2));
    for j = 1:1:s
        t = hutchIndSample_list(j);                                  
        Tr = Tr + W{t}/p(t); 
    end
    X = Tr/s;
