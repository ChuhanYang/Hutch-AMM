function [Tr,m_list] = AMM_coarse_uni(A,B,m,index_p)
% function AMM (Approximate Matrix Multiplication) but uses hucthplusplus/simple_hutchinson
% to estimate sampling probability. AMM_toep_uni Samples column groups/corase partion  
% by uniform sampling probability instead of fineset partition 
%  Inputs:
%   A,B: target matrix AB should be estimated
%   m: number of sampling/outer product in AMM
%   index_p: cell arrays of index partition, eg. {[1,3,5],2,[2]}

k = length(index_p);
Tr = zeros(size(A,1),size(B,2));
% Calculate probability p
p = ones(k,1);
p = p/sum(p);
W = cell(1,k); % create a cell array to store Submatrix group Product
for j = 1:1:k
    W{j} = A(:,index_p{j})*B(index_p{j},:);
end
m_list = randsample(k,m,true,p);
for j = 1:1:m
    t = m_list(j);
    Tr = Tr + W{t}/p(t);
end
Tr = Tr/m;