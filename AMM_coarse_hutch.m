function [Tr,m_list] = AMM_coarse_hutch(A,B,m,c,index_p)
% function AMM (Approximate Matrix Multiplication) but uses hucthplusplus/simple_hutchinson
% to estimate sampling probability. AMM_toep_coarse Samples column groups/corase partion 
% instead of fineset partition 
%  Inputs:
%   A,B: target matrix AB should be estimated
%   m: number of sampling/outer product in AMM
%   c: hyperparameter from hucthplusplus/simple_hutchinson
%   index_p: cell arrays of index partition, eg. {[1,3,5],2,[2]}

R_tran = A;
R = B;
k = length(index_p);
Tr = zeros(size(A,1),size(B,2));
% Calculate probability p
p = zeros(k,1);
W = cell(1,k); % create a cell array to store Submatrix group Product
for j = 1:1:k
    W{j} = R_tran(:,index_p{j})*R(index_p{j},:);
    p(j) = simple_hutchinson(transpose(W{j})*W{j},c);
    p(j) = sqrt(p(j));
end
p = p/sum(p);
m_list = randsample(k,m,true,p);
for j = 1:1:m
    t = m_list(j);
    Tr = Tr + W{t}/p(t);
end
Tr = Tr/m;