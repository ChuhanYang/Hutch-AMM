function [Tr,m_list] = AMM_coarse_uni_ver2(A,B,m,index_p,replacement)
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
if nargin < 5
    replacement = true;
end
m_list = randsample(k,m,replacement);
for j = 1:1:m
    t = m_list(j);
    Tr = Tr + (A(:,index_p{t})*B(index_p{t},:))/p(t);
end
Tr = Tr/m;