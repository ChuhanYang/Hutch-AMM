% code by RaphaelArkadyMeyerNYU from https://github.com/RaphaelArkadyMeyerNYU/HutchPlusPlus

function T = hutchplusplus(A,m)
    S = 2*randi(2,size(A,1),m/3);
    G = 2*randi(2,size(A,1),m/3);
    [Q,~]=qr(A*S,0);
    G = G - Q*(Q'*G);
    T = trace(Q'*A*Q)+1/size(G,2)*trace(G'*A*G);
end