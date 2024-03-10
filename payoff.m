function[pof] = payoff(k,ncr,alpha,beta,c,act)
R = alpha*(k^beta); % synergy factor
temp = R*c*ncr/k;
if act == 1
    pof = temp-c;
elseif act == 0
    pof = temp;
end
end