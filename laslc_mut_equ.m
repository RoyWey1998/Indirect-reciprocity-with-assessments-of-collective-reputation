%% selection-mutation equilibrium for lambda evolve
% selection-mutation equilibrium for 3-lambda assessment evolving
% with significant mutations
function[PopCoop,AvFr] = laslc_mut_equ(N,nPopstc,mu,s,Popstc,SMP,allCR)
W = zeros(nPopstc,nPopstc); % init the transition matrix
% calculate the transition probability
for nPs = 1:nPopstc
    num1 = Popstc(nPs,1); num2 = Popstc(nPs,2); num3 = Popstc(nPs,3); % select a population structure
    % calculate transition probabilities, for e.g. wa2b indicating the probability
    % that a player with lambda_1 transit to lambda_2
    wa2b = num1/N * (mu/2+(1-mu)*num2/(N-1)/(1+exp(-s*(SMP(nPs,2)-SMP(nPs,1)))));
    wa2c = num1/N * (mu/2+(1-mu)*num3/(N-1)/(1+exp(-s*(SMP(nPs,3)-SMP(nPs,1)))));
    wb2a = num2/N * (mu/2+(1-mu)*num1/(N-1)/(1+exp(-s*(SMP(nPs,1)-SMP(nPs,2)))));
    wb2c = num2/N * (mu/2+(1-mu)*num3/(N-1)/(1+exp(-s*(SMP(nPs,3)-SMP(nPs,2)))));
    wc2a = num3/N * (mu/2+(1-mu)*num1/(N-1)/(1+exp(-s*(SMP(nPs,1)-SMP(nPs,3)))));
    wc2b = num3/N * (mu/2+(1-mu)*num2/(N-1)/(1+exp(-s*(SMP(nPs,2)-SMP(nPs,3)))));
    % fit the results into right entries of the matrix
    if num1>0 % lambda_1 players can be reduced
        W(nPs,Popstc(:,1) == num1-1 & Popstc(:,2) == num2+1) = wa2b; % become lambda_2
        W(nPs,Popstc(:,1) == num1-1 & Popstc(:,2) == num2) = wa2c; % become lambda_3       
    end
    if num2>0 % lambda_2 players can be reduced
        W(nPs,Popstc(:,1) == num1+1 & Popstc(:,2) == num2-1) = wb2a; % become lambda_1
        W(nPs,Popstc(:,1) == num1 & Popstc(:,2) == num2-1) = wb2c; % become lambda_3       
    end
    if num3>0 % lambda_3 players can be reduced
        W(nPs,Popstc(:,1) == num1+1 & Popstc(:,2) == num2) = wc2a; % become lambda_1
        W(nPs,Popstc(:,1) == num1 & Popstc(:,2) == num2+1) = wc2b; % become lambda_2       
    end
    W(nPs,nPs) = 1-(wa2b+wa2c+wb2a+wb2c+wc2a+wc2b); % the remaining probability representing no change 
end
v = null(W'-eye(nPopstc)); % left eigenvector corresponding to 1
if size(v,2) == 1 % if the dimension of the solution space is 1
    SMEq = v/sum(v); % normalization to get the selection-mutation equilibrium
    PopCoop = sum(Popstc'.*allCR')/N*SMEq; % average cooperation rate of the population
    AvFr = Popstc'*SMEq/N; % average frequency of each strategy
else
    disp('error');
end
end