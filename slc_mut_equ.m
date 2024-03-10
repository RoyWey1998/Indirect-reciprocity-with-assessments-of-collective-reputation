%% selection-mutation equilibrium 
% selection-mutation equilibrium for 3 strategies with significant
% mutations
function[PopCoop,AvFr] = slc_mut_equ(N,nPopstc,mu,s,Popstc,SMP,LCR)
W = zeros(nPopstc,nPopstc); % init the transition matrix
% calculate the transition probability
for nPs = 1:nPopstc
    Lnum = Popstc(nPs,1); Cnum = Popstc(nPs,2); Dnum = Popstc(nPs,3); % select a population structure
    % calculate transition probabilities, for e.g. wL2C indicating the probability
    % that a Li player transit to ALLC
    wL2C = Lnum/N * (mu/2+(1-mu)*Cnum/(N-1)/(1+exp(-s*(SMP(nPs,2)-SMP(nPs,1)))));
    wL2D = Lnum/N * (mu/2+(1-mu)*Dnum/(N-1)/(1+exp(-s*(SMP(nPs,3)-SMP(nPs,1)))));
    wC2L = Cnum/N * (mu/2+(1-mu)*Lnum/(N-1)/(1+exp(-s*(SMP(nPs,1)-SMP(nPs,2)))));
    wC2D = Cnum/N * (mu/2+(1-mu)*Dnum/(N-1)/(1+exp(-s*(SMP(nPs,3)-SMP(nPs,2)))));
    wD2L = Dnum/N * (mu/2+(1-mu)*Lnum/(N-1)/(1+exp(-s*(SMP(nPs,1)-SMP(nPs,3)))));
    wD2C = Dnum/N * (mu/2+(1-mu)*Cnum/(N-1)/(1+exp(-s*(SMP(nPs,2)-SMP(nPs,3)))));
    % fit the results into right entries of the matrix
    if Lnum>0 % Li players can be reduced
        W(nPs,Popstc(:,1) == Lnum-1 & Popstc(:,2) == Cnum+1) = wL2C; % become ALLC
        W(nPs,Popstc(:,1) == Lnum-1 & Popstc(:,2) == Cnum) = wL2D; % become ALLD       
    end
    if Cnum>0 % ALLC players can be reduced
        W(nPs,Popstc(:,1) == Lnum+1 & Popstc(:,2) == Cnum-1) = wC2L; % become Li
        W(nPs,Popstc(:,1) == Lnum & Popstc(:,2) == Cnum-1) = wC2D; % become ALLD       
    end
    if Dnum>0 % ALLD players can be reduced
        W(nPs,Popstc(:,1) == Lnum+1 & Popstc(:,2) == Cnum) = wD2L; % become Li
        W(nPs,Popstc(:,1) == Lnum & Popstc(:,2) == Cnum+1) = wD2C; % become ALLC       
    end
    W(nPs,nPs) = 1-(wL2C+wL2D+wC2L+wC2D+wD2L+wD2C); % the remaining probability representing no change    
end
v = null(W'-eye(nPopstc)); % left eigenvector corresponding to 1
if size(v,2) == 1 % if the dimension of the solution space is 1
    %nv = v./sum(v); 
    SMEq = v/sum(v); % normalization to get the selection-mutation equilibrium
    %LiCoop = LCR(:,1)'*SMEq; % average cooperation rate of Li players
    PopCoop = (Popstc(:,1)'.*LCR(:,1)'+Popstc(:,2)')/N*SMEq; % average cooperation rate of the population
    AvFr = Popstc'*SMEq/N; % average frequency of each strategy
else
    disp('error');
end
end