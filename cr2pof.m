%% from cooperation rate to payoff 
function[SMP] = cr2pof(Popstc,Gamestc,LCRspt,alpha,beta,c,pk,kmin,Lstcprob,Cstcprob,Dstcprob)
nPopstc = size(Popstc,1); % number of population structures
nGamestc = size(Gamestc,1); % number of game states
SMP = zeros(nPopstc,3); % init
LCRspt(isnan(LCRspt)) = 0; % replace invalid value with 0
k = sum(Gamestc,2)'; % vector of numbers of gamers under each game state
k1 = Gamestc(:,1); k2 = Gamestc(:,2); % vector of numbers of Li and ALLC gamers
for nPs = 1:nPopstc    
    cr = LCRspt(nPs,:)'; % cooperation rate of Li under each population structure
    pi = zeros(nGamestc,4); % init payoff statistic, 4 columns for Li cooperate, defect, ALLC, and ALLD
    for nGs = 1:nGamestc
        lnc = 0:k1(nGs);% list of possible numbers of cooperating Li gamers
        bino1 = binopdf(lnc(1:k1(nGs)),k1(nGs)-1,cr(nGs)); % length k1
        bino2 = binopdf(lnc,k1(nGs),cr(nGs)); % length k1+1
        payoff1 = payoff(k(nGs),lnc+k2(nGs)+1,alpha,beta,c,1); 
        % for Li players only,actually only payoff1(1:k1(nGs)) here are
        % useful, which corresponding to k2+1:k1+k2 cooperators
        payoff2 = payoff(k(nGs),lnc+k2(nGs),alpha,beta,c,1); 
        % for Li and ALLC players, length k1+1
        pi(nGs,1) = bino1*payoff1(1:k1(nGs))'; % cooperation payoff for Li player 
        pi(nGs,2) = bino1*(payoff2(1:k1(nGs))+c)'; % defection payoff for Li player 
        pi(nGs,3) = bino2*payoff2'; % payoff for ALLC
        pi(nGs,4) = bino2*(payoff2+c)'; %payoff for ALLD
    end
    SMP(nPs,1) = (pk(k-kmin+1).*Lstcprob(nPs,:))*(cr.*pi(:,1)+(1-cr).*pi(:,2));
    SMP(nPs,2) = (pk(k-kmin+1).*Cstcprob(nPs,:))*pi(:,3);
    SMP(nPs,3) = (pk(k-kmin+1).*Dstcprob(nPs,:))*pi(:,4);
end
end