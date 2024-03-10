%% from cooperation rate to payoff in assessment criterion evolution
function[SMP] = lacr2pof(Popstc,Gamestc,CRspt1,CRspt2,CRspt3,alpha,beta,c,pk,kmin,Lstcprob1,Lstcprob2,Lstcprob3)
nPopstc = size(Popstc,1); % number of population structures
nGamestc = size(Gamestc,1); % number of game states
SMP = zeros(nPopstc,3); %init
CRspt1(isnan(CRspt1)) = 0;CRspt2(isnan(CRspt2)) = 0;CRspt3(isnan(CRspt3)) = 0; % replace invalid value with 0
k = sum(Gamestc,2)'; % vector of numbers of gamers under each game state
k1 = Gamestc(:,1); k2 = Gamestc(:,2); k3 = Gamestc(:,3);% vector of numbers of 3 kinds of gamers
for nPs = 1:nPopstc
    cr = [CRspt1(nPs,:)',CRspt2(nPs,:)',CRspt3(nPs,:)']; % cooperation rate of 3 kinds gamers under each population structure
    pi = zeros(nGamestc,6); % init payoff statistic, 6 columns for C and D of each type
    for nGs = 1:nGamestc
        lnc1 = 0:k1(nGs);lnc2 = 0:k2(nGs);lnc3 = 0:k3(nGs); % lists of possible numbers of cooperating gamers
        % number of cooperators and corresponding probabilities without focal player
        binoex1 = binopdf(lnc1(1:k1(nGs)),k1(nGs)-1,cr(nGs,1)); % length k1
        binoex2 = binopdf(lnc2(1:k2(nGs)),k2(nGs)-1,cr(nGs,2)); % length k2
        binoex3 = binopdf(lnc3(1:k3(nGs)),k3(nGs)-1,cr(nGs,3)); % length k3
        % number of cooperators and corresponding probabilities with focal player
        binoin1 = binopdf(lnc1,k1(nGs),cr(nGs,1)); % length k1+1
        binoin2 = binopdf(lnc2,k2(nGs),cr(nGs,2)); % length k2+1
        binoin3 = binopdf(lnc3,k3(nGs),cr(nGs,3)); % length k3+1
        payoffc = payoff(k(nGs),1:k(nGs),alpha,beta,c,1); % cooperating payoff (with 1 to k cooperators), length k
        payoffd = [0,payoffc(1:end-1)+c]; % defecting payoff (with 0 to k-1 cooperators), length k
        % calculate expect payoff of 3 types upon cooperate and defect
        % type 1
        for j1 = 0:k1(nGs)-1
            for j2 = 0:k2(nGs)
                for j3 = 0:k3(nGs)
                    pi(nGs,1) = pi(nGs,1)+binoex1(j1+1)*binoin2(j2+1)*binoin3(j3+1)*payoffc(j1+j2+j3+1);
                    pi(nGs,2) = pi(nGs,2)+binoex1(j1+1)*binoin2(j2+1)*binoin3(j3+1)*payoffd(j1+j2+j3+1);
                end
            end
        end
        % type 2
        for j1 = 0:k2(nGs)-1
            for j2 = 0:k1(nGs)
                for j3 = 0:k3(nGs)
                    pi(nGs,3) = pi(nGs,3)+binoex2(j1+1)*binoin1(j2+1)*binoin3(j3+1)*payoffc(j1+j2+j3+1);
                    pi(nGs,4) = pi(nGs,4)+binoex2(j1+1)*binoin1(j2+1)*binoin3(j3+1)*payoffd(j1+j2+j3+1);
                end
            end
        end
        % type 3
        for j1 = 0:k3(nGs)-1
            for j2 = 0:k1(nGs)
                for j3 = 0:k2(nGs)
                    pi(nGs,5) = pi(nGs,5)+binoex3(j1+1)*binoin1(j2+1)*binoin2(j3+1)*payoffc(j1+j2+j3+1);
                    pi(nGs,6) = pi(nGs,6)+binoex3(j1+1)*binoin1(j2+1)*binoin2(j3+1)*payoffd(j1+j2+j3+1);
                end
            end
        end
    end
    % calculate the final expect SMP
    SMP(nPs,1) = (pk(k-kmin+1).*Lstcprob1(nPs,:))*(cr(:,1).*pi(:,1)+(1-cr(:,1)).*pi(:,2));
    SMP(nPs,2) = (pk(k-kmin+1).*Lstcprob2(nPs,:))*(cr(:,2).*pi(:,3)+(1-cr(:,2)).*pi(:,4));
    SMP(nPs,3) = (pk(k-kmin+1).*Lstcprob3(nPs,:))*(cr(:,3).*pi(:,5)+(1-cr(:,3)).*pi(:,6));
end
end