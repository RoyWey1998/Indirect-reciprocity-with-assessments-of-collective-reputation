tic
load('Gamestc10_3.mat');
load('Popstc60_3.mat');
default;
err = 0.05; % error rate
strNum = 1;
str = strategy_list(strNum,:); % type of strategy
coreNum = 12;
par = parpool('local',coreNum);
for l = 1:9 % parfor endpoint must be internal
    lambda = 0.1*l; % assessment criteria
    disp(lambda); % monitoring flag
    LCR = zeros(nPopstc,1);
    LCRspt = zeros(nPopstc,ngs); %init
    parfor nPs = 1:nPopstc
        Lnum = Popstc60(nPs,1); Cnum = Popstc60(nPs,2); Dnum = Popstc60(nPs,3); % numbers of players with 3 strategies
        [Lcr,Lcrspt] = repevo2cr(str,Lnum,Cnum,Dnum,N,repItnum,kmin,kmax,Gamestc10,lambda,q,err);
        LCR(nPs) = Lcr;
        LCRspt(nPs,:) = Lcrspt;
    end
    filename = strcat();
    parsave(filename,LCR,LCRspt);
end
delete(par)
toc
