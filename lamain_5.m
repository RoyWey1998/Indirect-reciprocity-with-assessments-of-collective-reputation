% main code for reputation evolving process with 3-lambda and ALLC ALLD version
tic
load('Gamestc10_5.mat');
load('Popstc30_5.mat');
ladefault_5;
err = 0.05; % error rate
strNum = 1; % select a strategy number
str = strategy_list(strNum,:); % type of strategy
lambda_list = [0.1,0.4,0.7]; % value list of assessment criterion
allCR = zeros(nPopstc,3);
CRspt1 = zeros(nPopstc,ngs);
CRspt2 = zeros(nPopstc,ngs);
CRspt3 = zeros(nPopstc,ngs); % init cooperation rate statistic for 3 types of Li players
coreNum = 12;
par = parpool('local',coreNum);
parfor nPs = 1:nPopstc
    if mod(nPs,4000)==1
        disp(nPs);
    end
    popstc = Popstc30_5(nPs,:); % pick one population structure
    % run for cooperation rate (whole as well as separate) under this
    % population structure
    [CR,crspt1,crspt2,crspt3] = larepevo2cr_5(str,popstc,N,ngs,repItnum,kmin,kmax,Gamestc10_5,q,err,lambda_list);
    allCR(nPs,:) = CR;
    CRspt1(nPs,:) = crspt1; CRspt2(nPs,:) = crspt2; CRspt3(nPs,:) = crspt3;
end
delete(par)
filename = strcat();
save(filename,'allCR','CRspt1','CRspt2','CRspt3');
toc