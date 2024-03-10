% tool code for lambda evolution (blocked)
%% 1 payoff to selection-mutation equilibrium (\mu \to 0) for lambda
% evolution (3-value version)
load('Gamestc10_3.mat');
load('Popstc60_3.mat');
load('GamestcProb10_3.mat');

% payoff to fixation probability (for fixed alpha)
par = parpool('local',11);
Popstc = Popstc60;
Gamestc = Gamestc10;
ladefault;
pk = [kmin:kmax]/sum(kmin:kmax);
Gametype = [1,2;1,3;2,3]; % 3-strategy
nGametype = size(Gametype,1);
parfor a = 1:26
    disp(a);
    alpha = 0.1*(a-1);
    FXP = zeros(8,2*nGametype); % 8 rows for 8 strategies, 6 columns for pairwise fixation probabilities among 3 strategies
    for strNum = 1:8
        % load cooperation rate CRspt1, CRspt2, CRspt3
        Dataname = strcat();
        data = matfile(Dataname);
        CRspt1 = data.CRspt1; CRspt2 = data.CRspt2; CRspt3 = data.CRspt3;
        SMP = lacr2pof(Popstc,Gamestc,CRspt1,CRspt2,CRspt3,alpha,beta,c,pk,kmin,Lstcprob,Cstcprob,Dstcprob);
        SMPsavepath = strcat();
        if ~exist(SMPsavepath)
            mkdir(SMPsavepath);
        end            
        SMPname = strcat(SMPsavepath,'\L',num2str(strNum),'_SMP.mat');
        parsave_SMP(SMPname,SMP);
        for gt = 1:nGametype
            stg1 = Gametype(gt,1);
            stg2 = Gametype(gt,2);
            fxp = fixprob(s,stg1,stg2,Popstc,SMP);
            FXP(strNum,2*gt-1:2*gt) = fxp;
        end
    end
    FXPsavepath = strcat();
    if ~exist(FXPsavepath)
        mkdir(FXPsavepath);
    end
    filename = strcat(FXPsavepath,'\fixprob.mat');
    parsave_FXP(filename,FXP);

    % from fixation probabilities to selection-mutation equilibrium and average
    % cooperation rate (rare mutations)
    erthd = 1e-14; % calculation error threshold, with large alpha, the 
    % eigenvalue corresponding to the null space of T'-eye(3)  may not be 
    % accurately 0 due to the limitation of calculation accuracy
    allSMEq = zeros(8,3);
    allcoop = zeros(8,1);
    % load fixation probability
    data1 = matfile(filename);
    FXP = data1.FXP;
    for strNum = 1:8
        tFXP = zeros(3,3); % init temp result matrix
        % load corresponding results of cooperation rate
        CRdataname = strcat(); 
        data2 = matfile(CRdataname);
        allCR = data2.allCR;
        lcr = [allCR(1,1),allCR(61,2),allCR(1891,3)]; % cooperation rate of 3 lambda players when they fully occupy the population
        tFXP(2,1) = FXP(strNum,1); % la2 trans2 la1
        tFXP(1,2) = FXP(strNum,2); % la1 trans2 la2
        tFXP(3,1) = FXP(strNum,3); % la3 trans2 la1
        tFXP(1,3) = FXP(strNum,4); % la1 trans2 la3
        tFXP(3,2) = FXP(strNum,5); % la3 trans2 la2
        tFXP(2,3) = FXP(strNum,6); % la2 trans2 la3
        T = tFXP/2;
        % get the transition matrix
        for gt = 1:3
            T(gt,gt) = 1-sum(T(gt,:));
        end
        [x,y] = eig(T');
        errdet = abs(diag(y)-1); 
        if min(errdet)<erthd % insurence for error casued by insufficient calculating accuracy
            v = real(x(:,errdet==min(errdet)));
            SMEq=v'/sum(v);
        else 
            SMEq = [0,0,0];
            disp([strNum,alpha]);
        end
        allSMEq(strNum,:) = SMEq;
        coop = SMEq*[lcr(1);lcr(2);lcr(3)]; % average cooperation rate for the whole population
        allcoop(strNum) = coop;
    end
    smeqsavepath = strcat();
    if ~exist(smeqsavepath)
        mkdir(smeqsavepath);
    end
    smeqfilename = strcat(smeqsavepath,'\SMEq&coop.mat');
    parsave_smeq(smeqfilename,allSMEq,allcoop);
end
delete(par)

%% 2 generate the structure list of population/game state
N = 10;
Gamestc10_5 = [];
par = parpool('local',6);
parfor i1 = 0:N
    for i2 = 0:N-i1
        for i3 = 0:N-i1-i2
            for i4 = 0:N-i1-i2-i3
                i5 = N-i1-i2-i3-i4;
                v = [i5,i4,i3,i2,i1];
                Gamestc10_5 = [Gamestc10_5;v];
            end
        end
    end
end
save('Gamestc10_5.mat','Gamestc10_5');
delete(par)

%% 3 
% Find possible game states under different population structures 
% and decide the probabilities for players of each strategy taking
% part in each situation (5-strategy version)
Popstc = Popstc30_5; 
Gamestc = Gamestc10_5;
nPopstc = size(Popstc,1); % number of population structures
nGamestc = size(Gamestc,1); % number of game states
L1stcprob = zeros(nPopstc,nGamestc);
L2stcprob = zeros(nPopstc,nGamestc);
L3stcprob = zeros(nPopstc,nGamestc);
Cstcprob = zeros(nPopstc,nGamestc);
Dstcprob = zeros(nPopstc,nGamestc); % init

par = parpool('local',12);
for nPs = 1:nPopstc
    for nGs = 1:nGamestc
        if Popstc(nPs,1)>=Gamestc(nGs,1) && Popstc(nPs,2)>=Gamestc(nGs,2) && Popstc(nPs,3)>=Gamestc(nGs,3) && Popstc(nPs,4)>=Gamestc(nGs,4) && Popstc(nPs,5)>=Gamestc(nGs,5)
            n1 = Popstc(nPs,1); n2 = Popstc(nPs,2); n3 = Popstc(nPs,3); n4 = Popstc(nPs,4); n5 = Popstc(nPs,5);
            k1 = Gamestc(nGs,1); k2 = Gamestc(nGs,2); k3 = Gamestc(nGs,3); k4 = Gamestc(nGs,4); k5 = Gamestc(nGs,5);
            k = sum(Gamestc(nGs,:));
            if k1>=1
                L1stcprob(nPs,nGs) = nchoosek(n1-1,k1-1)*nchoosek(n2,k2)*nchoosek(n3,k3)*nchoosek(n4,k4)*nchoosek(n5,k5)/nchoosek(N-1,k-1); 
            end
            if k2>=1
                L2stcprob(nPs,nGs) = nchoosek(n1,k1)*nchoosek(n2-1,k2-1)*nchoosek(n3,k3)*nchoosek(n4,k4)*nchoosek(n5,k5)/nchoosek(N-1,k-1); 
            end
            if k3>=1
                L3stcprob(nPs,nGs) = nchoosek(n1,k1)*nchoosek(n2,k2)*nchoosek(n3-1,k3-1)*nchoosek(n4,k4)*nchoosek(n5,k5)/nchoosek(N-1,k-1); 
            end
            if k4>=1
                Cstcprob(nPs,nGs) = nchoosek(n1,k1)*nchoosek(n2,k2)*nchoosek(n3,k3)*nchoosek(n4-1,k4-1)*nchoosek(n5,k5)/nchoosek(N-1,k-1); 
            end
            if k5>=1
                Dstcprob(nPs,nGs) = nchoosek(n1,k1)*nchoosek(n2,k2)*nchoosek(n3,k3)*nchoosek(n4,k4)*nchoosek(n5-1,k5-1)/nchoosek(N-1,k-1); 
            end
        end
    end
end
filename = 'GamestcProb10_5.mat';
save(filename,'L1stcprob','L2stcprob','L3stcprob','Cstcprob','Dstcprob');
delete(par)

%% 4 payoff to average frequency for lambda evolution(\mu > 0)
% from payoff to average frequency and cooperation rate at 
% selection-mutation equilibrium with significant mutation rate (3-strategy version)
Popstc = Popstc60;
Gamestc = Gamestc10;
pk = [kmin:kmax]/sum(kmin:kmax);
mu = 0.05;
par = parpool('local',11);
parfor a = 1:26
    alpha = 0.1*(a-1);
    allAvFr = zeros(8,3); 
    allcoop = zeros(8,1); % init
    for strNum = 1:8
        % load SMP data
        Dataname1 = strcat();
        data1 = matfile(Dataname1);
        SMP = data1.SMP;        
        % load cooperation rate allCR
        Dataname2 = strcat();
        data2 = matfile(Dataname2);
        allCR = data2.allCR;
        allCR(isnan(allCR)) = 0;
        [PopCoop,AvFr] = laslc_mut_equ(N,nPopstc,mu,s,Popstc,SMP,allCR);
        allcoop(strNum) = PopCoop;
        allAvFr(strNum,:) = AvFr';
    end
    avfrsavepath = strcat();
    avfrfilename = strcat(avfrsavepath,'_AvFr&coop.mat');
    parsave_avfr(avfrfilename,allAvFr,allcoop);
end
delete(par)

%% 5-1 payoff to selection-mutation equilibrium (\mu \to 0) for lambda
% evolution (3-value with ALLC/ALLD version)
load('Gamestc10_5.mat');
load('Popstc30_5.mat');
load('GamestcProb10_5.mat');

% payoff to fixation probability (for fixed alpha)
par = parpool('local',8);
Popstc = Popstc30_5;
Gamestc = Gamestc10_5;
ladefault_5;
pk = [kmin:kmax]/sum(kmin:kmax);
Gametype = [1,2;1,3;1,4;1,5;2,3;2,4;2,5;3,4;3,5;4,5]; % 5 different types of strategies
nGametype = size(Gametype,1);
parfor strNum = 1:8
    disp(strNum);
    FXP = zeros(26,2*nGametype);% init
    Dataname = strcat(); % load cooperation rate
    data = matfile(Dataname);
    CRspt1 = data.CRspt1; CRspt2 = data.CRspt2; CRspt3 = data.CRspt3;
    for a = 1:26
        alpha = 0.1*(a-1);
        SMP = lacr2pof_5(Popstc,Gamestc,CRspt1,CRspt2,CRspt3,alpha,beta,c,pk,kmin,L1stcprob,L2stcprob,L3stcprob,Cstcprob,Dstcprob); % 待写
        SMPsavepath = strcat();
        if ~exist(SMPsavepath)
            mkdir(SMPsavepath);
        end            
        SMPname = strcat(SMPsavepath,'/α=',num2str(alpha),'_SMP.mat');
        parsave_SMP(SMPname,SMP);
        for gt = 1:nGametype
            stg1 = Gametype(gt,1);
            stg2 = Gametype(gt,2);
            fxp = fixprob(s,stg1,stg2,Popstc,SMP);
            FXP(a,2*gt-1:2*gt) = fxp;
        end
    end
    FXPsavepath = strcat();
    if ~exist(FXPsavepath)
        mkdir(FXPsavepath);
    end
    filename = strcat(FXPsavepath,'/fixprob.mat');
    parsave_FXP(filename,FXP);    
end
delete(par)

%% 5-1-1 FXP rearrange (once)
for strNum = 1:8
    allFXP = zeros(25,20);
    for i = 1:5
        dataname = strcat();
        data = matfile(dataname);
        FXP = data.FXP;
        allFXP(i*5-4:i*5,:) = FXP(i*5-4:i*5,:);
    end
    filename = strcat();
    save(filename,'allFXP');
end

%% 5-2
% from fixation probabilities to selection-mutation equilibrium and average
% cooperation rate (rare mutations)
erthd = 1e-14; % calculation error threshold, with large alpha, the 
% eigenvalue corresponding to the null space of T'-eye(3)  may not be 
% accurately 0 due to the limitation of calculation accuracy
Gametype = [1,2;1,3;1,4;1,5;2,3;2,4;2,5;3,4;3,5;4,5]; % 5 different types of strategies
nGametype = size(Gametype,1);
%par = parpool('local',6);
for strNum = 1:8
    disp(strNum);
    alpha = 0;
    allSMEq = zeros(26,5);
    allcoop = zeros(26,1); % init
    % load fixation probability
    FXPname = strcat();
    data1 = matfile(FXPname);
    FXP = data1.allFXP;
    % load corresponding results of cooperation rate
    CRdataname = strcat(); 
    data2 = matfile(CRdataname);
    allCR = data2.allCR;
    lcr = [allCR(1,1),allCR(31,2),allCR(496,3),1,0]; % cooperation rate of each typy of players when they fully occupy the population
    for a = 1:25
        alpha = 0.1*(a-1);
        tFXP = zeros(5,5); % init temp result matrix
        for nGt = 1:nGametype
            tFXP(Gametype(nGt,1),Gametype(nGt,2)) = FXP(a,2*nGt);
            tFXP(Gametype(nGt,2),Gametype(nGt,1)) = FXP(a,2*nGt-1);
        end
        T = tFXP/2;
        % get the transition matrix
        for gt = 1:5
            T(gt,gt) = 1-sum(T(gt,:));
        end
        [x,y] = eig(T');
        errdet = abs(diag(y)-1); 
        if min(errdet)<erthd % insurence for error casued by insufficient calculating accuracy
            v = real(x(:,errdet==min(errdet)));
            SMEq=v'/sum(v);
        else 
            SMEq = [0,0,0,0,0];
            disp([strNum,alpha]);
        end
        allSMEq(a,:) = SMEq;
        coop = SMEq*lcr'; % average cooperation rate for the whole population
        allcoop(a) = coop;
    end
    smeqsavepath = strcat();
    if ~exist(smeqsavepath)
        mkdir(smeqsavepath);
    end
    smeqfilename = strcat(smeqsavepath,'\SMEq&coop.mat');
    parsave_smeq(smeqfilename,allSMEq,allcoop);
end
%delete(par)

