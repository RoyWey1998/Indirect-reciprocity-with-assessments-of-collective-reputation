% tool code (blocked)
%% 1 
% Find possible game states under different population structures 
% and decide the probabilities for players of each strategy taking
% part in each situation (3-strategy version)
Popstc = Popstc60; 
Gamestc = Gamestc2_10;
nPopstc = size(Popstc,1); % number of population structures
nGamestc = size(Gamestc,1); % number of game states
Lstcprob = zeros(nPopstc,nGamestc);
Cstcprob = zeros(nPopstc,nGamestc);
Dstcprob = zeros(nPopstc,nGamestc); % init
for nPs = 1:nPopstc
    for nGs = 1:nGamestc
        if Popstc(nPs,1)>=Gamestc(nGs,1) && Popstc(nPs,2)>=Gamestc(nGs,2) && Popstc(nPs,3)>=Gamestc(nGs,3)
            n1 = Popstc(nPs,1); n2 = Popstc(nPs,2); n3 = Popstc(nPs,3);
            k1 = Gamestc(nGs,1); k2 = Gamestc(nGs,2); k3 = Gamestc(nGs,3);
            k = sum(Gamestc(nGs,:));
            if k1>=1
                Lstcprob(nPs,nGs) = nchoosek(n1-1,k1-1)*nchoosek(n2,k2)*nchoosek(n3,k3)/nchoosek(N-1,k-1); 
            end
            if k2>=1
                Cstcprob(nPs,nGs) = nchoosek(n1,k1)*nchoosek(n2-1,k2-1)*nchoosek(n3,k3)/nchoosek(N-1,k-1);
            end
            if k3>=1
                Dstcprob(nPs,nGs) = nchoosek(n1,k1)*nchoosek(n2,k2)*nchoosek(n3-1,k3-1)/nchoosek(N-1,k-1);
            end
        end
    end
end
filename = 'GamestcProb.mat';
save(filename,'Lstcprob','Cstcprob','Dstcprob');

%% 2 payoff to selection-mutation equilibrium (\mu \to 0)
% from payoff to fixation probability and selection-mutation equilibrium
% under the limitation of rare mutations (3-strategy version)
load('Popstc60_3.mat');
%   load('Gamestc10_3.mat');
%   load('GamestcProb10_3.mat');
%this two lines below are for unfixed k
load("Gamestc2_10_3.mat");
load('GamestcProb2_10_3.mat');

% payoff to fixation probability (for fixed alpha)
par = parpool('local',13);
Popstc = Popstc60;
%Gamestc = Gamestc10;
Gamestc = Gamestc2_10; % unfixed_k
pk = [kmin:kmax]/sum(kmin:kmax);
Gametype = [1,2;1,3;2,3]; % 3-strategy
nGametype = size(Gametype,1);
parfor a = 1:26
    disp(a);
    %alpha = 1*(a-1);%beta=0
    %alpha = 0.1*(a-1);%beta=1
    alpha = 0.01*(a-1);%beta=2
    for strNum = 1:8
        FXP = zeros(9,2*nGametype); % 9 rows for 9 lambdas, 6 columns for pairwise fixation probabilities among 3 strategies
        for l = 1:9
            lambda = 0.1*l;
            % load LCR and LCRspt
            Dataname = strcat();
            %load(Dataname);  
            data = matfile(Dataname);
            LCRspt = data.LCRspt;
            % calculate average payoff with cooperation rate
            SMP = cr2pof(Popstc,Gamestc,LCRspt,alpha,beta,c,pk,kmin,Lstcprob,Cstcprob,Dstcprob);
            SMPsavepath = strcat();
            if ~exist(SMPsavepath)
                mkdir(SMPsavepath);
            end            
            SMPname = strcat(SMPsavepath,'\L',num2str(strNum),'_λ=',num2str(lambda),'_SMP.mat');
            parsave_SMP(SMPname,SMP);
            % the four lines below skip the SMP calculations 
%             SMPsavepath = strcat();
%             SMPname = strcat(SMPsavepath,'\L',num2str(strNum),'_λ=',num2str(lambda),'_SMP.mat');
%             SMPdata = matfile(SMPname);
%             SMP = SMPdata.SMP;
            for gt = 1:nGametype
                stg1 = Gametype(gt,1);
                stg2 = Gametype(gt,2);
                fxp = fixprob(s,stg1,stg2,Popstc,SMP);
                FXP(l,2*gt-1:2*gt) = fxp;
            end
        end
        FXPsavepath = strcat();
        if ~exist(FXPsavepath)
            mkdir(FXPsavepath);
        end
        filename = strcat(FXPsavepath,'\L',num2str(strNum),'fixprob.mat');
        parsave_FXP(filename,FXP);
    end
    
    % from fixation probabilities to selection-mutation equilibrium and average
    % cooperation rate
    erthd = 1e-14; % calculation error threshold, with large alpha, the 
    % eigenvalue corresponding to the null space of T'-eye(3)  may not be 
    % accurately 0 due to the limitation of calculation accuracy
    for strNum = 1:8
        allSMEq = zeros(9,3);
        allcoop = zeros(9,1); % init
        % load fixation probabilities
        filename = strcat(FXPsavepath,'\L',num2str(strNum),'fixprob.mat');
        %load(filename); 
        data1 = matfile(filename);
        FXP = data1.FXP;
        for l = 1:9
            tFXP = zeros(3,3); % init temp result matrix
            lambda = 0.1*l;
            % load corresponding results of cooperation rate
            CRdataname = strcat(); 
            %load(CRdataname); 
            data2 = matfile(CRdataname);
            LCR = data2.LCR;
            lcr = LCR(1); % the cooperation rate applied here should be the one within all-Li population
            tFXP(2,1) = FXP(l,1); % C trans2 Li
            tFXP(1,2) = FXP(l,2); % Li trans2 C
            tFXP(3,1) = FXP(l,3); % D trans2 Li
            tFXP(1,3) = FXP(l,4); % Li trans2 D
            tFXP(3,2) = FXP(l,5); % D trans2 C
            tFXP(2,3) = FXP(l,6); % C trans2 D
            T = tFXP/2;
            % get the transition matrix
            for gt = 1:3
                T(gt,gt) = 1-sum(T(gt,:));
            end
            [x,y] = eig(T');
            errdet = abs(diag(y)-1);
            if min(errdet)<erthd
                v = real(x(:,errdet==min(errdet)));
                SMEq=abs(v'/sum(v));
            else 
                SMEq = [0,0,0];
                disp([strNum,lambda]);
            end
            %v=null(T'-eye(3)); % left eigenvector corresponding to 1
            %SMEq=v'/sum(v); % normalization
            allSMEq(l,:) = SMEq;
            coop = SMEq*[lcr;1;0]; % average cooperaion rate under the set of parameter
            allcoop(l) = coop;
        end
        smeqsavepath = strcat();
        if ~exist(smeqsavepath)
            mkdir(smeqsavepath);
        end
        smeqfilename = strcat(smeqsavepath,'\L',num2str(strNum),'_SMEq&coop.mat');
        parsave_smeq(smeqfilename,allSMEq,allcoop);
    end
end
delete(par)
%% 3 even distribution reputation result
% for the special case where 3 strategis are even distributed in the
% population
tic
load('Gamestc10_3.mat');
load('Popstc60_3.mat');
default;
err = 0.05; % error rate
coreNum = 8;
par = parpool('local',coreNum);
parfor strNum = 1:8
    disp(strNum);
    str = strategy_list(strNum,:); % type of strategy
    for l = 1:9 % parfor endpoint must be internal
        lambda = 0.1*l; % assessment criteria
        nPs = 1051;
        Lnum = Popstc60(nPs,1); Cnum = Popstc60(nPs,2); Dnum = Popstc60(nPs,3); % numbers of players with 3 strategies
        [Lcr,Lcrspt,M] = repevo2cr(str,Lnum,Cnum,Dnum,N,repItnum,kmin,kmax,Gamestc10,lambda,q,err);        
        filename = strcat();
        parsaveM(filename,Lcr,M);
    end
end
delete(par)
toc
%% 4 payoff to average frequency (\mu > 0)
% from payoff to average frequency and cooperation rate at 
% selection-mutation equilibrium with significant mutation rate (3-strategy version)
Popstc = Popstc60;
Gamestc = Gamestc10;
pk = [kmin:kmax]/sum(kmin:kmax);
mu = 0.05;
par = parpool('local',8);
parfor a = 1:16
    disp(a);
    alpha = 0.1*(a-1);
    for strNum = 1:8    
        allAvFr = zeros(9,3); % init
        allcoop = zeros(9,1);
        for l = 1:9
            lambda = 0.1*l;
            Dataname1 = strcat(); % load payoff data SMP
            data1 = matfile(Dataname1);
            SMP = data1.SMP;
            Dataname2 = strcat(); % load cooperation rate LCR
            data2 = matfile(Dataname2);
            LCR = data2.LCR;
            LCR(isnan(LCR)) = 0; % replace invalid value with 0 
            [PopCoop,AvFr] = slc_mut_equ(N,nPopstc,mu,s,Popstc,SMP,LCR); % calculate s-m equilibrium
            allcoop(l) = PopCoop;
            allAvFr(l,:) = AvFr';
        end
        avfrsavepath = strcat();
        if ~exist(avfrsavepath)
            mkdir(avfrsavepath);
        end
        avfrfilename = strcat(avfrsavepath,'\L',num2str(strNum),'_AvFr&coop.mat');
        parsave_avfr(avfrfilename,allAvFr,allcoop);
    end
end
delete(par)
%% 5 calculate alpha_c for unfixed k with different beta
beta=2;
sum = 0;
for k=2:10
    sum = sum+pk(k-1)*k^(beta-1);
end
alpha_c = 1/sum;

%% 6 cooperation rate converge
tic
load('Gamestc10_3.mat');
load('Popstc60_3.mat');
default;
err = 0.05; % error rate
coreNum = 4;
par = parpool('local',coreNum);
for l = 2:3:8 
    lambda = 0.1*l; % assessment criteria
    disp(lambda); % monitoring flag
    LCR = [];
    parfor strNum = 1:8 % parfor endpoint must be internal
        str = strategy_list(strNum,:); % type of strategy       
        nPs = 1051;
        Lnum = Popstc60(nPs,1); Cnum = Popstc60(nPs,2); Dnum = Popstc60(nPs,3); % numbers of players with 3 strategies
        [Lcr, Lcrspt, tcr] = repevo2cr(str,Lnum,Cnum,Dnum,N,repItnum,kmin,kmax,Gamestc10,lambda,q,err);
        LCR = [LCR;tcr];
    end
    filename = strcat();
    save(filename,'LCR');
end
delete(par)
toc

%% 7 run for payoff to verify cr2pof
tic
load('Gamestc10_3.mat');
load('Popstc60_3.mat');
default;
err = 0.05; % error rate
coreNum = 5;
par = parpool('local',coreNum);
Pof = [];
alpha = 0.5;
strNum = 6; 
str = strategy_list(strNum,:); % type of strategy       
nPs = 1276;
parfor l = 1:9 
    lambda = 0.1*l; % assessment criteria
    disp(lambda); % monitoring flag    
    Lnum = Popstc60(nPs,1); Cnum = Popstc60(nPs,2); Dnum = Popstc60(nPs,3); % numbers of players with 3 strategies
    [strmeanpof] = repevo2cr_pofvrf(str,Lnum,Cnum,Dnum,N,repItnum,kmin,kmax,Gamestc10,lambda,q,err,alpha,beta);
    Pof = [Pof;strmeanpof];
end
filename = strcat();
save(filename,'Pof');
delete(par)
toc
