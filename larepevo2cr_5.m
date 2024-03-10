%% Reputation evolution process for assessment criterion (3-lambda with ALC/ALLD version)
% 2024.1.28
function[CR,crspt1,crspt2,crspt3] = larepevo2cr_5(str,popstc,N,ngs,repItnum,kmin,kmax,Gamestc,q,err,lambda_list)
% this version of reputation iteration contains 3 types of Li players along
% with ALLC and ALLD players

% distribution of strategies
Lnum1 = popstc(1); Lnum2 = popstc(2); Lnum3 = popstc(3); % number of players with 3 values of lambda
Lnum = Lnum1+Lnum2+Lnum3; Cnum = popstc(4); Dnum = popstc(5); % number of players with 3 strategies
Lpop1 = 1:Lnum1; Lpop2 = Lnum1+1:Lnum1+Lnum2; Lpop3 = Lnum1+Lnum2+1:Lnum; % population of 3 lambda each
Lpop = 1:Lnum; Cpop = Lnum+1:Lnum+Cnum; Dpop = Lnum+Cnum+1:N;

% reputation init
%M = [randi([0,1],Lnum,N);ones(Cnum,N);zeros(Dnum,N)]; % initial reputation matrix(random)
M = [ones(Lnum,N);ones(Cnum,N);zeros(Dnum,N)]; % initial reputation matrix(good)
%M = [zeros(Lnum,N);ones(Cnum,N);zeros(Dnum,N)]; % initial reputation matrix(bad)

% game time and coopration statistic init
Gs = 1:ngs; %index list of all game state
GameTime = zeros(1,N); % init of game time statistic
CoopTime = zeros(1,N); % init of cooperate time statistic
GameTimespt = zeros(ngs,N); % init of separate game time 
CoopTimespt = zeros(ngs,N); % init of separate cooperate time

% game process
k = 10; % for fixed group size
for it = 1:repItnum
    %k = randi([kmin,kmax]); % for random group size
    temp = randperm(N);
    agtlist = sort(temp(1:k)); % select gamers
    % find players of 3 lambda in-group
    Lagt = intersect(agtlist,Lpop);
    Lagt1 = intersect(agtlist,Lpop1); 
    Lagt2 = intersect(agtlist,Lpop2); 
    Lagt3 = intersect(agtlist,Lpop3);
    Cagt = intersect(agtlist,Cpop); 
    Dagt = intersect(agtlist,Dpop);
    nAgtin = [length(Lagt1),length(Lagt2),length(Lagt3),length(Cagt),length(Dagt)]; % number of each type of players in game
    % decide game state
    gamestc = Gs(Gamestc(:,1) == length(Lagt1) & Gamestc(:,2) == length(Lagt2) & Gamestc(:,3) == length(Lagt3) & Gamestc(:,4) == length(Cagt) & Gamestc(:,5) == length(Dagt));
    % agents outside
    nAgtout = popstc-nAgtin;
    restagt = setdiff(Lpop,agtlist); % Li 
    % acting stage & inner reputation update
    [acts,inre] = lainner_game_5(agtlist,M,lambda_list,Lagt,nAgtin,str,k);
    M(Lagt,agtlist) = inre; % update reputations within group
    % observation and reputation update outside
    exre = laouter_eval_5(agtlist,M,lambda_list,acts,restagt,nAgtout,str,k,q,err);
    M(restagt,agtlist) = exre; % update reputation outside group

    % cooperation rate record
    if it > repItnum/4 % after one forth of iterations
        GameTime(agtlist) = GameTime(agtlist)+1; % update game times
        CoopTime(agtlist) = CoopTime(agtlist)+acts; % update cooperate times
        GameTimespt(gamestc,agtlist) = GameTimespt(gamestc,agtlist)+1; % update separating game times
        CoopTimespt(gamestc,agtlist) = CoopTimespt(gamestc,agtlist)+acts; % update separating cooperate times
    end
end

% cooperation rate statistic
popcooprat = CoopTime(Lpop)./GameTime(Lpop); % cooperation rate of each Li player
CR = [mean(popcooprat(Lpop1)),mean(popcooprat(Lpop2)),mean(popcooprat(Lpop3))]; % average cooperation rate of 3 lambda-populations
popcoopratspt = CoopTimespt(:,Lpop)./GameTimespt(:,Lpop); % separating cooperation rate of each Li player
popcoopratspt(isnan(popcoopratspt)) = 0; % replace invalid values with 0
crspt1 = mean(popcoopratspt(:,Lpop1),2)';
crspt2 = mean(popcoopratspt(:,Lpop2),2)'; 
crspt3 = mean(popcoopratspt(:,Lpop3),2)'; % average separating cooperation rate of 3 lambda-populations
end