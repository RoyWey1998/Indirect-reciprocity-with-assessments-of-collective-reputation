%% Reputation evolution process for assessment criterion (3-lambda version)
% 2024.1.27
function[CR,crspt1,crspt2,crspt3] = larepevo2cr(str,popstc,N,ngs,repItnum,kmin,kmax,Gamestc,q,err,lambda_list)
% the 3-lambda version contains only one type of leading-eight strategy 
% with3 values of lambda 

% distribution of 3 lambda values
num1 = popstc(1); num2 = popstc(2); %num3 = popstc(3); % number of players with 3 values of lambda
pop = 1:N; % the whole population
pop1 = 1:num1; pop2 = num1+1:num1+num2; pop3 = num1+num2+1:N; % population of 3 lambda each

% reputation init
%M = randi([0,1],N,N); % initial reputation matrix(random)
M = ones(N,N); % initial reputation matrix(good)
%M = zeros(N,N); % initial reputation matrix(bad)

% game time and coopration statistic init
Gs = 1:ngs; %index list of all game state
GameTime = zeros(1,N); % init of game time statistic
CoopTime = zeros(1,N); % init of cooperate time statistic
GameTimespt = zeros(ngs,N); % init of separate game time 
CoopTimespt = zeros(ngs,N); % init of separate cooperate time

% game process
k = 10; % for fixed group size
for it  = 1:repItnum
    %k = randi([kmin,kmax]); % for random group size
    temp = randperm(N);
    agtlist = sort(temp(1:k)); % select gamers
    % find players of 3 lambda in-group 
    agt1 = intersect(agtlist,pop1); 
    agt2 = intersect(agtlist,pop2); 
    agt3 = intersect(agtlist,pop3);
    nAgtin = [length(agt1),length(agt2),length(agt3)];
    % decide game state
    gamestc = Gs(Gamestc(:,1) == length(agt1) & Gamestc(:,2) == length(agt2) & Gamestc(:,3) == length(agt3));
    % agents out of group
    nAgtout = popstc-nAgtin;
    restagt = setdiff(pop,agtlist);
    % acting stage and inner reputation update
    [acts,inre] = lainner_game(agtlist,M,lambda_list,nAgtin,str,k);
    M(agtlist,agtlist) = inre; % update reputations within group
    % observation and reputation update outside
    exre = laouter_eval(agtlist,M,lambda_list,acts,restagt,nAgtout,str,k,q,err);
    M(restagt,agtlist) = exre; % update reputation outside group

    % cooperation rate record
    if it >= repItnum/4 % after one forth of iterations
        GameTime(agtlist) = GameTime(agtlist)+1; % update game times
        CoopTime(agtlist) = CoopTime(agtlist)+acts; % update cooperate times
        GameTimespt(gamestc,agtlist) = GameTimespt(gamestc,agtlist)+1; % update separating game times
        CoopTimespt(gamestc,agtlist) = CoopTimespt(gamestc,agtlist)+acts; %update separating cooperate times
    end
end

% cooperation rate statistic
popcooprat = CoopTime(pop)./GameTime(pop); % cooperation rate of each player
CR = [mean(popcooprat(pop1)),mean(popcooprat(pop2)),mean(popcooprat(pop3))]; % average cooperation rate of 3 lambda-populations
popcoopratspt = CoopTimespt./GameTimespt; % separating cooperation rate of each player
popcoopratspt(isnan(popcoopratspt)) = 0; % replace invalid values with 0
crspt1 = mean(popcoopratspt(:,pop1),2)';
crspt2 = mean(popcoopratspt(:,pop2),2)'; 
crspt3 = mean(popcoopratspt(:,pop3),2)'; % average separating cooperation rate of 3 lambda-populations
end