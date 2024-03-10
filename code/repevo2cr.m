%% Reputation Evolution
% 2024.1.25
function[Lcr, Lcrspt] = repevo2cr(str,Lnum,Cnum,Dnum,N,repItnum,kmin,kmax,Gamestc,lambda,q,err)
% Leading-eight strategy, number of Li ALLC ALLD, total number, number of game states,
% iteration number, group size bound, collection of game states, assessment
% criteria, observaton probability, error rate

% distribution of 3 strategy
Lpop = 1:Lnum;
Cpop = Lnum+1:Lnum+Cnum;
Dpop = Lnum+Cnum+1:N;

% reputation init
%M = [randi([0,1],Lnum,N);ones(Cnum,N);zeros(Dnum,N)]; % initial reputation matrix(random)
M = [ones(Lnum,N);ones(Cnum,N);zeros(Dnum,N)]; % initial reputation matrix(good)
%M = [zeros(Lnum,N);ones(Cnum,N);zeros(Dnum,N)]; % initial reputation matrix(bad)

% game time and coopration statistic init
ngs = size(Gamestc,1); % number of all game state
Gs = 1:ngs; %index list of all game state
GameTime = zeros(1,N); % init of game time statistic
CoopTime = zeros(1,N); % init of cooperate time statistic
GameTimespt = zeros(ngs,N); % init of separate game time 
CoopTimespt = zeros(ngs,N); % init of separate cooperate time
%tcr = zeros(1,3*repItnum/4); % test for convergence

% game process
%k = 10; % for fixed group size
for it = 1:repItnum
    k = randi([kmin,kmax]); % for random group size
    temp = randperm(N);
    agtlist = sort(temp(1:k)); % select gamers
    % find gamers of 3 strategies
    Lagt = intersect(agtlist,Lpop); 
    Cagt = intersect(agtlist,Cpop); 
    Dagt = intersect(agtlist,Dpop);
    kL = length(Lagt); % number of Li gamer
    kC = length(Cagt); % number of ALLC gamer
    kD = length(Dagt); % number of ALLD gamer
    % decide game state
    gamestc = Gs(Gamestc(:,1) == kL & Gamestc(:,2) == kC & Gamestc(:,3) == kD);
    % Li agents outside
    restagt = setdiff(Lpop,agtlist);
    % acting stage & inner reputation update
    [acts,inre] = inner_game(agtlist,M,lambda,Lagt,kL,kC,kD,str,k);
    M(Lagt,agtlist) = inre; % update reputations within group
    % observation and reputation update outside
    exre = outer_eval(agtlist,M,lambda,acts,restagt,str,k,q,err);
    M(restagt,agtlist) = exre; % update reputation outside group

    % cooperation rate record
    if it > repItnum/4 % after one forth of iterations
        GameTime(agtlist) = GameTime(agtlist)+1; % update game times
        CoopTime(agtlist) = CoopTime(agtlist)+acts; % update cooperate times
        GameTimespt(gamestc,agtlist) = GameTimespt(gamestc,agtlist)+1; % update separating game times
        CoopTimespt(gamestc,agtlist) = CoopTimespt(gamestc,agtlist)+acts; % update separating cooperate times
        % test for convergence
        %tpopcr = CoopTime(Lpop)./GameTime(Lpop);
        %tcr(it-repItnum/4) = mean(tpopcr);
    end
end

% cooperation rate statistic
popcooprat = CoopTime(Lpop)./GameTime(Lpop); % cooperation rate of each Li player
Lcr = mean(popcooprat); % average cooperation rate of Li population
popcoopratspt = CoopTimespt(:,Lpop)./GameTimespt(:,Lpop); % separating cooperation rate of each Li player
popcoopratspt(isnan(popcoopratspt)) = 0; 
Lcrspt = mean(popcoopratspt,2)'; % separating average cooperation rate of Li population

% test for convergece
%figure(1);
%plot(tcr);
%ylim([-0.05, 1.05]);

% image matrix
%figure(2);
%imagesc(M);
end