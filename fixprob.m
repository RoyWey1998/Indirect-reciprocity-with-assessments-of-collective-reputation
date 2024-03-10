%% from payoff to fixation probability
% this function concerns the population situations with only two strategies
% using the corresponding payoff to calculate the fixation probabilities
function[fxp] = fixprob(s,stg1,stg2,Popstc,SMP)
fxp = zeros(1,2); % init
% find situations with only strategy 1 and 2 
N = max(Popstc(1,:));
ind = Popstc(:,stg1)+Popstc(:,stg2)==N & Popstc(:,stg1)~=0 & Popstc(:,stg2)~=0;
nSMP = SMP(ind,[stg1,stg2]);
payoff1 = nSMP(:,1);
payoff2 = nSMP(:,2); % select corresponding payoff
% calculate fixation probabilities
z = exp(-s*(payoff1(end:-1:1)-payoff2(end:-1:1)));
fxp(1) = 1/(sum(cumprod(z))+1);
z = exp(-s*(payoff2-payoff1));
fxp(2) = 1/(sum(cumprod(z))+1);
end