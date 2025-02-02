%% parameter setup for lambda evolution (5-strategy version)
% game parameter
repItnum = 2e5; % iteration times
N = 30; % population size
nPopstc = size(Popstc30_5,1); % number of population structures
ngs = size(Gamestc10_5,1); % number of game structures
%err = 0.05; % error rate
q = 0.9; % observation probability
kmax = 10;
kmin = 10;
% payoff parameter
%alpha = 1;
beta = 1;
c = 1;
% evolution parameter
s = 1;
%mu = 1;
% strategy setup
L_1 = [0,0,1,1,1,0,1,1,1,1,0,1];
L_2 = [0,0,1,1,1,0,0,1,1,1,0,1];
L_3 = [1,0,1,1,1,0,1,1,0,1,0,1];
L_4 = [1,0,0,1,1,0,1,1,0,1,0,1];
L_5 = [1,0,1,1,1,0,0,1,0,1,0,1];
L_6 = [1,0,0,1,1,0,0,1,0,1,0,1];
L_7 = [0,0,0,1,1,0,1,1,0,1,0,1];
L_8 = [0,0,0,1,1,0,0,1,0,1,0,1];
ALLC = [1,1,1,1,1,1,1,1,1,1,1,1];
ALLD = [0,0,0,0,0,0,0,0,0,0,0,0];
strategy_list = [L_1; L_2; L_3; L_4; L_5; L_6; L_7; L_8; ALLC; ALLD];