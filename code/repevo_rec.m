%% Reputation evolution for recovery analysis
function[recst,recprob,rectime] = repevo_rec(str,N,repTime,repItnum,k,lambda)
% this function use the private assessment procedure to conduct the
% reputation evolving process for the special case of recovery from single
% disagreement, inner game is the traditional version

recst = zeros(1,repTime); % init the recovery result statistic
rectime = zeros(1,repTime); % init the recovery time statistic

% game process
parfor rt = 1:repTime
    %if mod(rt,100)==1
        disp(rt);
    %end
    % reputation init
    M = ones(N,N); M(1,2) = 0; % initial reputation matrix(all good with one disagreement)
    for it = 1:repItnum
        rpgre = zeros(N,1); % init the vector for reputations of the recipient group
        temp = randperm(N);
        donor = temp(1); % select one player as the donor
        rplist = sort(temp(2:k)); % select recipients
        % act stage
        act = sum(M(donor,rplist))>=lambda*(k-1);
        % evaluate stage
        rpgre(sum(M(:,rplist),2)>=lambda*(k-1))=1; % evaluate the reputation
        idx = M(:,donor)*4+act*2+rpgre; % decide situation index
        donre = str(idx+1);
        % reputation update
        M(:,donor) = donre;
        if ismember(0,M)==0 % if already recover
            recst(rt) = 1; % record the result for recover
            rectime(rt) = it; % record the time for recover
            break
        end
    end
end
% calculate the rate of recovery
recprob = sum(recst)/repTime;
end