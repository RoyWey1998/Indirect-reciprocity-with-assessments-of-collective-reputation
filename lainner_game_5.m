% inner game (with traditional perspective) for assessment criterion game process 
% 3-lambda with ALLC/ALLD version
function[acts,inre] = lainner_game_5(agtlist,RM,lambda_list,Lagt,nAgtin,str,k)
kL = length(Lagt); % number of Li gamers
% lambda list for Li gamers
Lagtla = [lambda_list(1)*ones(1,nAgtin(1)),lambda_list(2)*ones(1,nAgtin(2)),lambda_list(3)*ones(1,nAgtin(3))];
% init
acts = [ones(1,k-nAgtin(5)),zeros(1,nAgtin(5))];
ingre = zeros(kL,k);
inre = zeros(kL,k);
% act stage & evaluate the corresponding recipient groups with respect to Li gamers
for j = 1:kL
    foc = Lagt(j); % focal gamer
    eval = agtlist(agtlist ~= foc); % exclude focal to get evaluated gamers
    resum = sum(RM(Lagt,eval),2); % respectively sum up the images of evaluated gamers from the view of all gamers
    ingre(resum>=Lagtla'*(k-1),j) = 1; % decide the group reputations
    idx = RM(foc,foc)*2 + ingre(j,j); % acting index
    acts(j) = str(idx+9); % act
end
% evaluate the corresponding recipient groups with respect to other gamers
for j = kL+1:k
    foc = agtlist(j); % focal gamer
    eval = agtlist(agtlist ~= foc); % exclude focal to get evaluated gamers
    resum = sum(RM(Lagt,eval),2); % respectively sum up the images of evaluated gamers from the view of all gamers
    ingre(resum>=Lagtla'*(k-1),j) = 1; % decide the group reputations
end
% reputation update inside group
for j = 1:kL
    foc = Lagt(j);
    idx1 = RM(foc,agtlist)*4+acts*2+ingre(j,:); % reputation index
    inre(j,:) = str(idx1+1);
end
end