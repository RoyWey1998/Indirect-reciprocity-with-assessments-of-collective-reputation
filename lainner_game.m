% inner game (with traditional perspective) for assessment criterion game process
function[acts,inre] = lainner_game(agtlist,RM,lambda_list,nAgtin,str,k)
% lambda list for gamers
agtla = [lambda_list(1)*ones(1,nAgtin(1)),lambda_list(2)*ones(1,nAgtin(2)),lambda_list(3)*ones(1,nAgtin(3))];
% init
acts = ones(1,k);
ingre = zeros(k,k); % here we use the traditional inner game model
inre = zeros(k,k);
% act stage & evaluate the corresponding recipient groups with respect to Li gamers
for j = 1:k
    foc = agtlist(j); % focal gamer
    eval = agtlist(agtlist ~= foc); % exclude focal to get evaluated gamers
    resum = sum(RM(agtlist,eval),2); % respectively sum up the images of evaluated gamers from the view of all gamers
    ingre(resum>=agtla'*(k-1),j) = 1; % decide the group reputations
    idx = RM(foc,foc)*2 + ingre(j,j); % acting index
    acts(j) = str(idx+9); % act
end
% reputation update inside group
for j = 1:k 
    foc = agtlist(j); % focal updater
    idx1 = RM(foc,agtlist)*4+acts*2+ingre(j,:); % reputation index
    inre(j,:) = str(idx1+1);
end
end