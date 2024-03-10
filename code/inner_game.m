% inner game process
% inside evalution with sympathetic perspective
% this function can degenerate to the classic model
function[acts,inre] = inner_game(agtlist,RM,lambda,Lagt,kL,kC,kD,str,k)
acts = [ones(1,kL+kC),zeros(1,kD)];
ingre = zeros(kL,k); 
inre = zeros(kL,k); %init
% act stage & evaluate the corresponding recipient groups with respect to Li gamers
for j=1:kL
    foc = Lagt(j); % focal gamer
    eval = agtlist(agtlist ~= foc); % exclude focal to get evaluated gamers
    resum = sum(RM(Lagt,eval),2); % respectively sum up the images of evaluated gamers from the view of all gamers
    ingre(resum>=lambda*(k-1),j) = 1; % decide the group reputations
    idx = RM(foc,foc)*2 + ingre(j,j);
    acts(j) = str(idx+9); % act
end
% evaluate the corresponding recipient groups with respect to other gamers
for j = kL+1:k
    foc = agtlist(j); % focal gamer
    eval = agtlist(agtlist ~= foc); % exclude focal to get evaluated gamers
    resum = sum(RM(Lagt,eval),2); % respectively sum up the images of evaluated gamers from the view of all gamers
    ingre(resum>=lambda*(k-1),j) = 1; % decide the group reputations
end
% reputation update inside group
for j = 1:kL 
    foc = Lagt(j);
    idx1 = RM(foc,agtlist)*4+acts*2+ingre(j,:); % find the corresponding index
    inre(j,:) = str(idx1+1);
end
end