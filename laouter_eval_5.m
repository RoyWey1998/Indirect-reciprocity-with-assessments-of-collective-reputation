% observation and reputation update outside for assessment criterion game process
% (3-lambda with ALLC/ALLD type)
function[exre] = laouter_eval_5(agtlist,RM,lambda_list,acts,restagt,nAgtout,str,k,q,err)
koL = length(restagt); % number of Li players outside
exgre1 = zeros(nAgtout(1),k); exgre2 = zeros(nAgtout(2),k);exgre3 = zeros(nAgtout(3),k);exre = zeros(koL,k); % init
% decide group reputation with resp to each gamer inside
for j = 1:k
    foc = agtlist(j); % focal gamer
    eval = agtlist(agtlist~=foc); % exclude focal to get evaluated gamers
    resum = sum(RM(restagt,eval),2); % respectively sum up the images of evaluated gamers from the view of all agents
    exgre1(resum(1:nAgtout(1))>=lambda_list(1)*(k-1),j) = 1;
    exgre2(resum(nAgtout(1)+1:nAgtout(1)+nAgtout(2))>=lambda_list(2)*(k-1),j) = 1;
    exgre3(resum(nAgtout(1)+nAgtout(2)+1:koL)>=lambda_list(3)*(k-1),j) = 1; % decide the group reputation of recipients
end
exgre = [exgre1;exgre2;exgre3];
% update reputation image for each Li player outside
for j = 1:koL
    ob = restagt(j); % observing player
    if rand(1)<q % if observe
        err1 = rand(1,k); 
        obacts = acts; % init observation results
        obacts(err1<err) = ~acts(err1<err); % observe with error
        idx = RM(ob,agtlist)*4+obacts*2+exgre(j,:); % decide situation index
        exre(j,:) = str(idx+1); % new reputations after observation
    else
        exre(j,:) = RM(ob,agtlist); % reputations unchanged without observation
    end
end
end