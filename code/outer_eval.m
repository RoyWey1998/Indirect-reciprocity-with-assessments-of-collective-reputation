% observation and reputation update outside
function[exre] = outer_eval(agtlist,RM,lambda,acts,restagt,str,k,q,err)
koL = length(restagt); % number of Li players outside
exgre = zeros(koL,k); exre = zeros(koL,k);% init
% decide group reputation with resp to each gamer inside
for j = 1:k 
    foc = agtlist(j); 
    eval = agtlist(agtlist~=foc); 
    resum = sum(RM(restagt,eval),2); 
    exgre(resum>=lambda*(k-1),j) = 1;
end
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