%% recovery probability calculating (blocked)
% this file contains code blocks for calculating the recovery probabilities
% of different strategies(mainly L4,6,7,8) under all kinds of settings

%% 1 L4/7 with lambda>(k-2)/(k-1) and lambda<=1/(k-1)
% this condition is equivalent to k=2
% for L4/7, this result is a lower bound
N = 60;
k = 2;
rho1 = 1-2/factorial(N-1);

%% 2 L6/8 with lambda>(k-2)/(k-1) and lambda<=1/(k-1)
% this condition is equivalent to k=2
% for L6/8, this result is the exact value
N = 60;
k = 2;
rho = 1-1/N;

%% 3 L7 with lambda>(k-2)/(k-1) and lambda>1/(k-1)
% for L7, this result is a lower bound
N = 60;
k = 10;
x = (exp(1)*(k-1)/(N-k))^(N-k);
rho1 = 1-x/(N+k-2);

%% 4-1 L8 with lambda>(k-2)/(k-1) and lambda>1/(k-1), lower bound
% for L8, rho1 is a lower bound
N = 60;
k = 10;
sum = 0;
for i = 1:N-k
    prod = 1;
    for j = 1:i
        prod = prod*j*nchoosek(N-j,k-1)/((N-j)*(nchoosek(N-1,k-1)-nchoosek(N-j-1,k-1)));
    end
    sum = sum+prod;
end
sum = sum+prod*(N-k+1)/((k-1)*nchoosek(N-1,k-1));
rho1 = 1-1/(1+sum);

%% 4-2 L8 with lambda>(k-2)/(k-1) and lambda>1/(k-1), upper bound
% for L8, rho2 is a upper bound
sum = 0;
for i = 2:N-k
    prod = 1;
    for j = 2:i
        prod = prod*j*nchoosek(N-j,k-1)/((N-j)*(nchoosek(N-2,k-1)-nchoosek(N-j-1,k-1)+nchoosek(N-j-1,k-2)));
    end
    sum = sum+prod;
end
sum = sum+prod*(N-k+1)/((k-1)*(nchoosek(N-2,k-1)+1));
rho2 = 1-1/(1+(sum+1)/(k-1));

%% L4 with lambda>(k-2)/(k-1) and lambda>1/(k-1), recovery time, upper bound
N = 60;
k = 10;
problist = zeros(2,N-k+2);
for i = 2:N-k+1
    problist(1,i) = (i-1)*(nchoosek(N-i,k-1)+nchoosek(N-2,k-2)-nchoosek(N-i,k-2));
    problist(2,i) = (N-i)*nchoosek(N-i-1,k-2); % p_i^+
end
problist(1,1) = nchoosek(N-1,k-1);
problist(2,1) = (N-1)*nchoosek(N-2,k-2); 
problist(1,N-k+2) = (N-k+1)*(nchoosek(N-2,k-2)-1); % p_i^- for s_{N-k+2}

sum1 = 0; 
for i = 1:N-k+1
    prod = 1;
    for j = 1:i
        prod = prod*(problist(2,j)/problist(1,j+1));
    end
    sum1 = sum1+prod;
end
tau1 = (1+sum1)/problist(1,1)*N*nchoosek(N-1,k-1);

sum2 = 0;
for i = 1:N-k
    prod = 1;
    for j = 1:i
        prod = prod*(problist(2,j+1)/problist(1,j+2));
    end
    sum2 = sum2+prod;
end
tau2 = (1+sum2)/problist(1,2)*N*nchoosek(N-1,k-1);
tau3 = tau1+tau2+N*nchoosek(N-1,k-1)/(nchoosek(N-1,k-1)-nchoosek(N-2,k-1));

tau_up = N*nchoosek(N-1,k-1)/(2*nchoosek(N-1,k-1)-nchoosek(N-2,k-1))+(nchoosek(N-1,k-1)-nchoosek(N-2,k-1))/(2*nchoosek(N-1,k-1)-nchoosek(N-2,k-1))*tau3;

%% L1 with lambda>(k-2)/(k-1) and lambda>1/(k-1), recovery time, upper bound
N = 60;
k = 10;
problist = zeros(2,N-k+2);
for i = 2:N-k+1
    problist(1,i) = (i-1)*nchoosek(N-2,k-1);
    problist(2,i) = (N-i)*nchoosek(N-i-1,k-2); % p_i^+
end
problist(1,1) = nchoosek(N-1,k-1);
problist(2,1) = (N-1)*nchoosek(N-2,k-2);
problist(1,N-k+2) = (N-k+1)*nchoosek(N-2,k-1); % p_i^- for s_{N-k+2}

sum1 = 0; 
for i = 1:N-k+1
    prod = 1;
    for j = 1:i
        prod = prod*(problist(2,j)/problist(1,j+1));
    end
    sum1 = sum1+prod;
end
tau1 = (1+sum1)/problist(1,1)*N*nchoosek(N-1,k-1);

sum2 = 0;
for i = 1:N-k
    prod = 1;
    for j = 1:i
        prod = prod*(problist(2,j+1)/problist(1,j+2));
    end
    sum2 = sum2+prod;
end
tau2 = (1+sum2)/problist(1,2)*N*nchoosek(N-1,k-1);
tau3 = tau1+tau2+N*nchoosek(N-1,k-1)/(nchoosek(N-1,k-1)-nchoosek(N-2,k-1));

tau_up = N*nchoosek(N-1,k-1)/(2*nchoosek(N-1,k-1)-nchoosek(N-2,k-1))+(nchoosek(N-1,k-1)-nchoosek(N-2,k-1))/(2*nchoosek(N-1,k-1)-nchoosek(N-2,k-1))*tau3;

%% L2 with lambda>(k-2)/(k-1) and lambda>1/(k-1), recovery time, lower bound
% this is also L5 with lambda>(k-2)/(k-1) and lambda>1/(k-1), recovery time, lower bound
N = 60;
k = 10;
problist = zeros(2,N);
for i = 1:N-k
    problist(1,i) = i*nchoosek(N-1,k-1); % p_i^-
    problist(2,i) = (N-i)*(nchoosek(N-2,k-1)-nchoosek(N-i-1,k-1)+nchoosek(N-i-1,k-2)); % p_i^+
end
problist(1,N-k+1) = (N-k+1)*nchoosek(N-1,k-1);
problist(2,N-k+1) = (k-1)*(nchoosek(N-2,k-1)+1);
for i = N-k+2:N-1
    problist(1,i) = i*nchoosek(N-1,k-1); % p_i^-
    problist(2,i) = (N-i)*nchoosek(N-2,k-1); % p_i^+
end
problist(1,N) = (N-1)*nchoosek(N-2,k-1); % p_i^-

sum = 0;
for i = 1:N-1
    prod = 1;
    for j = 1:i
        prod = prod*(problist(2,j)/problist(1,j+1));
    end
    sum = sum+prod;
end
tau_low = (1+sum)/problist(1,1)*N*nchoosek(N-1,k-1);

%% L3 with lambda>(k-2)/(k-1) and lambda>1/(k-1), recovery time, upper bound
N = 60;
k = 10;
problist = zeros(2,N-k+2);
for i = 2:N-k+1
    problist(1,i) = (i-1)*(nchoosek(N-1,k-1)-nchoosek(N-i,k-2)); % p_i^-
    problist(2,i) = (N-i)*nchoosek(N-i-1,k-2); % p_i^+
end
problist(1,1) = nchoosek(N-1,k-1);
problist(2,1) = (N-1)*nchoosek(N-2,k-2);
problist(1,N-k+2) = (N-k+1)*(nchoosek(N-1,k-1)-1); % p_i^- for s_{N-k+2}

sum1 = 0; 
for i = 1:N-k+1
    prod = 1;
    for j = 1:i
        prod = prod*(problist(2,j)/problist(1,j+1));
    end
    sum1 = sum1+prod;
end
tau1 = (1+sum1)/problist(1,1)*N*nchoosek(N-1,k-1);

sum2 = 0;
for i = 1:N-k
    prod = 1;
    for j = 1:i
        prod = prod*(problist(2,j+1)/problist(1,j+2));
    end
    sum2 = sum2+prod;
end
tau2 = (1+sum2)/problist(1,2)*N*nchoosek(N-1,k-1);
tau3 = tau1+tau2+N*nchoosek(N-1,k-1)/(nchoosek(N-1,k-1)-nchoosek(N-2,k-1));

tau_up = N*nchoosek(N-1,k-1)/(2*nchoosek(N-1,k-1)-nchoosek(N-2,k-1))+(nchoosek(N-1,k-1)-nchoosek(N-2,k-1))/(2*nchoosek(N-1,k-1)-nchoosek(N-2,k-1))*tau3;

%% L6 with lambda>(k-2)/(k-1) and lambda>1/(k-1), recovery time, lower bound
N = 60;
k = 10;
problist = zeros(2,N);
for i = 2:N-k
    problist(1,i) = (i-1)*(nchoosek(N-i,k-1)+nchoosek(N-2,k-2)-nchoosek(N-i,k-2))+nchoosek(N-i,k-1); % p_i^-
    problist(2,i) = (N-i)*(nchoosek(N-2,k-1)-nchoosek(N-i-1,k-1)+nchoosek(N-i-1,k-2)); % p_i^+
end
problist(1,1) = nchoosek(N-1,k-1);
problist(2,1) = (N-1)*nchoosek(N-2,k-2);
problist(1,N-k+1) = (N-k)*(1+nchoosek(N-2,k-2)-nchoosek(k-1,k-2))+1;
problist(2,N-k+1) = (k-1)*(nchoosek(N-2,k-1)+1);
for i = N-k+2:N
    problist(1,i) = (i-1)*nchoosek(N-2,k-2);
    problist(2,i) = (N-i)*nchoosek(N-2,k-1);
end 
problist(1,N-k+2) = problist(1,N-k+2)-(N-k+1);

sum = 0;
for i = 1:N-1
    prod = 1;
    for j = 1:i
        prod = prod*(problist(2,j)/problist(1,j+1));
    end
    sum = sum+prod;
end
tau_low = (1+sum)/problist(1,1)*N*nchoosek(N-1,k-1);

%% L7 with lambda>(k-2)/(k-1) and lambda>1/(k-1), recovery time, upper bound (L7')
% given that the trace will go back to (1,N-1)
N = 60;
k = 10;
problist = zeros(2,N-k+2);
for i = 2:N-k+1
    %problist(1,i) = (i-1)*nchoosek(N-i,k-1); % p_i^-
    problist(1,i) = (i)*nchoosek(N-i,k-1); % p_i^-
    %problist(2,i) = (N-i)*nchoosek(N-i-1,k-2); % p_i^+
    problist(2,i) = max([(N-i)*nchoosek(N-i-1,k-2),nchoosek(N-1,k-1)-nchoosek]); % p_i^+
end
problist(1,1) = nchoosek(N-1,k-1);
problist(2,1) = (N-1)*nchoosek(N-2,k-2);
problist(1,N-k+2) = 2*nchoosek(N-1,k-1)-nchoosek(N-2,k-1); % p_i^- for s_{N-k+2}

sum1 = 0; 
for i = 1:N-k+1
    prod = 1;
    for j = 1:i
        prod = prod*(problist(2,j)/problist(1,j+1));
    end
    sum1 = sum1+prod;
end
tau1 = (1+sum1)/problist(1,1)*N*nchoosek(N-1,k-1);

sum2 = 0;
for i = 1:N-k
    prod = 1;
    for j = 1:i
        prod = prod*(problist(2,j+1)/problist(1,j+2));
    end
    sum2 = sum2+prod;
end
tau2 = (1+sum2)/problist(1,2)*N*nchoosek(N-1,k-1);
tau3 = tau1+tau2+N*nchoosek(N-1,k-1)/(nchoosek(N-1,k-1)-nchoosek(N-2,k-1));

tau_up = N*nchoosek(N-1,k-1)/(2*nchoosek(N-1,k-1)-nchoosek(N-2,k-1))+(nchoosek(N-1,k-1)-nchoosek(N-2,k-1))/(2*nchoosek(N-1,k-1)-nchoosek(N-2,k-1))*tau3;

%% L7 with lambda>(k-2)/(k-1) and lambda>1/(k-1), recovery time, upper bound (L7'')
N = 60;
k = 10;
problist = zeros(2,N-k+2);
for i = 2:N-k
    problist(1,i) = (i)*nchoosek(N-i,k-1); % p_i^-
    problist(2,i) = max([(N-i)*nchoosek(N-i-1,k-2),nchoosek(N-1,k-1)-nchoosek(N-i-1,k-1)]); % p_i^+
end
problist(1,1) = nchoosek(N-1,k-1);
problist(2,1) = nchoosek(N-1,k-1)-nchoosek(N-2,k-1);
problist(1,N-k+1) = N-k+1;
problist(2,N-k+1) = nchoosek(N-1,k-1);
problist(1,N-k+2) = 2*nchoosek(N-1,k-1)-nchoosek(N-2,k-1); % p_i^- for s_{N-k+2}

sum = 0;
for i = 1:N-k+1
    prod = 1;
    for j = 1:i
        prod = prod*(problist(2,j)/problist(1,j+1));
    end
    sum = sum+prod;
end
tau_up = (1+sum)/problist(1,1)*N*nchoosek(N-1,k-1);