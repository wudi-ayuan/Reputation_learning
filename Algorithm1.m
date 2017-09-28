clear all
%In this Monte Carlo simulation, we are given a fixed initial network and
%fixed true quality distribution means and variance. 
%We test this network M times with differnt settings of its users.
%Different settings include users' true quality which leads to differnt
%distribution of unscaled hitting time.
%_____I/O______
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  INPUT LIST: %%%%%%%%%%%%%%%%%%%%%%%
%the number of draws of the Monte Carlo simulation(M), a vector of the agent's
%prior quality distribution means (mu)
%, a vector of the agent's prior quality distribution variances (var)
%, a vector of the agent's base precisions(tau), the linking cost (c), the
%discount rate (p), and an initial network structure(A)
%Monte Carlo M, N, E, Graph A, base precision, hitting time hitT


%during the trial, the parameters should be carefully decided. 
%G(N,sum(k))
N=60;               %number of player
M=1000;             %number of draws, Monte Cralo simulation
c=9;              %cost 
mu=c+4*rand(1,N);%true quality distribution means, should be greater than c (9,13)
var=10+rand(1,N);%true quality distribtuion variance, sigma^2
p=0.1;              %discount rate 
%base precision, uniform distribution [2,7]
tau_min=0.001;
tau=tau_min+0.001*rand(1,N);
%Graph adjacent matrix A
A=randi([0,1],N,N);
B=triu(A,1);
A=B+B';% undirected graph, symmetric adjacent matrix
for i=1:N, 
   k(i)=sum(A(i,:));    %player's rank
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  OUTPUT LIST: %%%%%%%%%%%%%%%%%%%%%%%
W=zeros(M,N);
avgW=[];

%_____INIT______
for i=1:M,                       %initiate quality matrix(M*N)
qual(i,:)=normrnd(mu,sqrt(var));  %BUG: var=>sqrt(var)
end

%______Main______
%j=5
hitT=0;
for i=1:M,
    for j=1:N,
        if qual(i,j)>c, 
            rnd = rand(1);
            survive=1-exp(-2*(mu(j)-c)*(qual(i,j)-c)/var(j));
            if rnd<survive,
                hitT(i,j)=Inf;   %also has some probability to have finite hitting time
            else
                t=[0:0.1:1000];
                phi=sqrt(t.*tau(j)).*(qual(i,j)-c)+(mu(j)-c)./(var(j).*sqrt(tau(j)*t));
                phi2=sqrt(t*tau(j)).*(qual(i,j)-c)-(mu(j)-c)./(var(j).*sqrt(tau(j)*t));
                y = (1/2)*(1+erf(phi/sqrt(2)))-exp(-2/var(j)*(mu(j)-c)*(qual(j)-c))*(1/2)*(1+erf(phi2/sqrt(2)));
                y=[y(1:min([length(y) find(y<0|y==0,1,'first')-1])) 0];  %if the length is too long, 
                if length(y)>length(t), 
                    t=[t 1000.1];
                else
                    t=t(1:length(y));
                end
                y=unique(y);
                y=-sort(-y');                %eliminate duplicate
                rnd = rand(length(y), 1);
                t=t(length(t)-length(y)+1:length(t))';  %keep matrix length the same. 
                Xtmp = interp1(y, t, rnd','linear')';
                hitT(i,j)=Xtmp(1);
            end

        else
t=[0:0.1:1000];
phi=sqrt(t.*tau(j)).*(qual(i,j)-c)+(mu(j)-c)./(var(j).*sqrt(tau(j)*t));   %x
phi2=sqrt(t*tau(j)).*(qual(i,j)-c)-(mu(j)-c)./(var(j).*sqrt(tau(j)*t));
y = (1/2)*(1+erf(phi/sqrt(2)))-exp(-2/var(j)*(mu(j)-c)*(qual(j)-c))*(1/2)*(1+erf(phi2/sqrt(2)));
%cdf=y;
%figure(j)
%plot(t,cdf);
y=[y(1:min([length(y) find(y<0|y==0,1,'first')-1])) 0];  %if the length is too long, 
if length(y)>length(t), 
    t=[t 1000.1];
else
    t=t(1:length(y));
end
y=unique(y);
y=-sort(-y');                %eliminate duplicate
rnd = rand(length(y), 1);
t=t(length(t)-length(y)+1:length(t))';  %keep matrix length the same. 
Xtmp = interp1(y, t, rnd','linear')';
hitT(i,j)=Xtmp(1);
        end
    end
end


for i=1:M,
d=tau.*hitT(i,:);
v=k.*tau;
Ind=find(hitT(i,:)~=inf);% find the indices that doesn't have infinite hitting time. N
hitTO(i,:)=inf*ones(1,N);    %initiate hitting time, if previous hitting time=inf, keep it, otherwise, set to 0
hitTO(i,Ind)=0;  



dv=d(Ind)./v(Ind);
while isempty(Ind)-1
    dv=d./v;
    [Ma,j]=min(dv(Ind));       %find the indices of minimum element in dv, BUG:the minimum among Ind is equal to the value outside Ind
    j=Ind(j(1));
    hitTO(i,Ind)=hitTO(i,Ind)+dv(j(1));    %update hitting time
    d=d-v*(d(j(1))/v(j(1)));          %update d
    Ind(find(Ind==j(1)))=[];                  %update N
    %iconnect=find(A(i)==1);     %find all players i connect to 
    %k(iconnect)=max(ones(length(iconnect)),k(iconnect)-1); %update k
    %A(:,i)=0;                                              %update graph, cut all i's connection
    %A(i,:)=0;                                              %update graph 
end
end
hitTO;                                                       %output hitting time
%{
tau=0.9;
phi=sqrt(t*tau)*(q-c)+(mu-c)./(sigma^2*sqrt(tau*t));
y=(mu-c)/(sigma^2*sqrt(tau))*t.^(-1.5).*exp(-0.5*phi.^2)/sqrt(2*pi);
for i=1:length(y),
cdf(i)=sum(y(2:i));
end
plot(t,cdf);

tau=[0.01:0.01:1];
for j=1:100,
phi=sqrt(t*tau(j))*(q-c)+(mu-c)./(sigma^2*sqrt(tau(j)*t));
y=(mu-c)/(sigma^2*sqrt(tau(j)))*t.^(-1.5).*exp(-0.5*phi.^2)/sqrt(2*pi);
for i=1:length(y),
cdf(i)=sum(y(2:i));
end
cdfm(j)=cdf(length(cdf));
end
plot(tau,cdfm)
%}
%%%Compute social welfare with scaled hitting time hitTO
%input graph A, 
for i=1:N, 
    for j=1:N,
        if A(i,j)==1,
            W(:,i)=W(:,i)+(qual(:,i)+qual(:,j)-2*c).*(1-exp(-p*min(hitTO(:,i),hitTO(:,j))));
        end
    end
end
%output social welfare W (1000,60). average welfare avgW with respect to
%player i 
avgW=mean(W);
