%% This section aim to investigate the core-periphery networks. Agents are divided into 2 gropus. 
%%
%delete(gcp)
close all
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
e=1;
Ni=3:1:10;               %number of players
M=20000;             %number of draws, Monte Carlo simulation
c=0;              %cost 
%mu=[0 0 0 100 100]; %true quality distribution means, should be greater than c
%var=[5 5 5 10 10]; %true quality distribution variance, sigma^2
p=1;              %discount rate 
%base precision, uniform distribution [2,7]
avgW=zeros(1,length(Ni));
avgWind=zeros(1,length(Ni));

%_____INIT______
tau=ones(1,max(Ni));
qual=ones(M,max(Ni));
mu=1*ones(1,max(Ni));
var=20*ones(1,max(Ni));
%%
for i=1:M,                       %initiate quality matrix(M*N)
qual(i,:)=normrnd(mu,var);
end

parpool('local',2);
syms t;

%    hitT=zeros(M,N);
parfor i=1:M,
    for j=1:10,
        if qual(i,j)>c, 
            rnd = rand(1);
            survive=1-exp(-2*(mu(j)-c)*(qual(i,j)-c)/var(j));
            if rnd<survive,
            hitT(i,j)=Inf;   %also has some probability to have finite hitting time
            else
            phi=sqrt(t.*tau(j)).*(qual(i,j)-c)+(mu(j)-c)./(var(j).*sqrt(tau(j)*t));
            phi2=sqrt(t*tau(j)).*(qual(i,j)-c)-(mu(j)-c)./(var(j).*sqrt(tau(j)*t));
            y = (1/2)*(1+erf(phi/sqrt(2)))-exp(-2/var(j)*(mu(j)-c)*(qual(i,j)-c))*(1/2)*(1+erf(phi2/sqrt(2)));

            Xtmp=[];
            br=0;
            srchR=10;
            while isempty(Xtmp)
                rnd = survive+(1-survive)*rand(1);   %since in this case, the agent will be sotracize in a finite time.
                                                            %The time can be long, but should be finite. Hence it should be drawn evenly from survive to 1     
                Xtmp = double(vpasolve(y==rnd,t,[0,srchR],'random',true));
                br=br+1;
                if br>20
                    srchR=1;
                elseif br>40
                    srchR=inf;
                elseif br>60
                    break
                end
            end
                hitT(i,j)=Xtmp;
            end

        else
phi=sqrt(t.*tau(j)).*(qual(i,j)-c)+(mu(j)-c)./(var(j).*sqrt(tau(j)*t));
phi2=sqrt(t.*tau(j)).*(qual(i,j)-c)-(mu(j)-c)./(var(j).*sqrt(tau(j)*t));
y = (1/2)*(1+erf(phi/sqrt(2)))-exp(-2*(mu(j)-c)*(qual(i,j)-c)/var(j))*(1/2)*(1+erf(phi2/sqrt(2)));
Xtmp=[];
br=0;
srchR=10;
while isempty(Xtmp)
rnd =rand(1);
Xtmp = double(vpasolve(y==rnd,t,[0,srchR],'random',true));
br=br+1;
if br>20
    srchR=1;
elseif br>40
    srchR=inf;
elseif br>60
    break
end
end
hitT(i,j)=Xtmp;
        end
    end
end
    
    %______Main______
%%%%%%%%%%%%%%%%%%%%%%%%% After hitting time %%%%%%%%%%%%%%%%%%%%%%   
%% 

for N=Ni,

    A=zeros(N);
    for x=1:N, %change the network
       if x==1, 
        A(x,x+1)=1;
        A(x,N)=1;
        elseif x==N, 
        A(x,1)=1;
        A(x,x-1)=1;
        else
        A(x,x+1)=1;
        A(x,x-1)=1;
       end 
    end
    %Graph adjacent matrix A

   B=triu(A,1);
   A=B+B';% undirected graph, symmetric adjacent matrix
    for i=1:N, 
       k(i)=sum(A(i,:));    %player's rank
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  OUTPUT LIST: %%%%%%%%%%%%%%%%%%%%%%%
    W=zeros(M,N);
    hitTO=zeros(M,N);

%%%%%%%%%%%%%%%%%%%%%%%% scale the hitting time %%%%%%%%%%%%%%%%%%%%%%   
    % No Parallel 

   hitTO=inf*ones(M,N);
    for i=1:M,
    A_scale=A;
    for ki=1:N, 
        k(ki)=sum(A_scale(ki,:));    %player's new ranks
    end
    d=tau(1:N).*hitT(i,1:N);
    v=k(1:N).*tau(1:N);
    Ind=find(hitT(i,1:N)~=inf);% find the indices that doesn't have infinite hitting time. N
    hitTO(i,:)=inf*ones(1,N);    %initiate hitting time, if previous hitting time=inf, keep it, otherwise, set to 0
    hitTO(i,Ind)=0;  


    dv=d(Ind)./v(Ind);
    while isempty(Ind)-1
        dv=d./v;
        [Ma,j]=min(dv(Ind));       %find the indices of minimum element in dv, BUG:the minimum among Ind is equal to the value outside Ind
        j=Ind(j(1));
        hitTO(i,Ind)=hitTO(i,Ind)+dv(j(1));    %update hitting time
        if sum(isnan(hitTO(i,:)))~=0,
            i;
        end
        d=d-v*(d(j(1))/v(j(1)));          %update d
        Ind(find(Ind==j(1)))=[];                 %update N
        A_scale(:,j)=0;                          %update graph, cut all i's connection
        A_scale(j,:)=0;                          %update graph 
        for ki=1:N, 
        k(ki)=sum(A_scale(ki,:));    %player's new ranks
        end
        v=k(1:N).*tau(1:N);
        if sum(v)==0 && isempty(Ind)-1,
            v;
        end
        if sum(k(Ind))==0,
            for kd=Ind, 
                hitTO(i,kd)=max(hitTO(i,A(kd,:)==1));
            end
            break
        end
    end
    end
    
    
    %%Compute social welfare with scaled hitting time hitTO
    %input graph A, 
    wel_add=zeros(M,N);%only about i
    survive=zeros(M,N);%only about i
    fepsilon=zeros(M,N);%only about i
    for j=1:N, 
        survive(:,j)=1-exp(-2*(mu(j)-c)*(qual(:,j)-c)/var(j));%(18) welfare
        survive(qual(:,j)<=c,j)=0; 
    end
    for i=1:N, 
        for j=1:N,
            if A(i,j)==1,
                indek=find(hitT(:,j)==inf);
                wel_add(indek,i)=wel_add(indek,i)+(mu(j)-c)./survive(indek,j);
            end
        end
        W(:,i)=(1-exp(-p*hitTO(:,i)))/p.*wel_add(:,i);   %(11) welfare

    end
    %output social welfare W (1000,60). average welfare avgW with respect to
    %player i 
    avgWind(e)=mean(W(:,2));
    avgW(e)=mean(mean(W));
    mean(hitTO(isfinite(hitTO(:,1)),1));
    mean(hitTO(isfinite(hitTO(:,2)),2));

    e=e+1;
end
%delete(gcp)
muHi=Ni;
figure(1)
plot(muHi,avgW)
xlabel('ring size')
ylabel('social welfare of ring network')
title('proposition 7 ')
[D,I]=max(avgW)