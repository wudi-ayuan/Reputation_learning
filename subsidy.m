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
Ni=3:1:4;               %number of players
N=5;
M=2000;             %number of draws, Monte Carlo simulation
c_orig=0;              %cost 
c_delta=0.5;           %subsidy
%mu=[0 0 0 100 100]; %true quality distribution means, should be greater than c
%var=[5 5 5 10 10]; %true quality distribution variance, sigma^2
p=1;              %discount rate 
%base precision, uniform distribution [2,7]
c_gonna_test=c_orig-[linspace(1,14,14) 15:0.05:25];
avgW=zeros(1,length(c_gonna_test));
avgWind=zeros(1,N);

%_____INIT______
tau=1*ones(1,max(N));
qual=ones(M,max(N));
mu=[5 4 3 2 1].*ones(1,max(N));
var=5*ones(1,max(N));
hitcost=ones(1,length(c_gonna_test));

for i=1:M,                       %initiate quality matrix(M*N)
qual(i,:)=normrnd(mu,var);
end

%parpool('local',2);
syms t;
%Random=rand(M,N);
for c=c_gonna_test,
    c
%    hitT=zeros(M,N);
for i=1:M,
    for j=1:max(N),
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
                if 40>br&&br>=20
                    srchR=1;
                elseif 60>br&&br>=40
                    srchR=inf;
                elseif br>=60
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
if 40>br&&br>=20
    srchR=1;
elseif 60>br&&br>=40
    srchR=inf;
elseif br>=60
    break
end
end
hitT(i,j)=Xtmp;
        end
    end
    
    hitcost(e)=mean(hitT(hitT~=inf));
end
    
    %______Main______
%%%%%%%%%%%%%%%%%%%%%%%%% After hitting time %%%%%%%%%%%%%%%%%%%%%%   



    A=ones(N);

    %Graph adjacent matrix A

    B=triu(A,1);
    A=B+B';% undirected graph, symmetric adjacent matrix
    for i=1:N, 
       k(i)=sum(A(i,:));    %player's rank
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  OUTPUT LIST: %%%%%%%%%%%%%%%%%%%%%%%
    W=zeros(M,N);
    hitTO=zeros(M,N);
    %avgW=[];

%%%%%%%%%%%%%%%%%%%%%%%%% scale the hitting time %%%%%%%%%%%%%%%%%%%%%%   
    %% No Parallel 

%    hitTO=inf*ones(M,N);
    for i=1:M,
        for ki=1:N, 
            k(ki)=sum(A(ki,:));    %reset rank
        end
    d=tau(1:N).*hitT(i,:);
    v=k(1:N).*tau(1:N);
    Ind=find(hitT(i,:)~=inf);% find the indices that doesn't have infinite hitting time. N
    hitTO(i,:)=inf*ones(1,N);    %initiate hitting time, if previous hitting time=inf, keep it, otherwise, set to 0
    hitTO(i,Ind)=0;  


    dv=d(Ind)./v(Ind);
    A_scale=A;
    while isempty(Ind)-1
        dv=d./v;
        [Ma,j]=min(dv(Ind));       %find the indices of minimum element in dv, BUG:the minimum among Ind is equal to the value outside Ind
        j=Ind(j(1));
        hitTO(i,Ind)=hitTO(i,Ind)+dv(j(1));    %update hitting time
        d=d-v*(d(j(1))/v(j(1)));          %update d
        Ind(find(Ind==j(1)))=[];                  %update N
        %iconnect=find(A(i)==1);     %find all players i connect to 
        %k(iconnect)=max(ones(length(iconnect)),k(iconnect)-1); %update k
        A_scale(:,j)=0;                                              %update graph, cut all i's connection
        A_scale(j,:)=0;                                              %update graph 
        for ki=1:N, 
        k(ki)=sum(A_scale(ki,:));    %player's new ranks
        end
        v=k.*tau;
        if sum(isnan(hitTO(i,:)))~=0,
            i
        end
    end
    end
    

    %%%Compute social welfare with scaled hitting time hitTO
    %input graph A, 
    wel_add=zeros(M,N);%only about i
    survive=zeros(M,N);%only about i
    fepsilon=zeros(M,N);%only about i
    for j=1:N, 
        survive(:,j)=1-exp(-2*(mu(j)-c)*(qual(:,j)-c)/var(j));%(18) welfare
        survive(qual(:,j)<=c,j)=0; 
    end
    %(9), in subsidy case, add additional term for agent not ostracized. 
    for i=1:N, 
        for j=1:N,
            if A(i,j)==1,
                indek=find(hitT(:,j)==inf);
                wel_add(indek,i)=wel_add(indek,i)+(mu(j)-c_orig)./survive(indek,j);
                %if c~=c_orig,
                %wel_add(hitT(:,j)~=inf,i)=wel_add(hitT(:,j)~=inf,i)+c;
                %end
                reallybad=find(hitT(:,j)<hitT(:,i));
                W(reallybad,i)=W(reallybad,i)+(1-exp(-p*hitTO(reallybad,j)))/p*c;
            end
        end
        W(:,i)=W(:,i)+(1-exp(-p*hitTO(:,i)))/p.*wel_add(:,i);   %(11) welfare
        
    end
    %output social welfare W (1000,60). average welfare avgW with respect to
    %player i 
    avgWind(e)=nanmean(W(:,2));
    avgW(e)=mean(mean(W));
    mean(hitTO(isfinite(hitTO(:,1)),1));
    mean(hitTO(isfinite(hitTO(:,2)),2));

    e=e+1;
end
%delete(gcp)
%%W*
Wstar_add=zeros(1,N);

for i=1:N, 
    for j=1:N,
        if A(i,j)==1,
            Wstar_add(:,i)=Wstar_add(:,i)+(mu(j)-c_orig);
        end
    end
    W_star(:,i)=Wstar_add(:,i)./p;   %(11) welfare

end
avgW_star=mean(mean((W_star)));

muHi=c_orig*ones(1,length(c_gonna_test))-c_gonna_test;
figure(1)
scatter(muHi,avgW)
xlabel('subsidy')
ylabel('social welfare')
title('theorem 6  ')
hold on
plot(muHi,avgW_star*ones(length(muHi)));
[D,I]=max(avgW)