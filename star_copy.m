%% This section aim to investigate the core-periphery networks. Agents are divided into 2 gropus. 
%%
%close all
%clear all
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
%delete(gcp)

e=1;
Ni=1:1:10;               %number of players
M=2000;             %number of draws, Monte Carlo simulation
c=0;              %cost 
%mu=[0 0 0 100 100]; %true quality distribution means, should be greater than c
%var=[5 5 5 10 10]; %true quality distribution variance, sigma^2
p=1;              %discount rate 
%base precision, uniform distribution [2,7]
N=6;
tau=ones(M,N);
qual=ones(M,N);
mu=ones(M,N);
var=2*ones(M,N);
muCI=1;
varCI=1;

A=zeros(N);
A(:,1)=ones(1,N);
A(1,:)=ones(1,N);
A(1,1)=0; 
varCRang=1:0.1:10;
varl=length(varCRang);
muCRang=1:0.1:10;
mul=length(muCRang);
%parpool('local',2);

syms T Tau Qual Mu Var Phi Phi2 t
t=sym('t',[M,N]);
%%%%%%define function used to find inverse of hitting time
phi(T,Tau,Qual,Mu,Var)=sqrt(T.*Tau).*(Qual-c)+(Mu-c)./(Var.*sqrt(Tau*T));
phi2(T,Tau,Qual,Mu,Var)=sqrt(T*Tau).*(Qual-c)-(Mu-c)./(Var.*sqrt(Tau*T));
y(Phi,Phi2,Mu,Qual,Var) = (1/2)*(1+erf(Phi/sqrt(2)))-exp(-2*(Mu-c)*(Qual-c)/Var)*(1/2)*(1+erf(Phi2/sqrt(2)));
                        
avgW=[];
for muC=1:mul,
    mu=ones(M,N);
    mu(:,1)=muCRang(muC);
    for varC=1:varl,
        var=2*ones(M,N);
        var(:,1)=varCRang(varC);
    %Graph adjacent matrix A

%    B=triu(A,1);
%    A=B+B';% undirected graph, symmetric adjacent matrix
    k=sum(A);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  OUTPUT LIST: %%%%%%%%%%%%%%%%%%%%%%%
    W=zeros(M,N);
    %avgW=[];

    %_____INIT______
    qual=normrnd(mu,var);

    %______Main______

    hitT=zeros(M,N);
    hitTO=zeros(M,N);
    survive=1-exp(-2*(mu-c).*(qual-c)./var);
    rnd=rand(M,N);
    %%%%%%three classes
    h_s=find(qual>c&rnd<survive);
    h_n_s=find(qual>c&rnd>=survive);
    low=find(qual<=c);
    %first class: survive
    hitT(h_s)=inf;
    %second class: good but not survive
    rnd_hit = survive(h_n_s)+(1-survive(h_n_s)).*rand(size(h_n_s));
    search=zeros(length(h_n_s),2);
    search(:,2)=inf;
    %set the search region
    pre_guess=double(y(phi(10,tau(h_n_s),qual(h_n_s),mu(h_n_s),var(h_n_s)),phi2(10,tau(h_n_s),qual(h_n_s),mu(h_n_s),var(h_n_s)),mu(h_n_s),qual(h_n_s),var(h_n_s)));
    search(rnd_hit>pre_guess,2)=10;
    pre_guess=double(y(phi(1,tau(h_n_s),qual(h_n_s),mu(h_n_s),var(h_n_s)),phi2(1,tau(h_n_s),qual(h_n_s),mu(h_n_s),var(h_n_s)),mu(h_n_s),qual(h_n_s),var(h_n_s)));
    search(rnd_hit>pre_guess,2)=1;
    Xtmp = vpasolve(y(phi(t(h_n_s),tau(h_n_s),qual(h_n_s),mu(h_n_s),var(h_n_s)),phi2(t(h_n_s),tau(h_n_s),qual(h_n_s),mu(h_n_s),var(h_n_s)),mu(h_n_s),qual(h_n_s),var(h_n_s))==rnd_hit,t(h_n_s),search,'random',true);
    for ti=h_n_s',
        hitT(ti)=getfield(Xtmp,char(t(ti)));
    end

    for i=1:M,
        for j=1:N,
            if qual(i,j)>c, 
                rnd = rand(1);
                survive=1-exp(-2*(mu(j)-c)*(qual(i,j)-c)/var(j));
                if rnd<survive,
                    hitT(i,j)=Inf;   %also has some probability to have finite hitting time
                else
                phi(t)=sqrt(t.*tau(j)).*(qual(i,j)-c)+(mu(j)-c)./(var(j).*sqrt(tau(j)*t));
                phi2(t)=sqrt(t*tau(j)).*(qual(i,j)-c)-(mu(j)-c)./(var(j).*sqrt(tau(j)*t));
                y(phi(t),phi2(t)) = (1/2)*(1+erf(phi(t)/sqrt(2)))-exp(-2*(mu(j)-c)*(qual(i,j)-c)/var(j))*(1/2)*(1+erf(phi2(t)/sqrt(2)));
                    y=[y(1:min([length(y) find(y<0|y==0,1,'first')-1])) 0];  %if the length is too long, 
                    if length(y)>length(t), 
                        t=[t 7.1];
                    else
                        t=t(1:length(y));
                    end
                    y=unique(y);
                    y=-sort(-y');                %eliminate duplicate
                    rnd = survive+(1-survive)*rand(length(y), 1);   %since in this case, the agent will be sotracize in a finite time.
                                                                    %The time can be long, but should be finite. Hence it should be drawn evenly from survive to 1        
                    t=t(length(t)-length(y)+1:length(t))';  %keep matrix length the same. 
                    Xtmp = interp1(y, t, rnd','pchip')';
                    hitT(i,j)=Xtmp(1);
                end

            else
    %t=[0:0.00001:2];
    phi(t)=sqrt(t.*tau(j)).*(qual(i,j)-c)+(mu(j)-c)./(var(j).*sqrt(tau(j)*t));
    phi2(t)=sqrt(t*tau(j)).*(qual(i,j)-c)-(mu(j)-c)./(var(j).*sqrt(tau(j)*t));
    y(phi(t),phi2(t)) = (1/2)*(1+erf(phi(t)/sqrt(2)))-exp(-2*(mu(j)-c)*(qual(i,j)-c)/var(j))*(1/2)*(1+erf(phi2(t)/sqrt(2)));
    y=[y(1:min([length(y) find(y<0|y==0,1,'first')-1])) 0];  %if the length is too long, 
    if length(y)>length(t), 
        t=[t 2.1];
    else
        t=t(1:length(y));
    end
    y=unique(y);
    y=-sort(-y');                %eliminate duplicate
    rnd =min(y)+(1-min(y))*rand(length(y), 1);
    t=t(length(t)-length(y)+1:length(t))';  %keep matrix length the same. 
    Xtmp = interp1(y, t, rnd','pchip')';
    hitT(i,j)=Xtmp(1);
            end
        end
    end


    for i=1:M,
    for ki=1:N, 
       k(ki)=sum(A(ki,:));    %player's rank
    end
    d=tau(1:N).*hitT(i,1:N);
    v=k(1:N).*tau(1:N);
    Ind=find(hitT(i,1:N)~=inf);% find the indices that doesn't have infinite hitting time. N
    hitTO(i,:)=inf*ones(1,N);    %initiate hitting time, if previous hitting time=inf, keep it, otherwise, set to 0
    hitTO(i,Ind)=0;  


    dv=d(Ind)./v(Ind);
    A_scale=A;
    while isempty(Ind)-1
        dv=d./v;
        [Ma,j]=min(dv(Ind));       %find the indices of minimum element in dv, BUG:the minimum among Ind is equal to the value outside Ind
        j=Ind(j(1));
        hitTO(i,Ind)=hitTO(i,Ind)+dv(j(1));    %update hitting time
        if sum(isnan(hitTO(i,:)))~=0,
            i;
        end
        d=d-v*(d(j(1))/v(j(1))) ;         %update d
        Ind(find(Ind==j(1)))=[];                  %update N
        A_scale(:,j)=0;                                              %update graph, cut all i's connection
        A_scale(j,:)=0;                                              %update graph 
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

    for i=1:N, 
        wel_add=zeros(M,1);%only about i
        survive=zeros(M,1);%only about i
        fepsilon=zeros(M,1);%only about i
        for j=1:N,
            survive(:,1)=1-exp(-2*(mu(j)-c)*(qual(:,j)-c)/var(j));%(11) welfare
            if A(i,j)==1,
                indek=find(hitT(:,j)==inf);
                %W(:,i)=W(:,i)+(qual(:,j)-c).*(1-exp(-p*min(hitTO(:,i),hitTO(:,j))))/p;
                %W(:,i)=W(:,i)+(qual(:,j)+qual(:,i)-2*c).*(1-exp(-p*min(hitTO(:,i),hitTO(:,j))))/p;
                %W(:,i)=W(:,i)+(qual(:,i)+qual(:,j)-2*c).*(1-exp(-p*min(hitTO(:,i),hitTO(:,j))));
                wel_add(indek,1)=wel_add(indek,1)+(mu(j)-c)./survive(indek,1);
            end
        end
        %phiinternal=sqrt(hitTO(:,i)*tau(i)).*(qual(:,i)-c)+(mu(i)-c)./(var(i)*sqrt(hitTO(:,i)*tau(i)));
        %phif=exp(-phiinternal.^2)/sqrt(pi);
        %fepsilon=((mu(i)-c)/(var(i)*sqrt(tau(i)))*hitTO(:,i).^(-1.5)).*phif;
        W(:,i)=(1-exp(-p*hitTO(:,i)))/p.*wel_add;   %(11) welfare

    end
    %output social welfare W (1000,60). average welfare avgW with respect to
    %player i 
    avgW(muC,varC)=nanmean(mean(W),2);


        
    end
   

end
figure(1)
contourf(muRang,varCRang,avgW)
xlabel('Mu')
ylabel('Var')
zlabel('social welfare')
title('theorem 5 ')
