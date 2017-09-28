%% This section aim to investigate the core-periphery networks. Agents are divided into 2 gropus. 
%%
close all

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
avgW=ones(2,5);
diff32=[];
diff64=[];
iscore_peri=[];
e=1;
N=5;               %number of players
M=3000;             %number of draws, Monte Carlo simulation
c=0;              %cost 
mui=1:0.5:10; %true quality distribution means, should be greater than c
var=[5 5 5 10 10]; %true quality distribution variance, sigma^2
p=1;              %discount rate 
%base precision, uniform distribution [2,7]
tau=ones(1,N);
qual=ones(M,N);
for muH=mui,
    e=1;   %new round
    mu=[1 1 1 10^muH 10^muH];
    avgW=ones(1,N);

for x=1:2, 

    if x==1, 
        A=[0 0 0 1 1; 0 0 0 1 1; 1 1 0 1 1; 1 1 1 0 1; 1 1 1 1 0];
    else 
        A=[0 1 1 1 1; 0 0 1 1 1; 0 0 0 1 1; 0 0 0 0 1; 0 0 0 0 0];
    end

    %Graph adjacent matrix A

    B=triu(A,1);
    A=B+B';% undirected graph, symmetric adjacent matrix
    for i=1:N, 
       k(i)=sum(A(i,:));    %player's rank
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  OUTPUT LIST: %%%%%%%%%%%%%%%%%%%%%%%
    W=zeros(M,N);
    %avgW=[];

    %_____INIT______
    for i=1:M,                       %initiate quality matrix(M*N)
    qual(i,:)=normrnd(mu,var);
    end

    %______Main______

    hitT=zeros;
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
                    y = (1/2)*(1+erf(phi/sqrt(2)))-exp(-2/var(j)*(mu(j)-c)*(qual(i,j)-c))*(1/2)*(1+erf(phi2/sqrt(2)));
                    y=[y(1:min([length(y) find(y<0|y==0,1,'first')-1])) 0];  %if the length is too long, 
                    if length(y)>length(t), 
                        t=[t 1000.1];
                    else
                        t=t(1:length(y));
                    end
                    y=unique(y);
                    y=-sort(-y');                %eliminate duplicate
                    rnd = survive+(1-survive)*rand(length(y), 1);   %since in this case, the agent will be sotracize in a finite time.
                                                                    %The time can be long, but should be finite. Hence it should be drawn evenly from survive to 1        
                    t=t(length(t)-length(y)+1:length(t))';  %keep matrix length the same. 
                    Xtmp = interp1(y, t, rnd','linear')';
                    hitT(i,j)=Xtmp(1);
                end

            else
    t=[0:0.1:1000];
    phi=sqrt(t.*tau(j)).*(qual(i,j)-c)+(mu(j)-c)./(var(j).*sqrt(tau(j)*t));
    phi2=sqrt(t*tau(j)).*(qual(i,j)-c)-(mu(j)-c)./(var(j).*sqrt(tau(j)*t));
    y = (1/2)*(1+erf(phi/sqrt(2)))-exp(-2*(mu(j)-c)*(qual(i,j)-c)/var(j))*(1/2)*(1+erf(phi2/sqrt(2)));
    y=[y(1:min([length(y) find(y<0|y==0,1,'first')-1])) 0];  %if the length is too long, 
    if length(y)>length(t), 
        t=[t 1000.1];
    else
        t=t(1:length(y));
    end
    y=unique(y);
    y=-sort(-y');                %eliminate duplicate
    rnd =min(y)+(1-min(y))*rand(length(y), 1);
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
                %W(:,i)=W(:,i)+(qual(:,j)-c).*(1-exp(-p*min(hitTO(:,i),hitTO(:,j))))/p;
                %W(:,i)=W(:,i)+(qual(:,j)+qual(:,i)-2*c).*(1-exp(-p*min(hitTO(:,i),hitTO(:,j))))/p;
                %W(:,i)=W(:,i)+(qual(:,i)+qual(:,j)-2*c).*(1-exp(-p*min(hitTO(:,i),hitTO(:,j))));
                wel_add(indek,i)=wel_add(indek,i)+(mu(j)-c)./survive(indek,j);
            end
        end
        %phiinternal=sqrt(hitTO(:,i)*tau(i)).*(qual(:,i)-c)+(mu(i)-c)./(var(i)*sqrt(hitTO(:,i)*tau(i)));
        %phif=exp(-phiinternal.^2)/sqrt(pi);
        %fepsilon=((mu(i)-c)/(var(i)*sqrt(tau(i)))*hitTO(:,i).^(-1.5)).*phif;
        W(:,i)=(1-exp(-p*hitTO(:,i)))/p.*wel_add(:,i);   %(11) welfare

    end
    %output social welfare W (1000,60). average welfare avgW with respect to
    %player i 
    avgW(e,:)=mean(W);
    Agoal=[0 0 0 1 1; 0 0 0 1 1; 1 1 0 1 1; 1 1 1 0 1; 1 1 1 1 0];
    if sum(sum(A==Agoal))==16, %note down the index which network structure is core-periphery
        avgWI(e)=1;
    else
        avgWI(e)=0;
    end
    if e==1, 
        W_32=W;
        hitT_32=hitTO; 
        row2_32=hitT_32(:,2);
        row1_32=hitT_32(:,1);
        ['the realization survive CP ',num2str((M-sum(isfinite(row2_32)))/(M))]
        ['the survive probaility CP from prop1 ',num2str(mean(survive(:,2))) ]
        diff32=[diff32 abs((M-sum(isfinite(row2_32)))/(M)-mean(survive(:,2)))];
    elseif e==2, 
        W_64=W;
        hitT_64=hitTO;
        row2_64=hitT_64(:,2);
        row1_64=hitT_64(:,1);
        ['the realization survive COMP ',num2str((M-sum(isfinite(row2_64)))/(M))]
        ['the survive probaility COMP from prop1 ',num2str(mean(survive(:,2))) ]
        diff64=[diff64 abs((M-sum(isfinite(row2_32)))/(M)-mean(survive(:,2)))];

    end
    e=e+1;
    
    end
[D,I]=max(mean(avgW,2));
iscore_peri=[iscore_peri I];%note down the maximum 
%sum(isfinite(row1_32))+sum(isfinite(row2_32))
%sum(isfinite(row1_64))+sum(isfinite(row2_64))
end
figure(1)
scatter(mui,iscore_peri)
xlabel('muH')
ylabel('optimal network is core-periphery')
title('theorem 5 ')
