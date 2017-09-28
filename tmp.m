e=1;
Ni=3:1:10;               %number of players
N=max(Ni);
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
%qual=ones(M,max(Ni));
%qual=qualout(1:M,:);
%hitT=hitTout(1:M,:);
mu=1*ones(1,max(Ni));
var=20*ones(1,max(Ni));

%for y=1:1:N,  %star network
for N=Ni,       %ring network

    
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
%     A=zeros(N);
%     for x=1:N, %change the network
%        if x==y, 
%         A(x,:)=ones(1,N);
%         A(x,x)=0;
%         else
%         A(x,:)=zeros(1,N);
%         A(x,y)=1;
%        end 
%     end
    %%%%%%%%%%%%%%
    %Graph adjacent matrix A

%    B=triu(A,1);
%    A=B+B';% undirected graph, symmetric adjacent matrix
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  OUTPUT LIST: %%%%%%%%%%%%%%%%%%%%%%%
    W=zeros(M,N);
    hitTO=zeros(M,N);
    %avgW=[];

%%%%%%%%%%%%%%%%%%%%%%%%% scale the hitting time %%%%%%%%%%%%%%%%%%%%%%   
    %% No Parallel 

%    hitTO=inf*ones(M,N);
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
        d=d-v*(d(j(1))/v(j(1)));          %update d
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
                wel_add(indek,i)=wel_add(indek,i)+(mu(j)-c)./survive(indek,j);
            end
        end
        W(:,i)=(1-exp(-p*hitTO(:,i)))/p.*wel_add(:,i);   %(11) welfare

    end
    %output social welfare W (1000,60). average welfare avgW with respect to
    %player i 
    avgWind(e)=nanmean(W(:,2));
    avgW(e)=nanmean(mean(W));
    mean(hitTO(isfinite(hitTO(:,1)),1));
    mean(hitTO(isfinite(hitTO(:,2)),2));

    e=e+1;
end
muHi=Ni;
figure(1)
plot(muHi,avgWind)
xlabel('ring size')
ylabel('social welfare of ring network')
title('proposition 7 ')
[D,I]=max(avgWind)