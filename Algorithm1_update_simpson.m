clear all
%_____I/O______
%input list:
%N, E, Graph A, base precision, hitting time hitT
%G(N,E)
N=5;               %number of player
%E=1000;             %number of links
%base precision, uniform distribution
tau_max=50;
tau_min=2;
tau=randi([tau_min,tau_max],1,N);
%Graph adjacent matrix A
A=randi([0,1],N,N);
B=triu(A,1);
A=B+B';% undirected graph, symmetric adjacent matrix
for i=1:N, 
   k(i)=sum(A(i,:));    %player's rank
end
%hitting time without network effect, uniform distribution
hit_max=900;
hit_min=50;
hitT=randi([hit_min,hit_max],1,N);
hitTc=unidrnd(2,1,N)-1;
hitT=hitT./hitTc;

%output list:
%ouput hitting time
hitTO=[];                                                       %output hitting time


%_____INIT______
d=hitT.*tau;
v=k.*tau;
Ind=find(hitT~=inf);% find the indices that doesn't have infinite hitting time. N
hitTO=inf*ones(1,N);    %initiate hitting time, if previous hitting time=inf, keep it, otherwise, set to 0
hitTO(Ind)=0;  

%______Main______
dv=d(Ind)./v(Ind);
while isempty(Ind)-1
    dv=d./v;
    i=find(dv==min(dv(Ind)));       %find the indices of minimum element in dv
    hitTO(Ind)=hitTO(Ind)+dv(i);    %update hitting time
    d=d-v*(d(i)/v(i));          %update d
    Ind(find(Ind==i))=[];                  %update N
    iconnect=find(A(i)==1);     %find all players i connect to 
    %k(iconnect)=max(ones(length(iconnect)),k(iconnect)-1); %update k
    A(:,i)=0;                                              %update graph, cut all i's connection
    A(i,:)=0;                                              %update graph 
    for i=1:N, 
        k(i)=sum(A(i,:));    %player's new ranks
    end
    v=k.*tau;
end
hitTO                                                       %output hitting time