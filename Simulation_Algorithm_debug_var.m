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
avgW=ones(11,4);
whatmaxdegree=[];
e=1;
N=4;               %number of players
M=1000;             %number of draws, Monte Carlo simulation
c=0;              %cost 
mu=1*ones(1,N); %true quality distribution means, should be greater than c
var=2*ones(1,N); %true quality distribution variance, sigma^2
p=100;              %discount rate 
%base precision, uniform distribution [2,7]
tau=ones(1,N);
for p=1:0.3:4,
    e=1;
for x=1:1:N,
    A=zeros(N,N);
    A(1,:)=[ones(1,x) zeros(1,N-x)];
%during the trial, the parameters should be carefully decided. 
%G(N,sum(k))
    for y=0:1:2^(N-2)-1,
        A(2,N-1)=fix(y/2);
        A(2,N)=mod(y,2);
        for z=1:1:2^(N-3),
            A(3,N)=fix(z/2);

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
avgW(e,:)=mean(W);
avgWI(e)=sum(k)/2;
e=e+1;
        end
    end
end

%figure(p)
%scatter(avgWI,mean(avgW,2))
%xlabel('totall number of edges')
%ylabel('social welfare')
%title(['social welfare to different number of edges, given discountrate=', num2str(p)])
[D,I]=max(mean(avgW,2));
whatmaxdegree=[whatmaxdegree avgWI(I)];
end

