delete(gcp)
%close all
%clear all

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

for i=1:M,                       %initiate quality matrix(M*N)
qual(i,:)=normrnd(mu,var);
end

parpool('local',2);
syms t;

%    hitT=zeros(M,N);
parfor i=13335:M,
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
 