function [nr,nrho,nS0,nmu,nT,nH,nLL]=MCMC(H,r,rho,S0,N,mu,startT,MaxT,Iter,Every)

Q=log(1:max(H)); FACT=zeros(length(Q)+1,1);
FACT(1+[1:length(Q)])=cumsum(Q);

nr=zeros(Iter,1); nrho=zeros(Iter,1);
nS0=zeros(Iter,1); nmu=zeros(Iter,1);
nT=zeros(Iter,1); nH=zeros(Iter,MaxT);
nLL=zeros(Iter,1);

m=find(H>=0);

while nr(1)<0 || nrho(1)<1 || nS0(1)<0.01 || nS0(1)>0.99 || nmu(1)<0 || nmu(1)>0.5 || nT(1)>=0
    nr(1)=r+randn(1,1)*0.1; nrho(1)=rho+randn(1,1)*1;
    nS0(1)=S0+randn(1,1)*0.1; nmu(1)=mu+randn(1,1)*0.1;
    nT(1)=startT+randn(1,1)*1;
end

BEST=-1e100;

for i=1:Iter
    for j=1:Every
        
        tr=nr(i)+randn(1,1)*0.01;
        trho=nrho(i)+randn(1,1)*0.1;
        tS0=nS0(i)+randn(1,1)*0.01;
        tmu=nmu(i)+randn(1,1)*0.01;
        tT=nT(i)+randn(1,1)*0.1;        
        
        if tr<0 || trho<1 || tS0<0.0 || tS0>1 || tmu<0 || tmu>1 || tT>=0
            LL=-1e100;
        else
            
            I=1; S=round(tS0*N); R=0;
            Gamma=1/trho;
            Beta=(tr+Gamma)/(tS0*N);
            
            [t, pop]=ode45(@Diff_2_1,[tT 0:(MaxT+1)],[S I R],[],[Beta Gamma]);
            S=pop(:,1); I=pop(:,2); R=[pop(:,3)];
            n=R(3:end)-R(2:(end-1));  % number of recoveries in the past day.
            
            if min(n)>0             
                LL=sum(LPoisson(H(m),n(m)*tmu,FACT));    % Likelihood
                
                LL=LL+log(gampdf(tr,2,0.1));   %Prior on r
                LL=LL-(trho/7);                 %Prior on rho
                LL=LL+log(betapdf(tS0,1.5,1.5));     %Prior on S0
                LL=LL+log(betapdf(tmu,1.05,1.5));    %Prior on mu
                
                if LL>BEST+log(rand(1,1))
                    nr(i)=tr; nrho(i)=trho;
                    nS0(i)=tS0; nmu(i)=tmu;
                    nT(i)=tT;  nLL(i)=LL;
                    nH(i,1:MaxT)=n(1:MaxT)*tmu;
                    BEST=LL;
                end
            end
        end
    end

    if i<Iter
        nr(i+1)=nr(i); nrho(i+1)=nrho(i);
        nS0(i+1)=nS0(i); nmu(i+1)=nmu(i);
        nT(i+1)=nT(i); nLL(i+1)=nLL(i);
        nH(i+1,:)=nH(i,:);
    end
    
    
end


end




% Calculates the differential rates used in the integration.
function dPop=Diff_2_1(t,pop, parameter)

beta=parameter(1); gamma=parameter(2);
S=pop(1); I=pop(2); R=pop(3);

dPop=zeros(3,1);

dPop(1)= -beta*S*I;
dPop(2)= beta*S*I - gamma*I;
dPop(3)= gamma*I;

end


function [LL]=LPoisson(D,M,FACT)

M=M(1:length(D));
LL=sum( log(M).*D - M - FACT(D+1) );

end
