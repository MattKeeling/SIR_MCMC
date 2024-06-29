function [H]=Generate_Simulations(r,rho,S0,N,mu, MaxT)

% Gellespie algorithm for simple SIR model. Calculating the number of
% hospital admissions in each one-day time period. 

gamma=1/rho;
beta=(r+gamma)/(S0*N);

I=ceil(N/1e5); S=round(S0*N);
H=zeros(MaxT,1);

h=0; tlast=0; t=0;
while t<MaxT
    Rates=beta*S*I + gamma*I;
    if Rates==0
        t=MaxT+1;
    else
        t=t-log(rand(1,1))/Rates;
        if t>=tlast+1 % gone past one day  - this might seem an odd way of doing things, but we want to catch it before the latest hospitalisation is added.
            H(tlast+1)=random('binomial',h,mu); h=0;
            tlast=floor(t);
        end
        if rand(1,1)<beta*S*I/Rates
            I=I+1; S=S-1;
        else     
            I=I-1; h=h+1;
        end
    end
end

H=H(1:MaxT);
