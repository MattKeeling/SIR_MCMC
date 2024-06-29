%% Show Parameters at Four Time points

clear

FigNum=0;

for N=1e6 %[1e5 1e6 1e7]
    
    FigNum=FigNum+1;
    r=0.2; rho=7; S0=0.5; mu=0.1; MaxT=100;
    Iter=1e4; Every=1000; k=0;
   
    %% Generate Simulation
    MX=0; count=0;
    while MX<1e-4
        if count>0 fprintf(1,'Re-simulating %d\n',count+1); end
        count=count+1;
        h=Generate_Simulations(r,rho,S0,N,mu, MaxT);
        MX=max(h)/N;
    end
    SM=movmean(h,7);   %Some finessing to get the approximate peak at 50 days.
    [y i]=max(SM);
    k=k+1; maxH(k)=y;
    H=zeros(MaxT,1);
    for t=1:MaxT
        j=t+i-50;
        if j<1 | j>MaxT
            H(t)=-1;
        else
            H(t)=h(j);
        end
    end
    figure(FigNum); clf;
    semilogy(1:MaxT,H); drawnow;
    
    %% Run MCMC
    
    Times=[20 30 40 50 60]; l=length(Times);
    
    Col=get(gca,'ColorOrder');
    set(gcf,'position',[62 140 1226 1197]); clf;
    
    for tt=1:l
        
        for LOOP=1:6
            
            
            [nr,nrho,nS0,nmu,nT,nH,nLL]=MCMC_mex(H(1:Times(tt)),r,rho,S0,N,mu,0,MaxT,Iter,Every);
            % Note - for optimal performance MCMC should be compiled and
            % replaced with MCMC_mex file.
            
            figure(1);
            subplot(5,l,tt);
            m=unique(round(length(nH)*[0.5:0.01:1]));
            h=plot(1:MaxT,nH(m,:)','-','Color',1-0.4*(1-Col(LOOP,:))); hold on;
            
            m=unique(round(length(nH)*[0.5:0.001:1]));
            subplot(5,l,tt+l);
            [f,xi]=ksdensity(nr(round(length(r)/2):end)); f=f/sum(f);
            plot(xi,f,'-','LineWidth',2,'Color',Col(LOOP,:)); hold on
            
            subplot(5,l,tt+2*l);
            plot(nS0(m),nmu(m),'.k','MarkerSize',4,'Color',Col(LOOP,:)); hold on
            
            subplot(5,l,tt+3*l);
            plot(nrho(m),nmu(m),'.k','MarkerSize',4,'Color',Col(LOOP,:)); hold on
            
            subplot(5,l,tt+4*l);
            plot(nrho(m),nS0(m),'.k','MarkerSize',4,'Color',Col(LOOP,:)); hold on
            drawnow;
        end
        
        figure(FigNum);
        subplot(5,l,tt); set(gca,'FontSize',12);
        h=plot(1:MaxT,H,'-k',1:Times(tt),H(1:Times(tt)),'-r'); hold off;
        h(1).LineWidth=2; h(end).LineWidth=3;
        xlabel('Time','FontSize',14); ylabel('Hospital Admissions','FontSize',14);
        
        subplot(5,l,tt+l); set(gca,'FontSize',12);
        set(gca,'XLim',[0.1 0.4]);
        xlabel('Growth rate, r','FontSize',14);
        
        subplot(5,l,tt+2*l); set(gca,'FontSize',12);
        X=get(gca,'XLim'); X(1)=min(X(1),0.4); X(2)=max(X(2),0.6);
        Y=get(gca,'YLim'); Y(1)=min(Y(1),0.0); Y(2)=max(Y(2),0.2);
        X=[0 1]; Y=[0 0.5];
        plot(X,[0 0]+mu,'-k',[0 0]+S0,Y,'-k','LineWidth',1); hold off
        axis([X Y]);
        xlabel('S_0','FontSize',14); ylabel('\mu','FontSize',14);
        
        subplot(5,l,tt+3*l); set(gca,'FontSize',12);
        X=get(gca,'XLim'); X(1)=min(X(1),0); X(2)=max(X(2),10);
        Y=get(gca,'YLim'); Y(1)=min(Y(1),0.0); Y(2)=max(Y(2),0.2);
        X=[0 15]; Y=[0 0.5];
        plot(X,[0 0]+mu,'-k',[0 0]+rho,Y,'-k','LineWidth',1); hold off
        axis([X Y]);
        xlabel('\rho'); ylabel('\mu');
        
        subplot(5,l,tt+4*l); set(gca,'FontSize',12);
        X=get(gca,'XLim'); X(1)=min(X(1),0); X(2)=max(X(2),10);
        Y=get(gca,'YLim'); Y(1)=min(Y(1),0.4); Y(2)=max(Y(2),0.4);
        X=[0 15]; Y=[0 1];
        plot(X,[0 0]+S0,'-k',[0 0]+rho,Y,'-k','LineWidth',1); hold off
        axis([X Y]);
        xlabel('\rho','FontSize',14); ylabel('S_0','FontSize',14);
        drawnow;
    end
    
    for tt=1:l
        subplot(5,l,tt);
        title(['T=' num2str(Times(tt))],'FontSize',14);
    end
    
    print(['Fits_' num2str(N) '.png'],'-dpng');
end






