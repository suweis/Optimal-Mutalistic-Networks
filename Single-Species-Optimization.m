%by Samir Suweis. MATLAB CODE
%Please if you use this code please cite: Suweis et al., Nature 500 (449) 2013 

%% DESCRIPTION OF THE PROCESS WE WANT TO SIMULATE
% This algorithm simulates the adpative evolution where each single species
% try to optimize its population. We start with a set of N nodes (species)
% in a random weighted birpartite networks with connectivity C. Then we
% pick a random species, say bee i, and one of its partner, plant j. We then pick
% another plant, l, that may (or may not) interact with that bee, and we swap
% (or rewire) the link from i-j to i-l. In order to seed up the code, what
% we actually will do is to pick two couple of interacting species, say i-j
% and k-l. Then switch the link among them. And finally we accept the swap if both i
% and k (or both j and l), have increased their species population. This is
% two steps equivalent to what said above.
%We then calculate the new stationary state x*', and if x*'_i>x*_i we accept the switch. Otherwise
% we reject it. We repeat the operation for T steps.

%% ALGORITHM

clear all
close all
matlabpool    % Use default parallel configuration. IF you do not want to
%have parallal code comment this and the following line set for instead of parfor

R=1;      %Number of realizations of the process (so to calculate average behaviour)


for r=1:R
    tic
    % Set here  paramters of the model.
    nAvec=[10];% number of animal pollinators - you can set a vector of different size
    nPvec=[10];% number of plants
    beta_vec=[1]; %competition with respect to mutualism, i.e. =1 means competition=mutualism. you can set a vector of different size
    for dim=1:length(nAvec)
        nP=nPvec(dim);
        nA=nAvec(dim);
        N=nP+nA;
        c=4./(N.^(0.8)); % connectivity. Here set as in real data
        d=1; % self interaction.
        
        for bb=1:length(beta_vec)
            beta=beta_vec(bb);
            
            sigma_vec=[0.05]; % mutualistic interaction strengths
            
            for s=1:length(sigma_vec)
                sigma=sigma_vec(s);
                
                xeq=1.*ones(N,1);%initial populations. You can also set it random, i.e. 0.5+(1.5-0.5).*rand(N,1);%
                b=zeros(N,N);% starting interaction matrix
                
                
                %% Fills matrix with direct competition & mutalistic interaction
                for i=1:nA
                    for j=1:nA
                        if i==j
                            b(i,i)=d;
                        else
                            if rand<c
                                b(i,j)=beta*abs(sigma.*randn)*d;
                                b(j,i)=beta*abs(sigma.*randn)*d;
                            end
                        end
                        
                    end
                    for j=nA+1:N
                        if rand<c
                            b(i,j)=-abs(sigma.*randn);
                            b(j,i)=-abs(sigma.*randn);
                        end
                    end
                    
                end
                for i=nA+1:N
                    for j=nA+1:N
                        if i==j
                            b(i,i)=d;
                        else
                            if rand<c
                                b(i,j)=beta*abs(sigma.*randn)*d;
                                b(j,i)=beta*abs(sigma.*randn)*d;
                            end
                        end
                        
                    end
                end
                bold=b;
                b_antA=b(1:nA,1:nA); % competition matrix among pollinators
                b_antP=b(nP+1:N,nP+1:N); % competition matrix among plants
                
                Aold=zeros(N,N);%linearization around equilibrium xeq
                for i=1:N
                    Aold(:,i)=-xeq.*bold(:,i);
                end
                
                alpha=b*xeq;% intrisic grotwh parameter set in order to have positive populations
                T=100^3/2; % total number of time step
                Xtot=zeros(T,1); 
                %% Core of the algorithm. I use logical programming.
                for t=1:T
                    
                    bnewA=b(1:nA,nA+1:N);
                    bnewP=b(nA+1:N,1:nA);
                    
                    ind1=bnewA~=0;
                    allind=1:length(bnewA(:));
                    die=ceil(rand*length(bnewA(ind1)));
                    [i1,j1]=ind2sub(size(bnewA),allind(ind1));
                    
                    ind0=true(size(bnewA));
                    ind0(i1(die),j1(die))=0;
                    new=ceil(rand*length(bnewA(ind0)));
                    [i0, j0]=ind2sub(size(bnewA),allind(ind0));
                    
                    picked1=bnewA(i0(new),j0(new));
                    picked2=bnewP(j0(new),i0(new));
                    bnewA(i0(new),j0(new))= bnewA(i1(die),j1(die));
                    bnewP(j0(new),i0(new))= bnewP(j1(die),i1(die));
                    bnewA(i1(die),j1(die))=picked1;
                    bnewP(j1(die),i1(die))=picked2;
                    bnew=[b_antA, bnewA; bnewP,b_antP];
                    
                    xeqE=bnew\alpha;
                    
                    bee1=i0(new)+nP;
                    bee2=i1(die)+nP;%bees involved in the swap
                    plant1=j0(new);
                    plant2=j1(die);%plant involved in the swap
                    spI=[bee1 plant1 bee2 plant2];
                    rI=randi(2,1);
                    spR=[spI(rI) spI(rI+2)];
                    
                    if sum(xeqE>0)-length(xeqE)==0 && xeqE(spR(1))>=xeq(spR(1)) &&  xeqE(spR(2))>=xeq(spR(2))%spR(1)==spR(2) if you wish to constrain that species i and k are the same
                        xeq=xeqE;
                        b=bnew;
                    end
                    Xtot(t)=sum(xeq);
                end
                A=zeros(N,N);%linearization around equilibrium xeq
                for i=1:N
                    A(:,i)=-xeq.*b(:,i);
                end
                saveb(sprintf('_r=%d',r),bold,b,A,N,sigma_vec(s),beta);
            end
        end
    end
    toc
end
matlabpool close
