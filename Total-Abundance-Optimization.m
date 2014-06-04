% By Samir Suweis. If you use the code please cite: Suweis et al., Nature 500 (449), 2013

% For detail on the code please look at the Single Species Optimization code.

clear all
close all
matlabpool % Use default parallel configuration

R=100;

parfor r=1:R
    tic
    nAvec=[25];
    nPvec=[25];
    beta_vec=[1];
    for dim=1:length(nAvec)
        nP=nPvec(dim);
        nA=nAvec(dim);
        N=nP+nA;
        c=4./(N.^(0.8));
        d=1;
        for bb=1:length(beta_vec)
            beta=beta_vec(bb);
           
            sigma_vec=[0.005];
            
            for s=1:length(sigma_vec)
                sigma=sigma_vec(s);
                
                xeq=ones(N,1);%0.5+(1.5-0.5).*rand(N,1);
                b=zeros(N,N);
                
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
                b_antA=b(1:nA,1:nA);
                b_antP=b(nP+1:N,nP+1:N);
                
                Aold=zeros(N,N);
                for i=1:N
                    Aold(:,i)=-xeq.*bold(:,i);
                end
                
                alpha=b*xeq;
                T=100*N;
                for t=1:T
                    bnewA=b(1:nA,nA+1:N);
                    bnewP=b(nA+1:N,1:nA);
                    ind0=bnewA==0;
                    ind1=bnewA~=0;
                    allind=1:length(bnewA(:));
                    die=ceil(rand*length(bnewA(ind1)));
                    new=ceil(rand*length(bnewA(ind0)));
                    [i1,j1]=ind2sub(size(bnewA),allind(ind1));
                    [i0, j0]=ind2sub(size(bnewA),allind(ind0));
                    
                    bnewA(i0(new),j0(new))= bnewA(i1(die),j1(die));
                    bnewP(j0(new),i0(new))= bnewP(j1(die),i1(die));
                    bnewA(i1(die),j1(die))=0;
                    bnewP(j1(die),i1(die))=0;
                    dim2=size(bnewA);
                    
                    bnew= [b_antA, bnewA; bnewP,b_antP];
                    xeqE=bnew\alpha;
                    
                    if sum(xeqE>0)-length(xeqE)==0 && sum(xeqE)>= sum(xeq)
                        xeq=xeqE;
                        b=bnew;
                    end

                end
                A=zeros(N,N);
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
