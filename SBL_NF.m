function [SBL_xhat,err] = SBL_NF(A,y,N,sig2e,Tau_p,Hsf,dft_h)
max_iter_sbl=100;M=size(y,2);
rho_0=.01;
lambda=1/M;
stopping_creterion=1e-2;
Gamma(:,:,1)=.1*eye(N);
SBL_SRR=zeros(1);
sbl_max_iter=1;

u=zeros(N,N);
c=zeros(N,N);
gamma=1./real(diag(Gamma(:,:,1)));%B=eye(N,N);
w=ones(N,N);l_k=1;
for t=1:max_iter_sbl
F1=Gamma(:,:,t)*A'/(A*Gamma(:,:,t)*A'+sig2e*eye(Tau_p));
    %Sigma_y=(A*Gamma(:,:,t)*A'+sig2e*eye(Tau_p));
    mu_x=F1*y;
    Sigma_x=Gamma(:,:,t)-F1*A*Gamma(:,:,t);
    for m=1:M
        St(:,:,m)=mu_x(:,m)*mu_x(:,m)';
    end
    S = real(sum(St,3)/M + Sigma_x);
    if t<sbl_max_iter
        gamma=real(diag(S));
%         for n=1:N
%             gamma(n)=real(S(n,n));
%         end
    elseif t>=sbl_max_iter
        
        rho=max(1e-4,.5*rho_0*(1+cos(t*pi/max_iter_sbl)));
        l_k=3;max(1,round(sqrt(M)-1));
        for i=1:N
           % kappa=100;
            for j=1:i-1
                % w(i,j)=exp(-.05*abs(log(1+kappa*gamma(i))/log(1+kappa)-log(1+kappa*gamma(j))/log(1+kappa))^2);
                mu=min(.01*t*M,2);
            w(i,j)=exp(-1*abs(log(gamma(i))-log(gamma(j)))^2);
           w(i,j)=exp(-1*abs(gamma(i)-gamma(j))^2);

            %  w(i,j)=exp(-1*abs(norm(mu_x(i))-norm(mu_x(j))));
              w(j,i)=w(i,j);
            end
        end
    
        for tep=1:2
            
            for i=1:N
                for j=1:i-1
                    gmbar=real(w(i,j)*((gamma(i))-(gamma(j)))-u(i,j)/rho);
                    c(i,j)=real(sign(gmbar)*max(0,abs(gmbar)-lambda/rho)); %%B(i,j)*
                    c(j,i)=-c(i,j);
                    u(i,j)=real(u(i,j)+rho*(c(i,j)-w(i,j)*((gamma(i))-(gamma(j)))));
                    u(i,j)=real(u(i,j)+rho*(c(i,j)-w(i,j)*(gamma(i)-gamma(j))));
                    u(j,i)=-u(i,j);
                end
            end
            bdiag=gamma;
            for i=1:N
                if i==1
                    se=[];
                else
                    se=max(1,i-l_k):i-1;
                end
                if i==N
                    se2=[];
                else
                    se2=i+1:min(i+l_k,N);
                end
                gm_bar=sum((c(i,se)).'+(w(i,se)).'.*(bdiag(se))+u(i,se).'/rho);
                gm_tilde=sum((c(se2,i))-(w(se2,i)).*(bdiag(se2))+u(se2,i)/rho);              % Nhat=length(se)+length(se2);
                Nhat=sum(w(i,[se,se2]),2);
                %gamma(i)=abs(real(S(i,i)/(1+rho*Nhat-rho*(gm_bar-gm_tilde))));
                tem=S(i,i)+rho*(gm_bar-gm_tilde);
                delta=tem^2+4*rho*Nhat;
                gamma(i)=(-tem+sqrt(delta))/(2*rho*Nhat);
 if gamma(i)>1e5
    gamma(i)=1e5;
 end
%                 gamma(i)=1e-9;
%             end
                %for i=1:N
%             if gamma(i)<1e-7
%                 gamma(i)=1e-9;
%             end
           % bdiag(i)=gamma(i);
            end
            
            GM(:,tep)=gamma;
            if tep>2&& (norm(GM(:,tep)-GM(:,tep-1))^2/norm(GM(:,tep))^2)<.01
                break;
            end
        end
    end
  
    Gamma(:,:,t+1)=(diag(1./gamma));
    xhat_sbl(:,:,t)=mu_x;
%    Err(t)=norm(mu_x - x,'fro')^2/norm(x,'fro')^2;
    Err(t)=norm(Hsf-dft_h* mu_x,'fro')^2/norm(Hsf,'fro')^2;
    Ey=(A*Gamma(:,:,t)*A'+sig2e*eye(Tau_p));
    ww=0;
    for i=2:N
       for j=i-1
        ww=ww+lambda*w(i,j)*abs(log(gamma(i)-log(gamma(j))));
    end
    end
WW(t)=ww;
    obj(t)=real(log(det(Ey))+trace(y*y'*inv(Ey)));
    %obj(t)=sum(diag(real(S))./gamma)+sum(log(gamma))+ww;
    %Err(t)=norm(x-dftmtx(N)* mu_x)^2/norm(x)^2;
    if mod(t,20)==0
        Err(t);
    end
    res_norm(t)=norm(y-A*mu_x);
    
    if t>30
       % norm(xhat_sbl(:,:,t)-xhat_sbl(:,:,t-1))^2;
        if norm(xhat_sbl(:,:,t)-xhat_sbl(:,:,t-1),'fro')<stopping_creterion || isnan(norm(gamma))
            break;
        end
    end
end
[~,I]=min(obj);%(res_norm);
%[~,I]=min(obj(2:end));

SBL_xhat=xhat_sbl(:,:,min(I+1,length(Err)));
err=min(Err);
% act=find(gamma>.0001);
% ch_est_ls = zeros(N,1);
%ch_est_ls(act) = pinv(A(:,act))*y;
%ErrLS(t)=norm(Hsf-dftmtx(N)* ch_est_ls,'fro')^2/norm(Hsf,'fro')^2;
%% 
%err=Err(end);
%% 
end



