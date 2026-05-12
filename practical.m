function [H_prac] = practical(y,psi_f,psi_n,psi,Uf,Un,L)
gamma=1;
Omega_f=[];
r=y;
Nt=size(psi_f,2);
[N_t,S]=size(psi_n);  

%% COarse Estimation
for l=1:L
   R(:,l)=r; hf_hat=zeros(Nt,1);
   [~,n_star]=max(vecnorm(psi_f'*r,2,2));
   Omega_f=[Omega_f n_star];
   hf_hat(Omega_f,1)=pinv(psi_f(:,Omega_f))*y;
   r=y-psi_f*hf_hat;
end
%% COarse Estimation

r_m=r;
hm_hat=hf_hat;
 h_hat=zeros(Nt+S,1);
for gammaL =(L-1):-1:0
    
    Omega=Omega_f(1:gammaL);
    r=R(:,gammaL+1);
    
    for l=1:(L-gammaL)
        [~,n_star]=max(vecnorm(psi_n'*r,2,2));
        Omega=[Omega n_star];
        h_hat(Omega)=pinv(psi(:,Omega))*y;
        r=y-psi*h_hat;
    end
    [gammaL norm(r)]
    if norm(r)<norm(r_m)
        r_m=r;
        hm_hat=h_hat;
        %break;
    end
end

        H_prac=Uf* h_hat(1:Nt)+Un* h_hat(Nt+1:Nt+S);
end

