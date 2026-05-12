clc; 
clear all
close all

SNR_dB = 10;  
SNR_linear=10.^(SNR_dB/10.);
sigma2=1/SNR_linear;
sample =500; 
sparsity=12;

%%% system parameters
N = 256*2; % number of beams (transmit antennas)
L =4; % number of all paths
gamma=0.5; 
Lf =ceil(L*gamma); % number of paths for far-field 
Ln = L-Lf; % number of paths for near-field
pilot = 128:32:32*8; % number of pilot overhead
len = length(pilot);

fc = 30e9; % carrier frequency
c = 3e8;
lambda_c = c/fc; % wavelength 
d = lambda_c / 2; % antenna space

% the far-field angle-domain DFT matrix
Uf =(1/sqrt(N))* exp(-1i*pi*[0:N-1]'*[-(N-1)/2:1:(N/2)]*(2/N));

% the near-field polar-domain transform matrix [5]
Rmin=10;
Rmax=50;
eta = 2.5; 
[Un, label, dict_cell, label_cell] = QuaCode(N, d, lambda_c, eta, Rmin, Rmax);
S=size(Un,2);
% error_omp_dft=zeros(sample,len);
% error_omp_qua=zeros(sample,len);
% error_homp=zeros(sample,len);
% error_LS=zeros(sample,len);
% error_MMSE=zeros(sample,len);
% energy=zeros(sample,1);
% 
Rh=zeros(N,N);


for s=1:sample
    s
    [h,hf,hn] = generate_hybrid_field_channel(N, Lf, Ln, d, fc,Rmin, Rmax);
    
    for iS=1:len
        M=pilot(iS)
        P=zeros(M,N);
        Noise=zeros(M,1);
                P = Generate_pilot(M,N);

%         P=sqrt(.5)*(1i*((rand(M,N)>0.5)*2-1)/sqrt(M)+((rand(M,N)>0.5)*2-1)/sqrt(M)); % pilot matrix
%        for m=1:N
%            P(:,m)=P(:,m)/norm(P(:,m));
%        end
%         P=P/norm(P);
               % P=sqrt(.5)*(1i*rand(M,N)+rand(M,N)); % pilot matrix

        noise = sqrt(sigma2)*(randn(M,1)+1i*randn(M,1))/sqrt(2);
        y=P*h+noise;
        
        %% the far-field OMP based scheme with DFT matrix
        hshat_omp_dft = OMP(y,P*Uf,sparsity*(Lf+Ln));
        hhat_omp_dft = Uf*hshat_omp_dft;
        error_omp_dft(s,iS)=sum(abs(hhat_omp_dft-h).^2)/norm(h,'fro')^2;
        
        %% the near-field OMP based scheme with polar-domain transform matrix
        [hshat_omp_qua,pos_xhat] = OMP(y,P*Un,sparsity*(Lf+Ln));
        hhat_omp_qua = Un*hshat_omp_qua;
        error_omp_qua(s,iS)=sum(abs(hhat_omp_qua-h).^2)/norm(h,'fro')^2;
        
        %% the proposed hybrid-field OMP based scheme
        hhat_homp=Hybrid_OMP(y,P,Uf,Un,sparsity*Lf,sparsity*Ln);
        error_homp(s,iS)=sum(abs(hhat_homp-h).^2)/norm(h,'fro')^2;
        
       %% the LS
       hhat_LS=h+sqrt(sigma2)*(randn(N,1)+1i*randn(N,1))/sqrt(2);
       error_LS(s,iS)=sum(abs(hhat_LS-h).^2)/norm(h,'fro')^2;
       
       %% My algorithm
               
       [SBL_xhat,err(s,iS)] =log_SBL_NF(P*Uf,y,N,sigma2,M,h,Uf);
       error_sbl(s,iS) =  norm(h-Uf* SBL_xhat,'fro')^2/norm(h,'fro')^2;
       %% Practical algorithm
       psi=[P*Uf P*Un];
       [Hhat_homp_1,hhat_f] = Hybrid_OMP_1( y , P*Uf,P*Un,      12,  L);
       H_prac=Uf* Hhat_homp_1(1:N)+Un* Hhat_homp_1(N+1:N+S);
        error_practical(s,iS) =  norm(h-H_prac ,'fro')^2/norm(h,'fro')^2;
        
     %% Far field Burst
     
     [hburst]=Burst_CE(y,P,N,N,1);
     error_burst(s,iS) =  norm(h-hburst ,'fro')^2/norm(h,'fro')^2;

    end
end
 nb_of_best_samples=ceil(.95*sample)


nmse_omp_dft = mean(mink(error_omp_dft,nb_of_best_samples))
nmse_omp_qua = mean(mink(error_omp_qua, nb_of_best_samples))
nmse_homp = mean(mink(error_homp,nb_of_best_samples))
nmse_burst=mean(mink(error_burst,nb_of_best_samples))
nmse_practical=mean(mink(error_practical,nb_of_best_samples))
nmse_proposed=mean(mink(err,nb_of_best_samples))
nmse_proposed2=mean(mink( error_sbl,nb_of_best_samples))

nmse_omp_dft=10*log10(nmse_omp_dft)
nmse_omp_qua=10*log10(nmse_omp_qua)
nmse_homp=10*log10(nmse_homp)
%nmse_LS=10*log10(nmse_LS)
%nmse_Hamza=10*log10(mean(err,1))
nmse_proposed=10*log10(nmse_proposed2)
nmse_practical=10*log10(nmse_practical)
nmse_burst=10*log10(nmse_burst)
X_axe=pilot;
figure('color',[1,1,1]); 
ha=gca;
plot(X_axe,nmse_proposed,'-.s',X_axe,nmse_burst,'--^',X_axe,nmse_practical,'--',X_axe,nmse_omp_dft,'<-',X_axe,nmse_omp_qua,'>-',...
        X_axe,nmse_homp,'rs-','linewidth',2.5);

legend('Proposed','Far-field-Unifom Burst [6]','hybrid-field OMP -unknwon $L$ [14]','OMP-angular domain',...
'OMP- Polar domain [3]','hybrid-field OMP [13]','Interpreter','latex','FontSize',10)
xlabel('\bf Pilot length $T$','Interpreter','latex','FontSize',14)
xlim([X_axe(1) X_axe(end)])
xticks(pilot)
ylabel('\bf NMSE [dB]','Interpreter','latex','FontSize',14)    
grid on

% semilogy(1+abs(Uf*h),'-','linewidth',1.5)
% 
% xlabel('\bf Pilot length $T$','Interpreter','latex','FontSize',14)
% xticks(pilot)
% ylabel('\bf NMSE [dB]','Interpreter','latex','FontSize',14)    
% grid on



% hold on
% plot(SNR_dB,nmse_omp_qua,'b>-','linewidth',1.5);
% hold on
% plot(SNR_dB,nmse_homp,'rs-','linewidth',1.5);
% hold on
% plot(SNR_dB,nmse_Hamza,'k--','linewidth',1.5);
% hold on
% grid on
% legend('Far-field OMP [3]','Near-field OMP [6]','Proposed hybrid-field OMP','MMSE')
% xlabel('SNR (dB)')
% ylabel('NMSE (dB)')
% hold off