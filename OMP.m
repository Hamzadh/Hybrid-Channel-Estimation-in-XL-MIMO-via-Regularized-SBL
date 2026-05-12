function[hhat,support]=OMP(y,A,t)
%y:measurement
%A:sensing matrix
%t:number of iteration

[M,N]=size(A);
r=y;  %Initialization of residual
support=[];

for i=1:t
	Product = A'*r;
    [~,pos]=max(abs(Product));
	support=[support,pos(end)];
	hhat=zeros(N,1);
	hhat(support)=pinv(A(:,support))*y;
    r=y-A*hhat;  
end

