function [P] = Generate_pilot(M,N)
P=sqrt(.5)*(1i*((rand(M,N)>0.5)*2-1)/sqrt(M)+((rand(M,N)>0.5)*2-1)/sqrt(M)); % pilot matrix
        for i=1:N %normalize!
            P(:,i)= P(:,i)/norm( P(:,i));
        end
end

