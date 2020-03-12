function [kd,kl] = new_gaussiansimilarity(A)
%gaussianssimilarity
%% ��������
interaction = A;
[nd,nl]=size(interaction);
yitad = 1;
yital =1;
lanbuda = 0.5;
gamadd = 0.0000001;
gamall = 0.00000001;

%calculate gamad for Gaussian kernel calculation
 for i=1:nd
   sd(i)=norm(interaction(i,:))^2;
 end
   gamad=nd/sum(sd')*gamadd;

    %calculate gamal for Gaussian kernel calculation
for i=1:nl
        sl(i)=norm(interaction(:,i))^2;
    end
    gamal=nl/sum(sl')*gamall;
    %calculate Gaussian kernel for the similarity between disease: kd


for i=1:nd
    for j=1:nd
        kd(i,j)=exp(-gamad*(norm(interaction(i,:)-interaction(j,:)))^2);
    end
end
    %calculate Gaussian kernel for the similarity between lnc-RNA: kl

for i=1:nl
    for j=1:nl
        kl(i,j)=exp(-gamal*(norm(interaction(:,i)-interaction(:,j)))^2);
    end
end

%��kd�Mһ�����������ռ�����������֮������ƶȾ��󣬶����е�C=-15��D=log��9999��

for i=1:nd
    for j=1:nd
        kd(i,j)=1/(1+exp(-15*kd(i,j)+log(9999)));
    end
end




