clc;clear;

A = load('D:\MATLAB\Bidirectional_label_propagation\WLDAP\WLDAP\know_dis_lncRNA_3.txt');

[nl,nd] = size(A);
A_ori = A;
% [kl,kd] = gaussiansimilarity(A_ori,nl,nd);
[kd,kl] = new_gaussiansimilarity(A_ori);
 lncRNAsimilarity = kl;
 diseasesimilarity = kd;


 [score_ori] = WLDAP(A_ori);
 score_ori_ori = score_ori;
index = find(1 == A_ori);
%%  5-fold
auc = zeros(1,100);
for k = 1:100
    k
    indices = crossvalind('Kfold', length(index), 5);
    A = A_ori;
    score_ori=score_ori_ori;
for cv = 1:5
       cv;
       index_2 = find(cv == indices);
       A(index(index_2)) = 0;
        % [kl,kd] = gaussiansimilarity(A,nl,nd);
        [kd,kl] = new_gaussiansimilarity(A);
        lncRNAsimilarity = kl;
        diseasesimilarity = kd;

    for u = 1:length(index_2)
    [lncRNA,disease]=ind2sub(size(A),index(index_2(u)));
      if sum(A(lncRNA,:))==0
          A(lncRNA,:)=Sim_lnc(A,lncRNAsimilarity,lncRNA);
          B = max(A(lncRNA,:))；
          C = min(A(lncRNA,:))；
          for j=1:nd
            A(lncRNA,j)=(A(lncRNA,j)-C)/(B-C)；
        end
      end
      if sum(A(:,disease))==0
          A(:,disease)=Sim_dis(A,diseasesimilarity,disease);
          B = max(A(:,disease));
          C = min(A(:,disease));
          for j=1:nl
            A(j,disease) = (A(j,disease)-C)/(B-C);
        end
      end
    end
% ��һ��
for i = 1:nl
    B = max(A(i,:));
    C = min(A(i,:));
    for j =1:nd
        A(i,j) = (A(i,j) - C) / (B-C);
    end
end
     %%%����÷־���
%%  ���ݴ���


%% ������
      [result] = WLDAP(A);
      score_ori(index(index_2)) = result(index(index_2));
       A = A_ori;
end
%% ��auc����
    pre_label_score = score_ori(:);
    label_y = A_ori(:);
    auc(k) = roc_1(pre_label_score,label_y,'red');

end
%%
auc_ave = mean(auc);
auc_std = std(auc);
% x(n) = gama;
% y(n) = auc_ave;
% n = n+1;
% end
