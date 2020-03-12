
clc,clear;
A = load('D:\MATLAB\Bidirectional_label_propagation\WLDAP\WLDAP\known_lncRNA_disease_interaction_1.txt');
[nl,nd] = size(A);
% [kl,kd] = gaussiansimilarity(A,nl,nd);
[kd,kl] = new_gaussiansimilarity(A);
 lncRNAsimilarity = kl;
 diseasesimilarity = kd;
A_ori = A;
index = find(1 == A_ori);
[score_ori] = WLDAP(A_ori);
 A = A_ori;
for i = 1:length(index)
    i
    A(index(i)) = 0;
    % [kl,kd] = gaussiansimilarity(A,nl,nd);
    [kd,kl] = new_gaussiansimilarity(A);
    lncRNAsimilarity = kl;
    diseasesimilarity = kd;
    [lncRNA,disease]=ind2sub(size(A),index(i));
      if sum(A(lncRNA,:))==0
          A(lncRNA,:)=Sim_lnc(A,lncRNAsimilarity,lncRNA);
          B = max(A(lncRNA,:));
          C = min(A(lncRNA,:));
          for j= 1:nd
            A(lncRNA,j)=(A(lncRNA,j)-C)/(B-C);
          end
      end
      if sum(A(:,disease))==0
          A(:,disease)=Sim_dis(A,diseasesimilarity,disease);
          B = max(A(:,disease));
          C = min(A(:, disease));
          for j= 1:nl
            A(disease,j)= (A(disease,j)-C)/(B-C)
          end
      end
    [result]=WLDAP(A);
    score_ori(index(i)) = result(index(i));
    A = A_ori;
end
pre_label_score = score_ori(:);
save score_ori;
label_y = A_ori(:);
auc=roc(pre_label_score(:),label_y,'b')
