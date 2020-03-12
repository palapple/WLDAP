clc,clear;
%% ��������
lncRNA_Disease=load('known_lncRNA_disease_interaction.txt');

pre_label_score2 = importdata('pre_label_score_RWR_5fold.mat');
pre_label_score4 = importdata('pre_label_score_SIMCLDA_5fold.mat');
pre_label_score1 = importdata('pre_label_score_WLDAP_5fold.mat');
pre_label_score3 = importdata('pre_label_score_LRLSLDA_5fold.mat');
pre_label_score5 = importdata('pre_label_score_RWRLDA_5fold.mat');

label_y = lncRNA_Disease(:);
roc_1(pre_label_score1(:),label_y,'red');
roc_1(pre_label_score2(:),label_y,'green');
roc_1(pre_label_score3(:),label_y,'blue');
roc_1(pre_label_score4(:),label_y,'m');
roc_1(pre_label_score5(:),label_y,'c');
% roc_1(pre_label_score1(:),label_y,'green');

title('5fold');

