function [B] = get_interaction(A)
A=textread('knowndiseaselncrnainteraction.txt');
% nd:the number of diseases
% nl:the number of lnc-RNAs
% pp:the number of known diseae-lncRNA associations
nd=max(A(:,1)); 
nl=max(A(:,2));
[pp,qq]=size(A);
%interaction: adjacency matrix for the disease-lncRNA association network
%interaction(i,j)=1 means lnc-RNA j is related to disease i
for i=1:pp
        interaction(A(i,1),A(i,2))=1;
end
B = interaction;
