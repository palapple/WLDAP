function [final_Rscore]=WLDAP(A)

%lncRNA_disease_interaction: adjacency matrix for the lncRNA_disease associations
%lncRNA_disease_interaction(i,j) means lncRNA i associated with disease j

% nlA:the number of lncRNAs in A
% ndA:the number of diseases in A

[nlA,ndA] = size(A);

% Initialize matrix Rsocre1_lnc_dis
Rscore1_lnc_dis=zeros(nlA,ndA);

%calculate the corresponding weight matrix W in A
for i=1:nlA
        q=bsxfun(@rdivide,repmat(A(i,:),nlA,1).*A,sum(A));
        W_A(i,1:nlA)=1./sum(A,2).*sum(q,2);
end
%calculate the level of consistency between the contribution of resource moved in both directions
W=W_A';
W=W./(repmat(sum(W),nlA,1));
W_A = (W_A+W);
%obtain the first level of resource score about lncRNA_disease_associations
Rscore1_lnc_dis= W_A*A;

% ndB:the number of diseases in A
% nlB:the number of lncRNAs in A

[nlB,ndB]= size(A);

final_Rscore=zeros(nlB,ndB);

%construct lncRNA_disease_gene tripartite graph
%calculate the final resource socre Rscore to infer potenial lncRNA-disease associations
final_Rscore= (Rscore1_lnc_dis);

end


