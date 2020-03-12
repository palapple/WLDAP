# WLDAP
## License
Copyright (C) 2019 LiFeng Wu(2431705149@mail2.gdut.edu.cn),Guobo Xie(xiegb@gdut.edu.cn)

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program; if not, see http://www.gnu.org/licenses/.

**Type:** Package Title: WLDAP: a computational model of weighted lncRNA-disease associations prediction

**Description:** This package implements the WLDAP algorithm, predicting lncRNA-disease associations. 	



## Requirement
8GB memory
MATLAB R2017a or later



## Related data informations need to first load 

**known_lncRNA_disease_interaction_1.txt:** The known lncRNA-disease associations is represented by adjacency matrix A in our method,which shows binary associations between diseases and lncRNAs.1 represents disease j is associated with lncRNA i,otherwise 0. This dataset includes 115 lncRNAs and 178 diseases, with 540 kinds of experimentally verified associations.

Then there's the comparison dataset, which is  **lncRNA_Disease_Matrix_2.txt**, **know_dis_lncRNA_3.txt**

The Related data informations are placed in the data folder.

In order to load related data informations, you should input the appropriate code in the matlab Command Window:

```
A = load('.\known_lncRNA_disease_interaction_1.txt');
```

##  Calculate the cosine similarity matrix and integrate 
new_gaussiansimilarity.m is used to calculate the similarity of diseases and lncRNA;
Sim_dis.m and Sim_lnc.m  are used to integrate disease similarity and lncRNA similarity, respectively.
you should input the appropriate code in the matlab Command Window:
```
 [kd,kl] = new_gaussiansimilarity(A);                 % Return the processed gaussiansimilarity matrix
inter_prof_vec=Sim_lnc(A,lncRNAsimilarity,lncRNA)  %calculate isolated lncRNA interaction profile vector

inter_prof_vec=Sim_dis(A,diseasesimilarity,disease) %calculate isolated disease interaction profile vector
```

## Run WLDAP to infer potential associations between lncRNAs and diseases

To analyze these data on WLDAP to further infer potential associations between lncRNAs and diseases, you should input the appropriate code in the matlab Command Window:

```
%% Calculated scoring matrix
[final_Rscore]=WLDAP(A)
%% result
allresult(disease,lncRNA,interaction,score);
```

**WLDAP.m:** WLDAP core algorithm to predict potential lncRNA-disease associations; it corresponds to the part of network  consistency projection in this paper.

**allresult.m:** the predicted results will be automatically saved in the excel table(allresult.xlsx)

## Cross validation

**WLDAP_LOOCV.m:** function Leave one for cross-validation. 

**pre_label_score_WLDAP_LOOCV.mat**the prediction score of Leave one for cross-validation.

**WLDAP_5fold.m:** function 5-fold cross-validation.

**pre_label_score_WLDAP_5fold.mat**the prediction score of 5-fold cross-validation.

Running  WLDAP_LOOCV.m get the AUC value and ROC chart after implementing Leave one for cross-validation.


Running WLDAP_5fold.m get the AUC value and ROC chart after implementing 5-fold cross-validation.

## Self comparison
We performed  Leave one for cross-validation in self comparison.

**pre_label_score_NCPHLDA.mat:** the prediction score in the disease space and lncRNA spase.


**duibi.m:** the code for self-comparison.
