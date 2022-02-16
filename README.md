# MLVar

This repository contains the code to reproduce results in:

Nicora, G., Zucca, S., Limongelli, I. et al. **A machine learning approach based on ACMG/AMP guidelines for genomic variant classification and prioritization**. Sci Rep 12, 2517 (2022). https://doi.org/10.1038/s41598-022-06547-3

### AIM
Classification and prioritization of genomic variants associated with inherited disorders. 

### Background: 
inherited variants are nowadays interpreted according to the ACMG/AMP guidelines [1], which define 5 different tiers of classification: *Pathogenic, Likely pathogenic, Benign, Likely benign, VUS (variant of uncertain significance)*. Therefore, when evidence is not enough for classification, a variant can be still uncertain (VUS). To interpret VUS variants, data-driven approaches, such as Machine Learning algorithms, can be used. 

### Methods
Variants are annotated and interpreted according to the ACMG/AMP guidelines [1] by the eVai software [2]. ACMG/AMP levels of evidence are used as features for Machine learning approaches (in particular, Logistic Regression). 

### Results
We applied our framework on different datasets. Results are shown in the Notebook.


#### References

[1] Richards, S., Aziz, N., Bale, S. et al. Standards and guidelines for the interpretation of sequence variants: a joint consensus recommendation of the American College of Medical Genetics and Genomics and the Association for Molecular Pathology. Genet Med 17, 405–423 (2015). https://doi.org/10.1038/gim.2015.30

[2] https://www.engenome.com/wp-content/uploads/2020/10/eVai_enGenome_WhitePaper_v0.7.pdf
