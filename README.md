# Covid19_systems-level_analysis
Systems-level analysis reveals repurposable drugs and compounds against SARS-CoV-2

Fangyuan Chen†, 1, 2, Qingya Shi†, 1, 2, Fen Pei†,1,3, Andreas Vogt†,1,3, Rebecca A. Porritt4,8 Gustavo Garcia Jr5, Angela C. Gomez4, Mary Hongying Cheng1, Mark Schurdak1,3, Bing Liu1, Stephen Y. Chan6,7, Vaithilingaraja Arumugaswami5, Andrew M Stern1,3, D Lansing Taylor1,3, Moshe Arditi4,8, and Ivet Bahar*, 1,3

1. Department of Computational and Systems Biology, School of Medicine, University of Pittsburgh, Pittsburgh, PA, 15213, USA; 2. School of Medicine, Tsinghua University, Beijing, 100084, China; 3. University of Pittsburgh Drug Discovery Institute, Pittsburgh, PA 15213; 4. Department of Pediatrics, Division of Pediatric Infectious Diseases and Immunology, Cedars-Sinai Medical Center, Los Angeles, CA, 90048; 5. Department of Molecular and Medical Pharmacology, David Geffen School of Medicine, University of California and Eli and Edythe Broad Center of Regenerative Medicine and Stem Cell Research, University of California, Los Angeles; 6. Pittsburgh Heart, Lung, Blood, and Vascular Medicine Institute, and 7. Division of Cardiology, Department of Medicine, University of Pittsburgh Medical Center, Pittsburgh, PA, 15217; 8. Biomedical Sciences, Infectious and Immunologic Diseases Research Center, Cedars-Sinai Medical Center, Los Angeles, CA, 90048

† These authors made equal contribution
* Correspondence: Dr. Ivet Bahar

==================================================================

Codes

DESeq2_antiviral.ipynb: code for evaluation of host-targeted antiviral signature from SARS-CoV-2-infected A549 cells.

DESeq2_anticytokine.ipynb: code for evaluation of host-targeted anti-hyperinflammatory signature from SARS-CoV-2-infected A549-ACE2 cells.

CTI_similarity.py: code for calculating interaction-pattern-based similarities between compounds

Clustering.py: code for clustering the compounds into different groups based on the CTI similarities calculated above. 

==================================================================

Websites

QuartataWeb: http://quartata.csb.pitt.edu/, used to retrieve the known drug-target interactions from DrugBank and STITCH datasets

Enrichr: https://maayanlab.cloud/Enrichr/, used for pathway enrichment on 348 host genes associated with SARS-CoV-2 infection
