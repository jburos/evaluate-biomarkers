
## References reviewed

#### [Key Concepts and Limitations of Statistical Methods for Evaluating Biomarkers of Kidney Disease](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC4116071/) 

Not cancer-specific but contains a useful overview of statistical methods used to evaluate potential biomarkers. Here the application is to evaluate biomarkers for diagnostic utility, rather than their ability to benefit from treatment. 

1. reviews traditional and novel measures, including a nice [table](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC4116071/table/t1/) summarizing advantages & disadvantages of each.
2. highlights benefit of some approaches which have easy clinical interpretatation. IE result is clinically intuitive. This is one disadvantage of AUC curve analysis, for example.
3. yields a nice summary of the "ideal" statistical method for biomarker evaluation: 
        1. Result has a clear interpretation that is clinically relevant
        2. Evaluates the relative performance of models, not just individual variables for treatment selection
        3. Can handle continuous & categorical variables. 
        4. Allows comparison of biomarkers that are differently distributed
        5. Can evaluate both the contribution of biomarker to baseline risk as well as differential response to treatment.

4. key quotations regarding AUC curve improvement: 
       * "The area under the ROC curve (AUC) is probably the most widely used summary index. "
       * "Furthermore, the AUC of a new biomarker is highly dependent on its comparison with the gold standard. In the presence of an imperfect gold standard, such as serum creatinine for the cases of AKI and CKD, the classification potential of the new biomarker may be falsely diminished."
       * "The widely used method by DeLong et al.18 is designed to nonparametrically compare two correlated ROC curves (clinical model with and without the biomarker); however, it has recently been shown that the test may be overly conservative and may occasionally produce incorrect estimates. Begg et al have used simulations to show that the use of same risk predictors from nested models while comparing AUCs with and without risk factors leads to grossly invalid inferences."
       * "Researchers have observed that ΔAUC depends on the performance of the underlying clinical model. For example, good clinical models are harder to improve on, even with markers that have shown strong association."

5. alternative approaches include Integrated Discrimination Improvement (IDI) & Net Reclassification Index (NRI).
       * there is a categorical & continuous version of NRI. Categorical is recommended where there is an existing model to classify risk categorically that is well defined. Otherwise, continuous is recommended. Continuous NRI is considered "inflated", however, and may not be reliable.
       * The IDI is recommended for continuous outcomes, and can be applied to risk assessment for both event & non-event scenarios. This is essentially the sum of improvement of sensitivity+(1-specificity) for event & non-event patients. 
           - There is not a well-established cutoff of IDI that indicates a useful biomarker
           - There is a graphical summary proposed which displays the components of the IDI.

#### [Biomarker-Based Prediction of Response to Therapy for Colorectal Cancer](http://ajcp.oxfordjournals.org/content/ajcpath/134/3/478.full.pdf) - 2010
       * uses HR for each biomarker, as potential independent predictor of survival as selection criteria
       * segments biomarkers into those that predict adverse events (toxicity) vs benefit from therapy (efficacy). Biomarkers are summarized by their prognostic utility for each of several common therapeutic approaches. Highlights the complexity of the application & accumulation of knowledge in this area.

#### [Biomarker-based selection of therapy for colorectal cancer](http://www.sci-hub.bz/10.2217/bmm.11.38) - 2011
 One year after the previous paper investigating potential biomarkers, many of those summarized above are now mentioned as being part of standard practice. Many of these are very specific to disease, germline mutations which predict toxicity to therapy, and specific mechanisms of biological action which may be impaired. 
 
Article contains a nice "executive summary" listing the many biomarkers in current use: 

##### The following methods are used for the prediction of response to primary treatment for colorectal cancer

* Pathologic and molecular staging
* Oncotype Dx® colon cancer test
      - Impact of adjuvant chemotherapy in higher risk primary CRC not as validated as for breast cancer.
* ColoPrint®
      - Requires fresh tissue, is microarray based and has limited prospective validation.
* Circulating tumor cells
      - Established independent prognostic factor for circulating tumor cells (CTCs). Enhanced methods are needed to improve capture efficiency. There are limited biomarker testing data on CTCs.
* Microsatellite instability
      - Associated both with predisposition for CRC and prognosis in sporadic disease, is a validated marker (low microsatellite instability associated with adverse outcome) and has achieved wide use in routine practice.

##### The following biomarkers are used for the prediction of efficacy and toxicity of therapeutic antibodies for colorectal carcinoma:

* EGFR expression, amplification and mutation
      - Limited use given lack of efficacy in CRC for anti-EGFR tyrosine kinase inhibitors.
* KRAS mutation
      - KRAS mutation is associated with resistance to anti-EGFR antibodies. KRAS mutation testing is now required by the US FDA for selecting patients for cetuximab/panitumumab. Traditional sequencing methods focus on hot-spot codons 12 and 13 in exon 2, but recent evidence may support including codon 61.
* BRAF mutation
      - BRAF mutations seen in approximately 10% of CRCs that are wild-type for KRAS. Early studies favored the hypothesis that V600E BRAF mutation was associated with cetuximab/panitumumab resistance. Recent evidence suggests that BRAF mutation is a significant adverse prognostic factor and not a specific predictive factor of drug resistance.
* PIK3CA mutation
* p53 mutation
      - Immunohistochemistry test not accurate. Approximately 1–5% of KRAS wild-type colorectal cancers will have a p53 gene mutation. KRAS wild-type p53 mutated CRCs have been associated with a positive response to cetuximab-based therapy.
* PTEN expression
      - PTEN expression loss is identified in approximately 40% of CRCs. PTEN expression loss in KRAS wild-type tumors has been associated with cetuximab resistance.
* IGF-1 receptor pathway
* EGFR ligands
      - In KRAS wild-type CRC, low expression of the EGFR activating ligands amphiregulin and epiregulin by tumor cells has been associated with resistance to anti-EGFR antibody therapy and adverse clinical outcome.
* Fcg receptor status and antibody-dependent cellular cytotoxicity (ADCC) response
      - Preliminary studies have linked both germline polymorphisms and post-translational modifications (glycosylation and fucosylation) of the Fcg receptor with impaired ADCC response to cetuximab and panitumumab.

##### The following biomarkers are used for the prediction of efficacy and toxicity of cytotoxic chemotherapy for colorectal carcinoma:

* Thymidylate synthase (TS)
      - Direct inhibition of TS is considered to be the major mechanism of 5-FU’s antitumor mechanism. Early reports indicated that increased TS expression was associated with resistance to 5-FU treatment. Numerous more recent studies using both slide-based immunohistochemistry and mRNA expression detected by RT-PCR have failed to achieve a consensus.
* Dehydropyrimidine dehydrogenase (DPD)
      - Reports of DPD expression as a predictor of 5-FU efficacy are generally less convincing than for TS. In certain patients with germline polymorphisms of the DPD gene, significant reduction in enzymatic activity is present resulting in severe toxicity when 5-FU is administered.
* Methylenetetrahydrofolate reductase
      - Large studies have not confirmed the prognostic or predictive value of methylenetetrahydrofolate reductase (MTHFR) genotyping. MTHFR polymorphisms have also been associated with capecitabine toxicity.
* Topoisomerase I
      - Topoisomerase I expression has been correlated with irinotecan response in several large studies.
* UGT1A1
      - Polymorphisms in the UGT1A1 gene have been strongly linked to toxicity of irinotecan-based treatment for CRC. One particular polymorphism, homozygosity for the seven-repeat allele, is associated with severe diarrhea when irinotecan is administered.
* Excision repair cross complementing (ERC)C1
      - ERCC1 expression may predict for resistance to FOLFOX treatment of CRC.
* Ribonucleotide reductase M1
      - Ribonucleotide reductase M1 (RRM1) expression has been linked to response of lung cancer to DNA damaging agents. 

       
#### [Mechanism-driven biomarkers to guide immune checkpoint blockade in cancer therapy](http://www.nature.com/nrc/journal/v16/n5/full/nrc.2016.36.html)  - 2016


       


    