
# predictive response to chemo biomarkers in UTUC

* 2009 - []()

* 2014 - [Comprehensive molecular characterization of urothelial bladder carcinoma](http://www.nature.com/nature/journal/v507/n7492/full/nature12965.html)

* 2011 - [Biomarkers for prognosis and treatment selection in advanced bladder cancer patients](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3257805/) Review paper on known prognostic & predictive biomarkers
	- describes known prognostic biomarkers in Table 1
	- describes known predictive biomarkers in [Table 2](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3257805/table/T2/)
		- ERCC1 & MDR1 (low expression -> likely benefit from cisplatin-based adjuvant chemotherapy)
		- ERCC1 (low expression -> benefit from chemoradiotherapy)
		- BRCA1 (low expression -> benefit from neoadjuvant cisplatin chemotherapy)
		- MRE11 (high expression -> benefit from radical radiotherapy than from cystectomy)
		- ET1 & ETaR (inhibitors of these compounds should be administered pre-metastasis rather than post-metastasis)
	- most biomarkers not yet validated in external cohorts 
		- should be validated since validation has proven useful in invalidating previously-identified biomarkers
			> A study testing published prognostic mRNA based gene signatures for survival found that none could be validated on external data, either in combined NMI and MI tumors or MI tumors alone [50*]. Furthermore, the published gene signatures for survival performed no better than randomly selected genes in external datasets, suggesting that genes which correlate with survival are difficult to identify and may vary between different patient populations and treatment regimens [50*]. Importantly, validating the prognostic biomarkers in multiple independent patient cohorts is essential for ensuring their clinical utility.
	- Utilities to help translate biomarkers to clinical application could be useful
		- discusses COXEN
			> An algorithm called CO-eXpression ExtrapolatioN (COXEN) was developed to predict response to chemotherapeutic agents in bladder, breast, and ovarian cancer patients based on drug response in cell lines [51, 52], and for in silico drug discovery [51, 53]. 
		- not sure if this is still being used. Interesting that it is unrelated to the above-described research.

* 2011 - [Upper urinary tract urothelial carcinoma: what have we learned in the last 4 years?](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3150069/)
	- literature review including results from the UTUCC (Upper Tract Urothelial Carcinoma Consortium) 

* 2012 - [ERCC1 as a biomarker for bladder cancer patients likely to benefit from adjuvant chemotherapy](http://bmccancer.biomedcentral.com/articles/10.1186/1471-2407-12-187)
	- retrospective analysis of 93 patients following cystectomy
		- 57 treated with cisplatin-based adjuvent chemotherapy
		- 36 were not treated
	- evaluated ECRR1 expression by IHC staining
	- ERCC1 as prognostic biomarker
		- overall, little prognostic value of ERCC1
	- ERCC1 is predictive biomarker
		- higher survival with ERCC1 expression in patients without chemotherapy
		- lower survival with ERCC1 expression in patients with chemotherapy
		- significant interaction term (p ~ 0.03) for differential association of ERCC1 expression with survival, depending on treatment

* 2016 - [BcCluster: A Bladder Cancer Database at the Molecular Level](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC4927921/)
	- consolidation of molecular data on bladder cancer from over 100 studies
	- over 1k protein-protein interactions

* 2016 - [Prognostic factors and predictive tools for upper tract urothelial carcinoma: a systematic review.](http://www.ncbi.nlm.nih.gov/pubmed/27101100)
	- due to low incidence rate, heterogenous population & small sample sizes of existing studies, it's been difficult to establish validated biomarkers for UTUC.
		> This review is the result of the International Consultation on Upper Tract Urothelial Carcinoma (ICUD-UTUC) endeavor, providing an overview of existing prognostic biomarkers and factors, and predictive tools.
	- prognostic biomarkers
		- smoking status
			- cessation time of >10y at time of diagnosis mitigates prognostic risk of smoking
		- age & gender are controversial prognostic indicators
			- age at diagnosis not independent correlate with risk, but is epidemiologically related to more aggressive tumor growth
			- gender shows similar divergent association in epidemological/clinical studies in bladder cancer
			- gender does not appear to be prognostic indicator for UTUC
		- obese patients have higher incidence of more advanced stage UTUC, and worse prognosis
		- tumor stage / grade
		- lymph node metastasis
		- LVI (lymphovascular invasion)
		- tumor architecture
		- tumor multifocality
		- tumor necrosis
		- high NLR
		- high fibrinogen
		- 
	- Conclusion re prognostic biomarkers:
		> This systematic review of the literature highlights the relatively low level of evidence of most of the studies which have been conducted on prognostic factors in UTUC. When strictly applying the OCEBM criteria [9], these retrospective case series should be considered as level 4 of evidence. The absence of level 1 of evidence studies constitutes a major limitation and reduces the thrust of any conclusion that can be established from these data. 
		> But the usage of highly selective inclusion criteria for this systematic review permitted us to isolate a few well-conducted studies from high-volume centers and international collaborations that could therefore be upgraded to level 3 of evidence.
	- Conclusion re predictive biomarkers:
		> Another key message is that predictive tools that have been developed so far seem to be accurate and well calibrated. 
		> However, only Yates’ model has been externally validated in a population different from the original one. Such approaches should be further performed before either of the existing tools could be put into widespread use.



## Methodological considerations (types of biomarkers; biomarker study design)

* 2012 - [Cancer biomarkers: selecting the right drug for the right patient](http://www.nature.com.sci-hub.bz/nrd/journal/v11/n3/full/nrd3651.html)
	- overview of state of biomarkers, including similar definitions to those defined above
	- section on clinical trial designs, appropriate for various types of predictive biomakers (imaging, genetic signature, etc)
	- TODO read more thoroughly

* 2015 - [Biomarker Development in the Context of Urologic Cancers](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC4521394/)
	- defines 4 types of biomarkers
		- **prognostic biomarker**
			- used to make decisions regarding clinical management of patient
			- related to risk or outcome but does not correlate with differential response to treatments
		- **predictive biomarker**
			- used to determine sensitivity to treatment options
			- typically assessed prior to treatment
			- is often / may also be a prognostic biomarker
			- gives targeted therapies as examples of scenarios where predictive biomarker evaluation is part of drug-development process
				- define **integral biomarker** as one that informs choice of therapy on an investigational protocol
				- often accompanied by FDA-approved companion diagnostic, often as strictly regulated as the drug itself
		- **response indicator** 
			- occur during & after treatment
			- show that there has been a change following treatment, but do not necessarily indicate benefit
			- e.g. pharmacodynamic attributes of drug metabolism, tumor regression, mRECIST scores, etc
			- e.g. in bladder cancer (for BCG treatment):
				- urinary IL-2, IL-8, IL-18, and TNF levels
		- **efficacy response biomarker**
			- also occur during & after treatment
			- indicate that the patient is benefiting from treatment. 
			- are robust enought to be considered surrogate endpoints in context of a clinical trial
			- e.g. clinical events related to disease progression
	- requirements and aspects of biomarker validation
		- analytical validation: 
			- assess performance characteristics of the biomarker measurement itself
			- consistency & reproducibility of the measurement
			- if seeking FDA validation, the assay used to measure the biomarker subject to same criteria as IND
		- clinical validity:
			- demonstration that the biomarker is fit for particular context of clinical use
			- gold standard here would be a randomized trial using the biomarker in clinical decision-making context vs not 
		- clinical utility:
			- demonstration that a decision made in light of biomarker improves a patient's risk/reward ratio for a clinical outcome
			- can relate to mitigating risk, minimizing cost, improving diagnostic accuracy, etc. 
		- clinical qualification: 
			- An additional step related to regulatory implications of the biomarker
			- e.g. if FDA approval is sought, or if the biomarker is critical for FDA approval of a related drug 

* 2016 - [Biomarker-Guided Adaptive Trial Designs in Phase II and Phase III: A Methodological Review.](http://www.ncbi.nlm.nih.gov/pubmed/26910238)
	- more recent review paper focusing on design of studies to validate biomarkers

