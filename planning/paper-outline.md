
## Background

* current approaches to biomarker detection
    - potential outcomes framework
    - AUC 
    - TPR ; FPR

* tradeoff between 
    1. interpretability
        * Are results summarized in a manner that decision-makers can understand? 
        * What is the likely outcome of a particular scenario?
            - IE given patient characteristics and a planned course of treatment?
    2. illuminating underlying mechanisms
        * for patients who did benefit, do we understand why they benefitted?
        * why are patients with this biomarker value likely to benefit from treatment?

* goal of this method is to identify which patients are likely to respond to treatment
    - Difficult because of
        1. Each study has a small sample size 
        2. Expect that many biomarkers contribute to outcomes / benefit, given complexity of disease
        3. Often have redundant biomarkers (ie several ways of capturing same type of information)
        4. Many biomarkers are measured imperfectly - e.g. ITH which is inferred from 
    - Analysis is confounded due to
        1. diversity of disease states among patients within a study
        2. diversity of study designs and treatment regimens across studies
        3. diversity of disease types/sites across studies
        3. Lack of uniformity of data across studies (same biomarkers not captured for every study)
    - Ideally, we would
        1. include data from all studies that speak to the process relevant to a particular biomarker
        2. include prior information given what we know about biological processes underlying disease
        3. include all measures of a biomarker to better estimate reliability of the measure

* context for biomarker analysis in area of cancer immunotherapy
    - difficult to ascribe causal relationship to biomarker*therapy interaction
    - cancer immunogram proposing broad framework for decision making informed by biomarkers
    - systems biology approaches often too complex
    - desire for knowledge to accumulate

## The proposed method 

* Intro proposed method 
    - mechanistic model, borrowing features from ecological modeling & from pharmacodynamics (Pk/Pd models)
        - compartment model 
        - latent terms for major components of immune-tumor interaction
        - subject-level frailty estimates 
    - provides a natural way to incorporate new potential biomarkers into the analysis

* Evaluation methodology
    - consider two different simulation processes 
        - process matching the proposed model
        - process according to complex systems-biological interactions
    - evaluate using two different analysis routines
        - standard survival analysis for mortality / illness 
        - generative model with limited covariate set 
        - generative model with expanded covariate set

* Summarizing results 
    - net impact of biomarker
    - identify subsets of patients likely to respond

## results - part 1

* Data simulation - matching proposed model

* Evaluation of straightforward survival model

* Evaluation of generative model

* Comparing generative model vs survival model

## results - part 2

* Data simulation - dynamic systems biology model

* Evaluation of straightforward survival model

* Evaluation of generative model

* Comparing generative model vs survival model

## summary & conclusions


