

## What is a biomarker?

Biomarker is general term that refers to a measure (a piece of data) of a putative or known biological process or state.

In the context of biomedicine, biomarkers are used in the clinical context for monitoring and diagnosis. Serum cholesterol (HDL, LDL) is an example of a cardiac biomarker.

According to the [Wikipedia page on biomarkers](https://en.wikipedia.org/wiki/Biomarker), 
    
    In 1998, the National Institutes of Health Biomarkers Definitions Working Group defined a biomarker as "a characteristic that is objectively measured and evaluated as an indicator of normal biological processes, pathogenic processes, or pharmacologic responses to a therapeutic intervention

## In what context is this analysis being conducted?

In this context, we are concerned with a scenario where a patient has already been diagnosed, and we are assessing the impact of treatment on outcomes.

## What types of cancer biomarkers are you concerned with?

'''Biomarker''' is a very general term. Biomarkers can be characterized by how they are measured (e.g. molecular, genetic, epigenetic, etc) and by their utility.

For this analysis, we characterize biomarkers as:

1. *prognostic*: correlated with outcome irrespective of treatment
2. *predictive*: correlates with benefit from treatment (ie outcome of treatment varies according to value of the biomarker). Is typically measured pre-treatment.
3. *surrogate*: early indicator that a patient is or is not benefiting from treatment. Is typically measured post-treatment.

There are other uses of biomarkers, such as diagnosis for example. In this case we are assuming diagnosis was made appropriately, so we are not considering those here.

***Note:***

We are primarily interested in the '''predictive''' biomarkers, but often these are indistinguishable from prognostic. For this reason, we include both.

Typically surrogate markers are considered an added bonus. We include surrogate here because they are useful & relevant to the method being proposed.

## How are predictive & prognostic biomarkers different?

  (insert graphic of predictive biomarker)


  (insert graphic of predictive biomarker)

## How are biomarkers typically assessed?

To be a valid biomarker, a measurement has to be (first and foremost) accurate, reliable, and related to the underlying process it purports to measure.

Having established this as a baseline, biomarkers are then assessed for their prognostic/predictive utility in the context of a clinical study, registry or meta-analysis. 

Typical measures for predictive biomarkers include

### Improvement in AUC 

This method uses the AUC or "area under the ROC curve" as a measure of the predictive value of the biomarker. AUC-comparison methods are one method used for model comparison -- while they are popular there is some disagreement over whether AUC is the best measure of a model's fit. 

But, the main question the method tries to answer is '''How much better can we predict the outcome given the biomarker, than we could without it'''.

Typical workflow for this type of analysis proceeds as follows:
    * Patients are scored as "benefit" or "not benefit" according to criteria determined by the investigators. Often, '''benefit''' = disease-free survival, or absence of disease progression.
    * AUC0: A baseline model is developed, given the best biomarker, or the set of validated biomarkers. The AUC (which reflects sensitivity/specificity) is estimated for the baseline case.
    * AUC1: A new model is developed, given the biomarker being considered. The AUC for this biomarker is established.
    * The difference between AUC0 & AUC1 is estimated. This is often done by bootstrap, the Venkatraman (parametric) comparison, or a non-parametric method such as that proposed by deLong.

It's worth noting that this approach is very general.

Its utility does not depend on the specifics of the model used (whether survival vs logistic or something else), and of the process for model comparison (AUC curve analysis vs something else).

This method also makes no distinction between predictive vs prognostic indicators. 

Very often, the analysis only includes patients who (a) received treatment, and who (b) would have otherwise had a poor prognosis.  

The assumption is then that any improvement in outcomes that is correlated with biomarker values can be attributed to treatment benefit.

### FPR / TPR 

Another common method used is to summarize the sensitivity / specificity of the biomarker. This is practically very useful since it provides information in a method that can be immediately applied to the clinical context. 

This method borrows from the diagnostic-biomarker use case, however instead of classifying the "true" state of underlying disease the biomarker is used to distinguish between patients who "benefit" from treatment (vs those who didn't).

### Regression analysis 

Another common approach to biomarker evaluation is to perform a regression analysis on survival (or cause-specific survival).

Regression coefficients are evaluated to see if the biomarker is independently and/or "significantly" correlated with the outcome.

This has the benefit of allowing for distinction between predictive vs prognostic biomarkers. A predictive biomarker, for example, would be expected to have a significant interaction with treatment.

### Potential outcomes framework 

An underlying problem with almost all of these models is that of '''causal inference'''. They all suffer from the problem inherent to most analyses - to what extent are observed correlations causal? To what extent are they generally applicable?

The classic approach to this problem (and from which all methods implicitly borrow) is the [potential outcomes framework]().  Some methods of biomarker analysis reference this framework explicitly, but all borrow from it in some way.

## What is a generative model?

Generative model is a model whose likelihood function is exposed, and that contains an explicit "data-generating process" from which one could simulate data.  A typical "GLM" model is a generative 

## Why is a generative model better than what we're doing? (embedded in this question: what are we doing, i.e. how are we currently validating putative biomarkers?)

First, 
## Why is a generative model better than a feature importance analysis in linear models/decision trees?

1. encodes/includes domain-specific knowledge & the state of cancer oncology
2. straightforward way to include data from external studies (either as informative priors or in a meta-analysis)
3. knowledge can accumulate. As our knowledge improves, the model will improve
4. greater variety 

## How can generative model answer open questions about putative biomarkers



Generative models 


