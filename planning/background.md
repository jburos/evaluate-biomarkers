
## Types of biomarkers 

* risk/prognostic indicators (correlated with outcome)
* predictive indicators (correlate with benefit from treatment)
* interim/proxy indicators (early signals that a patient is benefiting from treatment)

We are interested in identifying and reporting on biomarkers of all three types.

Note that these are not mutually exclusive, since a biomarker (like PD-L1 expression) can for example can be *both* a prognostic & a predictive indicator.

## Primary goals of biomarker analysis 

(in an ideal sense).. 

1. Basic science: To know why a treatment regimen failed/succeeded
2. Translational medicine: To inform treatment decisions in practice

In practice, both are difficult. We will never know the first & can only approximate it by using inferences drawn from a generative model. Even the second goal is difficult for several reasons, mostly related to small sample sizes, which I will summarize below.

Even though these goals serve different users, and each imposes its own set of requirements, they are not inconsistent. 

## How these two goals are different

Let's think a bit more about how these two goals are different from one another.

#### The basic science use case 

For the basic science use case, we care more about **parameter values** than we do about **prediction of outcomes**. We are happy to to **include all data available**, including redundant signals where two biomarkers are essentially measuring the same process, insofar as this supports our goal. Results can be **hypothesis-generating**, and so nuanced descriptions of patterns are tolerated. 

The critical constraint on this analysis is that we want to know _why_ a particular outcome was observed, and so a ''black-box'' model is less ideal.  A result that averages over many models, or whose data-generating-process is opaque will not help promote this goal.

#### The translational medicine use case

For the translational medicine use case, in some ways, we are more about the **likely outcomes** under different treatment scenarios for a particular patient.  By comparison with the basic science use case, we care less about _why_ the outcomes are worse under treatment A than B, so long as we can reliably identify the patients for whom this will be true & summarize the magnitude of the difference in outcome under both treatment regimens. 

The most critical constraint on this use case is that the summary be **relevant** to the decision-making process. A secondary goal is that it should be simple to apply. This secondary goal implies that, all else being equal, we would prefer a result that requires fewer measurements. And, the method should enable summaries of outcomes where less data is available.

## What the two goals have in common

The argument put forth here is that the two goals -- even though they require different summaries of data -- are not that different from one another. 

They are both looking at the same data, and thus they suffer from the same sets of problems.

Specifically:

1. Most clinical studies including biomarker data have very small sample sizes (~10-30)
2. Many only include treated patients 
3. Many include either a variety of treatments, or a single treatment
     - in any case, each patient only receives one treatment
4. Many include patients across a diversity of disease states and/or types

There is a philosophical argument re: why one should use a generative model for translational medicine, but we will not consider that here. We will instead focus on the practical issue of how to overcome these difficulties in the method of analysis. 

## The two goals in detail

To see how similar these two goals are in practice, let's consider an analysis of a particular biomarker (PD-L1 expression) in the context of immunotherapy.

The basic science goal of this analysis may be to see how PD-L1 expression can help us understand why some patients did or did not respond to treatment. We may gain some useful insight from this analysis, even if PD-L1 expression is not directly correlated with overall survival. This may be one of many features used to characterize the tumor microenvironment for each patient, which may be used to determine the importance of other biomarkers under consideration.

The translational medicine use case, by comparison, has a goal of summarizing each patient's likely outcome under competing therapies. The fundamental problem is that most patients only receive a single treatment. Most analyses resolve this problem by using a potential outcomes framework -- that is, they use a model to estimate the expected outcome under an alternative therapy based on how similar patients fared under alternative treatments. 

This method is only as good as the model which is used to estimate alternate treatments, and to identify "similar" patients. 

There is also a need to determine how representative the study population is, since the outcomes of these summaries should (ideally) be generalizable to a target treatment population. 

Both methods suffer from the same problems of having a small, diverse & possibly non-representative patient population. In order to achieve the goal of deriving value from the data which these patients have put themselves at risk to generate, it would behoove us to include as many external sources of information as possible in our analysis. 

## Defining the generative model

Broadly speaking, generative models are preferred when the goal of an analysis is to draw inferences about parameters of interest, in addition to generating predictions.

A generative model is any model where the data-generating process can be described, and where this is used to define the likelihood function for the model. Technically, most linear regressions are generative models since they expose a likelihood function and one can use that likelihood to simulate data according to the model. 

In the context of cancer, I am referring to a specific subset of generative models -- namely, those that replicate the data-generating process we believe operates in the context of cancer therapies.

## Benefits of a generative model

1. Include or encode what we know about cancer progression & interaction with treatments
    - specifically, include known biological processes relevant to the treatment 
2. Gain strength from including redundant biomarkers -- IE biomarkers that measure the same process
3. Account for the diversity of patients and/or disease states in a study
4. Provide a reasonable way to include data from other studies
    - can use other-study data to inform priors on parameters 
    - or, can include similar studies in a form of a meta-analysis
5. Can benefit from including data for redundant biomarkers
6. Can incorporate intermediate or proxy outcomes

Taken in concert, the model can address many of the the concerns 

## Best practices in biomarker analysis 

1. Define a core set of biomarkers which should be evaluated on most if not all studies
2. Define the target population for the translational medicine use case
3. Irrespective of how the model is defined, summarize the result in terms that are relevant for decision-making
4. Include data from external studies where possible
5. Use the most sensitive outcome measure available 
    - time to mortality/censor if available
    - semi-competing-risks analysis including progression and/or mortality
    - joint growth model of tumor burden over time if available
6. Clinical data should include any inclusion/exclusion criteria used to define the study population
    - e.g. don't filter on age without collecting age as a covariate




