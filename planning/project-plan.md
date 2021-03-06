Overview
--------

Recently, several clinical trials have been proposed to evaluate
potential for checkpoint blockade to improve treatment of cancer and
eradicate disease. Concurrently, researchers have also identified
several [biomarkers]() of therapeutic efficacy. While there is a
pressing need for such research, there is also a potential for
confusion.

Specifically, because of the following factors:

1.  Not all biomarkers are measured in every cohort
2.  Variety of cancers being evaluated
3.  Variety of therapies under investigation
4.  **Small sample sizes** - always make inferences difficult

The goal of this project is to consider how best to evaluate potential
biomarkers in this context.

Standard approach
-----------------

Given a clinical trial with study data on a set of patients, there is a
standard approach to the analysis of a potential biomarker for clinical
efficacy.

The typical approach to this kind of question is to :

1.  Develop a model for clinical outcome given all known data (e.g.
    demographics, stage of disease, etc).
2.  Given the model in (1), estimate portion of clinical benefit due to
    the therapy.
3.  Given the model in (2), estimate whether clinical benefit varies
    according to measured value of biomarker.
4.  Revise model in (2), to include features (e.g. other biomarkers)
    known to correlate with benefit from disease.
    -   After accounting for known biomarkers, known correlates of
        clinical outcome, and the therapeutic intervention, how much of
        the variance in outcome is explained by the biomarker?

Improved standard approach
--------------------------

There are some ways we can improve upon the standard approach.

For example:

1.  Use a more informative prior on expected outcome, by incorporating
    survival rates typical according to published literature
2.  Include populations where biomarkers were used to screen patients
    for enrollment
    -   e.g. some studies screen for patients with particular molecular
        attributes believed to be more likely to benefit from therapy

Various roles biomarkers play
-----------------------------

Putative biomarkers can be separated into several groups, depending on
their purpose / correlation. This also speaks to the type of decision
the biomarker would be used to inform.

Biomarkers may be:

1.  Correlated with outcome (survival), but do not modify treatment
    effect on outcome
2.  Used to select patients likely to benefit from treatment
3.  Used as an indicator of treatment efficacy - ie measured during or
    following therapy but strongly indicates likelihood of positive
    outcome

Best case scenario
------------------

Sometimes it can be helpful to consider what a best case scenario would
look like.

Ideally, this research would

1.  Evaluate each new biomarker in the context of *all* known and
    putative biomarkers
2.  Adjust their particular analysis for the ways their sample cohort
    varies from those studied to date
3.  Update or summarize the state of knowledge, to identify how well the
    set of all biomarkers is performing

A base model for disease
------------------------

For this to work, we need to arrive at a plausible base model for
disease progression & how this is modified by the therapeutic agent.

This model won't be perfect - it's meant to be a work in progress, but
the idea is to build on this over time, to incorporate impact of
potential biomarkers on the therapeutic agent.

### The tumor microenvironment

One idea is to model disease (tumor progression) as an outcome of a
dynamic system.

We will start with a simplified model of the [tumor microenvironment]().

There are a few key players of interest:

1.  Tumor cells (the tumor)
    -   have a number (N)
    -   have an inherent rate of growth/decay (change in N)

2.  Immune cells (the immune system)
    -   have a number (N)
    -   have an inherent rate of growth/decay (change in N)

3.  Cytokines
    -   number may not be that important
    -   are secreted by both Tumor cells and Immune cells
    -   often cancel each other out (ie one cytokine has effect of
        inhibiting another)
    -   have net effect of modifying the parameters of interest

There are other attributes, but for now let's just consider these.

We will address each in turn.

### Immune activation

Thinking simplistically, the primary method of therapy is to *enable the
immune system*. Modeling the level of immune-system activation is
therefore pretty important to our goal.

Let's think about the primary way that these two players interact.

1.  **immune activation** Some portion of immune cells are activated
    against (primed to kill) tumor cells
2.  **tumor resistance** Some portion of tumor cells resistent to
    immune-cell-mediated cell death

From a very simplistic level, we can combine these two values to yield a
**net activation** of the immune cells.

IE:

net\_activation = immune\_activation\*(1 - tumor\_resistance)

So, if we have 80% of immune cells activated to kill the tumor, but 60%
of tumor cells are resistant to immune-cell-mediated cell death, then
the net activation proportion is 0.8\*(1-0.6) or 0.32.

Tumor cells obviously interact with the immune system via cytokines, and
we'll get to that in a bit. For now, let's try to keep things simple.

### Tumor size

The primary effect of immune system activation that we care about is to
reduce the size of the tumor.

Now, there are some limits (possibly) for the number of cells in the
tumor at any time, but for now let's start with a simple exponential
growth model.

IE: at any time t, the tumor size (N) can be modeled as a function of
the size in the previous time (t-1) & the amount of tumor cell death:

tumor\_N\[t\] = tumor\_N\[t-1\] \* tumor\_growth\_rate\[t\] -
tumor\_cell\_death\[t\]

(note that if we have tumor heterogeneity, we may have a diversity of
growth rates. Here we're thinking of an average growth rate over all
cells in the tumor).

### Tumor cell death

The next step is to think about how to determine the amount of
**tumor\_cell\_death** at any point in time.

This will be a function of the number of immune cells and their net
level of activation (all at time t).

tumor\_cell\_death\[t\] = immune\_N\[t\] \* net\_activation\[t\]

### Immune "size" (N cells)

Finally, we will need to think about a model for how many immune cells
have infiltrated the tumor microenvironment.

Let's start by assuming that the immune cells mediate their own "growth"
via cytokines. We can call this rate of growth immune\_rate

immune\_N\[t\] = immune\_N\[t\] \* immune\_growth\_rate\[t\]

This **immune\_growth\_rate** is, in turn, heavily influenced by the
cytokine environment -- we'll hold off on informing this for now. We
have enough here to start thinking about how to put the model together.

### Patient-level hazard model

Often the primary outcome of interest is cancer-free days, or event-free
survival. The final component of this model will be to include an
underlying **hazard** indicating the patient's risk for an adverse event
at a given time.

It's worth noting that there are likely other factors that contribute to
that underlying hazard -- some of which may be biomarkers! -- but for
now we want to start with the minimal number of components.

hazard\[t\] ~ tumor\_size\[t\] + noise

Simulating data
---------------

Before trying to fit a model to real data, let's simulate some data
according to the hypothesized model to see how the components work
together.

First pass at a model
---------------------

Known biomarker: CD8+ T cells
-----------------------------

One example biomarker is the proportion of T cells which are CD8+. This
is a useful biomarker to start with -- it is known to correlate with
positive outcomes and has a proposed mechanism of action.

Patients with higher proportions of CD8+ T cells should have higher
**immune\_activation**.

Let's walk through how this could be incorporated into the model for
disease proposed here.

Potential biomarker: mutational burden
--------------------------------------

Let's illustrate the potential of this approach by considering a
potential biomarker & attempting to do analysis of its impact on
event-free survival following treatment.

### What are the hypothesized mechanisms of action?

Lingering questions
-------------------

What about a biomarker whose impact isn't mediated by one of the
mechanisms shown here, but whose correlation with treatment benefit is
nonetheless robust?
