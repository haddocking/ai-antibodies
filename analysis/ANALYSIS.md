# Analyse results

In this folder there are a few scripts that are useful to analyse the results of our HADDOCK3 runs.

Before starting, be sure to have the aiabs python environment active.

```
conda activate aiabs
```

## 1 generate_ref_data

first of all, we want to generate our reference data, that will be accessed by the majority of the scripts present here

```
python generate_ref_data.py
```

## 2 bound-af2 comparison graphs

To create a comparison graph, such as the one in Fig.1 of the paper, please run:

```
python protocols_comparison.py
```

## 3 clustering

To generate the same plots, but focusing on the clustering stage

```
python clustering.py
```

## 4 epitope-paratope-comparison

To generate plots comparing epitope-paratope accuracy

```
python epitope_paratope_comparison.py
```

## 5 zdock comparison

To generate a plot illustrating the comparison between the performances of HADDOCK3 and those of ZDOCK, such as that in Figure XXX of the paper:

```
python zdock_comparison.py
```

## 6 violin plot

To generate the violin plots comparing CDR H3 RMSD and Paratope RMSD:

```
python violinplot.py
```

## 7 energy minimisation assessment

To assess the impact of energy minimisation on rigidbody models:

```
python rigid_emref_assessment.py
```

## 8 refinement 

To assess the impact of the refinement on the different protocols
```
python analyse_refinement_step.py
```
