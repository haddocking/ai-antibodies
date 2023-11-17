# AI-antibodies

This is the repo with all the material for generating and analysing HADDOCK3 files for machine learning-modelled antibodies.

# 1. create the environment with conda
```
conda create -n aiabs python=3.10.6
```

# 2. install the package with pip
```
pip install .
```
# 3. run the setup for your benchmark

for each pdb in your benchmark set, you can run the following command:
```
aiabs ${PDB} --input_dir=benchmark_haddock_23_May_2023/${PDB} --output_dir=${path_to_output_dir} --act_act_path=./bin/generate-act-act.sh
```

# 4. analyse the results
```
aiabs-analysis --ref_folders=haddock_runs_caprieval
```

You can exclude some antibodies from the analysis
```
aiabs-analysis --ref_folders=haddock_runs_caprieval --exclude_pdbs=7seg,7kn4
```

for plotting the results you can access the script available in the [analysis](https://github.com/mgiulini/aiabs/tree/main/analysis) folder.

# Notebooks
Additional plots for correlating DockQ vs confidence measures and against interface/paratope/epitope RMSD can be created
using the notebooks in `analysis/notebooks`.