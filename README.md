# KnowEnG's Gene Signature Pipeline 
This is the Knowledge Engine for Genomics (KnowEnG), an NIH BD2K Center of Excellence, Signature Analysis Pipeline.

This pipeline performs network-based **signature analysis** on the columns of a given spreadsheet, where spreadsheet's columns correspond to sample-labels and rows correspond to gene-labels.  The signature is based on correlating gene expression data (network smoothed) against known gene signature data.

There are four clustering methods that one can choose from:


| **Options**                                      | **Method**                           | **Parameters** |
| ------------------------------------------------ | -------------------------------------| -------------- |
| Signature                                        | cosine                               | cos            |
| Consensus Signature                              | bootstrapping with cos               | cc_cos         |
| Signature  with network regularization           | network-based cos                    | net_cos        |
| Consensus Signature  with network regularization | bootstrapping with network-based cos | cc_net_cos     |


Note: all of the signture methods mentioned above use the cosine similarity (cos) as the main similarity algorithm.


* * * 
## How to run this pipeline with Our data
* * * 
### 1. Clone the Gene_Signature_Pipeline Repo
```
 git clone https://github.com/KnowEnG-Research/Gene_Signature_Pipeline.git
```
 
### 2. Install the following (Ubuntu or Linux)
  ```
 pip3 install pyyaml
 pip3 install knpackage
 pip3 install scipy==0.18.0
 pip3 install numpy==1.11.1
 pip3 install pandas==0.18.1
 pip3 install matplotlib==1.4.2
 pip3 install scikit-learn==0.17.1
 
 apt-get install -y python3-pip
 apt-get install -y libfreetype6-dev libxft-dev
 apt-get install -y libblas-dev liblapack-dev libatlas-base-dev gfortran
```

### 3. Change directory to Gene_Signature_Pipeline

```
cd Gene_Signature_Pipeline
```

### 4. Change directory to test

```
cd test
```
 
### 5. Create a local directory "run_dir" and place all the run files in it
```
make env_setup
```

### 6. Use one of the following "make" commands to select and run a clustering option:


| **Command**         | **Option**                                              | 
|:------------------- |:------------------------------------------------------- | 
| make run_cos        | cosine similarity                                       |
| make run_net_cos    | cosine similarity with network regularization           |
| make run_cc_cos     | Consensus cosine similarity                             |
| make run_cc_net_cos | Consensus cosine similarity with network regularization |

 
* * * 
## How to run this pipeline with Your data
* * * 

__***Follow steps 1-3 above then do the following:***__

### * Create your run directory

 ```
 mkdir run_directory
 ```

### * Change directory to the run_directory

 ```
 cd run_directory
 ```

### * Create your results directory

 ```
 mkdir results_directory
 ```
 
### * Create run_paramters file  (YAML Format)
 ``` 
 Look for examples of run_parameters in the Gene_Signature_Pipeline/data/run_files zTEMPLATE_cc_net_cos.yml
 ```
### * Modify run_paramters file  (YAML Format)
Change processing_method to one of: serial, parallel depending on your machine.
```
processing_method: serial
```

set the data file targets to the files you want to run, and the parameters as appropriate for your data.


### * Run the Gene Signature Pipeline:

  * Update PYTHONPATH enviroment variable
   ``` 
   export PYTHONPATH='../src':$PYTHONPATH    
   ```
   
  * Run
   ```
  python3 ../src/samples_signature.py -run_directory ./run_dir -run_file zTEMPLATE_cc_net_cos.yml
   ```

* * * 
## Description of "run_parameters" file
* * * 

| **Key**                   | **Value** | **Comments** |
| ------------------------- | --------- | ------------ |
| method                    | **cos**, **cc_cos**, **net_cos** or **cc_net_cos**  | Choose clustering method |
| gg_network_name_full_path | directory+gg_network_name |Path and file name of the 4 col network file |
| spreadsheet_name_full_path | directory+spreadsheet_name|  Path and file name of user supplied gene sets |
| phenotype_data_full_path | directory+phenotype_data_name| Path and file name of user supplied phenotype data |
| threshold | 10 | cluster eval - catagorical vs continuous cut off level |
| results_directory | directory | Directory to save the output files |
| tmp_directory | directory | Directory to save the intermediate files |
| rwr_max_iterations | 100| Maximum number of iterations without convergence in random walk with restart |
| rwr_convergence_tolerence | 1.0e-8 | Frobenius norm tolerence of spreadsheet vector in random walk|
| rwr_restart_probability | 0.7 | alpha in `V_(n+1) = alpha * N * Vn + (1-alpha) * Vo` |
| rows_sampling_fraction| 0.8| Select 80% of spreadsheet rows|
| cols_sampling_fraction| 0.8| Select 80% of spreadsheet columns|
| number_of_bootstraps| 4 | Number of random samplings |
| number_of_clusters| 3 | Estimated number of clusters |
| cos_conv_check_freq| 50 | Check convergence at given frequency |
| cos_max_invariance| 200 | Maximum number of invariance |
| cos_max_iterations| 10000 | Maximum number of iterations |
| cos_penalty_parameter| 1400 | Penalty parameter |
| top_number_of_genes| 100 | Number of top genes selected |
| processing_method| serial or parallel or distribute | Choose processing method |
| parallelism| number of cores to use in parallel processing | Set number of cores for speed or memory |

gg_network_name = STRING_experimental_gene_gene.edge</br>
spreadsheet_name = ProGENI_rwr20_STExp_GDSC_500.rname.gxc.tsv</br>
phenotype_data_name = UCEC_phenotype.txt

* * * 
## Description of Output files saved in results directory
* * * 

* Output files of all four methods save genes by sample heatmap variances per row with name **genes_variance_{method}_{timestamp}_viz.tsv**.</br>

 |  |**variance**|
 | :--------------------: |:--------------------:|
 | **gene 1**|float|
 |...|...|
 | **gene m**| float|

* Output files of all four methods save genes by samples heatmap with name **genes_by_samples_heatmp_{method}_{timestamp}_viz.tsv**.</br>

 |  |**sample 1**|...|**sample n**|
 | :--------------------: |:--------------------:|:--------------------:|:--------------------:|
 | **gene 1**|float|...|float|
 |...|...|...|...|
 | **gene m**|float|...|float|

* Output files of all four methods save samples by samples heatmap with name **consensus_matrix_{method}_{timestamp}_viz.tsv**.</br>

 |  |**sample 1**|...|**sample n**|
 | :--------------------: |:--------------------:|:--------------------:|:--------------------:|
 | **sample 1**|float|...|float|
 |...|...|...|...|
 | **sample n**|float|...|float|
 
* Output files of all four methods save patients to cluster map with name **samples_labeled_by_cluster_{method}_{timestamp}_viz.tsv**.</br>

 |    |**cluster**|
 | :--------------------: |:--------------------:|
 | **sample 1** | int|
 |...|...|
 | **sample n** |int|
 
* Output files of all four methods save gene scores by cluster with name **genes_averages_by_cluster_{method}_{timestamp}_viz.tsv**.</br>

 |  |**cluster 1**|...|**cluster k**|
 | :--------------------: |:--------------------:|:--------------------:|:--------------------:|
 | **gene 1**|float|...|float|
 |...|...|...|...|
 | **gene m**|float|...|float|
 
* Output files of all four methods save spreadsheet with top ranked genes per sample with name **top_genes_by_cluster_{method}_{timestamp}_download.tsv**.</br>

 |  |**cluster 1**|...|**cluster k**|
 | :--------------------: |:--------------------:|:--------------------:|:--------------------:|
 | **gene 1**|1/0|...|1/0|
 |...|...|...|...|
 | **gene m**|1/0|...|1/0|
  
* All four methods save **silhouette number of clusters** and **corresponding silhouette score** with name silhouette_average\_{method}\_{timestamp}\_viz.tsv.</br>
 ```
 File Example: 
 silhouette number of clusters = 3, corresponding silhouette score = 1
 ```

* Output files of all four methods save patients to cluster map with name **phenotypes_labeled_by_cluster_{method}_{timestamp}_viz.tsv**.</br>

 | **sample id** |**cluster**|**phenotype 1**|...|**phenotype k**|
 | :--------------------: |:--------------------:|:--------------------:|:--------------------:|:--------------------:|
 | **sample 1**|int|mixed type|...|mixed type|
 |...|...|...|...|...|
 | **sample n**|int|mixed type|...|mixed type|
 
 
 * The clustering evaluation output file has the name 
 **clustering_evaluation_result_{timestamp}.tsv**.</br>

 |  |**Measure**|**Trait_length_after_dropna**| **signature_after_dropna**|**chi/fval**|**pval**|
 | :--------------------: |:--------------------:|:--------------------:|:--------:|:-------:|:--------------------:|
 | **sample 1**|f_oneway|int(more than threshold)|int|float|float|
 |...|...|...|...|...|...|
 | **sample m**|chisquare|int(less than threshold)|int|float|float|
