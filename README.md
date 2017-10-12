# KnowEnG's Gene Signature Pipeline 
This is the Knowledge Engine for Genomics (KnowEnG), an NIH BD2K Center of Excellence, Signature Analysis Pipeline.

This pipeline performs network-based **signature analysis** on the columns of a given spreadsheet, where spreadsheet's columns correspond to sample-labels and rows correspond to gene-labels.  The signature is based on correlating gene expression data (network enriched) against known gene signature data.

There are four similarity "signature"  methods that one can choose from:

- similarity        (traditional method) 
- net_similarity    (with network enrichment)
- cc_similarity     (with bootstraps)
- cc_net_similarity (with bootstraps and network enrichment)


and two correlation measures:
- spearman 
- cosine 


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

### 6. Use one of the following "make" commands to select and run a similarity  option:


| **Command**              | **Option**                                              | 
|:-------------------      |:------------------------------------------------------- | 
| make run_spearman        | spearman similarity                                     |
| make run_net_spearman    | spearman similarity with network enrichment             |
| make run_cc_spearman     | spearman similarity with bootstraps                     |
| make run_cc_net_spearman | spearman similarity with bootstraps & network enrichment|

 
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

| **Key**                   | **Value**                                           | **Comments** |
| ------------------------- | --------------------------------------------------- | ------------ |
| method                    | **cos**, **cc_cos**, **net_cos** or **cc_net_cos**  | Choose similarity  method |
| gg_network_name_full_path | directory+gg_network_name                           | Path and file name of the 4 col network file |
| spreadsheet_name_full_path| directory+spreadsheet_name                          | Path and file name of user supplied gene sets |
| signature_name_full_path  | directory+signature_data_name                       | Path and file name of user supplied signature data |
| results_directory         | directory                                           | Directory to save the output files |
| tmp_directory             | directory                                           | Directory to save the intermediate files |
| rwr_max_iterations        | 100                                                 | Maximum number of iterations without convergence in random walk with restart |
| rwr_convergence_tolerence | 1.0e-8                                              | Frobenius norm tolerence of spreadsheet vector in random walk|
| rwr_restart_probability   | 0.7                                                 | alpha in `V_(n+1) = alpha * N * Vn + (1-alpha) * Vo` |
| rows_sampling_fraction    | 0.8                                                 | Select 80% of spreadsheet rows|
| number_of_bootstraps      | 4                                                   | Number of random samplings |
| processing_method         | serial or parallel or distribute                    | Choose processing method |
| parallelism               | number of cores to use in parallel processing       | Set number of cores for speed or memory |

gg_network_name = STRING_experimental_gene_gene.edge</br>
spreadsheet_name = ProGENI_rwr20_STExp_GDSC_500.rname.gxc.tsv</br>
signature_data_name = 

* * * 
## Description of Output files saved in results directory
* * * 

* Output files of all four methods save samples by signature heatmap with name **similarity_matrix_{method}_{timestamp}_viz.tsv**.</br>

 |                        |**signature 1**          |...                   |**signature m**       |
 | :--------------------: |:-----------------------:|:--------------------:|:--------------------:|
 | **sample 1**           |float                    |...                   |float                 |
 |...                     |...                      |...                   |...                   |
 | **sample n**           |float                    |...                   |float                 |
 
