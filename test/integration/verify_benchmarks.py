"""
sobh@illinois.edu
"""

import os
import filecmp
import pandas as pd

verification_dir = '../data/verification'
results_dir      = '../test/run_dir/results'

def verify_benchmark(measure, BENCHMARK_name_list, BENCHMARK_YML) :

    run_command  = 'python3 ../src/gene_signature.py -run_directory ./run_dir -run_file ' + BENCHMARK_YML
    os.system(run_command)

    All_files_in_results_dir = os.listdir(results_dir)
    for f in All_files_in_results_dir:
        
        for BENCHMARK_name in BENCHMARK_name_list[:1]:
            # run checks on result_xxx files
            BENCHMARK_result = 'result_' + BENCHMARK_name
            if BENCHMARK_result in f :
                RESULT    = os.path.join(results_dir     , f                      )
                BENCHMARK_rst = os.path.join(verification_dir, BENCHMARK_result + '.tsv')

                if filecmp.cmp(RESULT, BENCHMARK_rst) == True:
                    print(BENCHMARK_rst,'\t\t', '______ PASS ______')
                else:
                    print(BENCHMARK_rst,'\t\t', '****** FAIL ******')

            # run checks on Gene_to_TF_Association_xxx files
            BENCHMARK_Gene_to_TF_Association = 'Gene_to_TF_Association_' + BENCHMARK_name
            if BENCHMARK_Gene_to_TF_Association in f:
                RESULT = os.path.join(results_dir, f)
                BENCHMARK_bmt = os.path.join(verification_dir, BENCHMARK_Gene_to_TF_Association + ".tsv")
                if filecmp.cmp(RESULT, BENCHMARK_bmt) == True:
                    print(BENCHMARK_bmt, '\t\t', '______ PASS ______')
                else:
                    print(BENCHMARK_bmt,'\t\t', '****** FAIL ******')
    return

def main():
    BENCHMARK = {'spearman': 
                           { 'similarity'      :[ 
                                      'BENCHMARK_1_GS_spearman.yml'
                                    , 'similarity_spearman'
                                    , 0.801527]
                           ,'net_similarity'   :[  
                                      'BENCHMARK_2_GS_net_spearman.yml'
                                    , 'net_similarity_spearman'
                                    , 0.793893]
                           ,'cc_similarity'    :[  
                                      'BENCHMARK_3_GS_cc_spearman.yml'
                                    , 'cc_similarity_spearman'
                                    , 0.839695]
                           ,'cc_net_similarity':[  
                                      'BENCHMARK_4_GS_cc_net_spearman.yml'
                                    , 'cc_net_similarity_spearman'
                                    , 0.793893]
                           }, 
                 'pearson':
                           { 'similarity'      :[
                                      'BENCHMARK_1_GS_pearson.yml'
                                    , 'similarity_pearson'
                                    , 0.0]
                           ,'net_similarity'   :[
                                      'BENCHMARK_2_GS_net_pearson.yml'
                                    , 'net_similarity_pearson'
                                    , 0.0]
                           ,'cc_similarity'    :[
                                      'BENCHMARK_3_GS_cc_pearson.yml'
                                    , 'cc_similarity_pearson'
                                    , 0.0]
                           ,'cc_net_similarity':[
                                      'BENCHMARK_4_GS_cc_net_pearson.yml'
                                    , 'cc_net_similarity_pearson'
                                    , 0.0]
                           },
                 'cosine': 
                           { 'similarity'      :[ 
                                      'BENCHMARK_1_GS_cos.yml'
                                    , 'similarity_cosine'
                                    , 0.648855]
                           ,'net_similarity'   :[  
                                      'BENCHMARK_2_GS_net_cos.yml'
                                    , 'net_similarity_cosine'
                                    , 0.641221]
                           ,'cc_similarity'    :[  
                                      'BENCHMARK_3_GS_cc_cos.yml'
                                    , 'cc_similarity_cosine'
                                    , 0.641221]
                           ,'cc_net_similarity':[  
                                      'BENCHMARK_4_GS_cc_net_cos.yml'
                                    , 'cc_net_similarity_cosine'
                                    , 0.671756]
                          }
                }

    os.system('make env_setup')
    for measure in BENCHMARK.keys():
        for key in BENCHMARK[measure].keys(): 
            BENCHMARK_list = BENCHMARK[measure][key]
            BENCHMARK_YML  = BENCHMARK_list[0]
            verify_benchmark(measure, BENCHMARK_list[1:], BENCHMARK_YML)
            os.system('rm ./run_dir/results/*')

if __name__ == "__main__":
    main()
