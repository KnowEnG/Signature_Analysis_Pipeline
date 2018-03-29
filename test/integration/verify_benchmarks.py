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

                accuracy = verify_accuracy(RESULT)
                if round(accuracy, 6) == BENCHMARK_name_list[1]:
                    print(BENCHMARK_rst,'\t\t', '______ Accuracy PASS ______' )
                else:
                    print(BENCHMARK_rst,'\t\t', '****** Accuracy FAIL ******' )

                if filecmp.cmp(RESULT, BENCHMARK_rst) == True:
                    print(BENCHMARK_rst,'\t\t', '______ PASS ______')
                else:
                    print(BENCHMARK_rst,'\t\t', '****** FAIL ******')

            # run checks on best_match_xxx files
            BENCHMARK_best_match = 'best_match_' + BENCHMARK_name
            if BENCHMARK_best_match in f:
                RESULT = os.path.join(results_dir, f)
                BENCHMARK_bmt = os.path.join(verification_dir, BENCHMARK_best_match + ".tsv")
                if filecmp.cmp(RESULT, BENCHMARK_bmt) == True:
                    print(BENCHMARK_bmt, '\t\t', '______ PASS ______')
                else:
                    print(BENCHMARK_bmt,'\t\t', '****** FAIL ******')
    return

def verify_accuracy(file_name):
    """ Calculate accuracy given similarity dataframe and benchmark result

    Args:
        file_name: result dataframe file_name from run_similarity methods
    """
    similarity_df        = pd.read_csv(file_name, index_col=0, header=0, sep='\t')
    similarity_names     =  similarity_df.columns
    similarity_names     = [i.split('.')[0] for i in similarity_names]
    similarity_df.columns= similarity_names
    result    = similarity_df.idxmax(axis=1, skipna=True)
    benchmark = pd.read_csv('../data/spreadsheets/label_validation.txt', index_col=None, header=None, sep='\t')
    ret_li    = result.values
    ben_li    = benchmark.values.reshape((1,-1))[0]

    common = ret_li==ben_li
    common[common==True] = 1
    common[common==False] = 0
    accuracy = sum(common)/len(ret_li)
    return accuracy

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
