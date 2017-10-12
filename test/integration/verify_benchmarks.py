"""
sobh@illinois.edu
"""

import os
import filecmp

verification_dir = '../data/verification'
results_dir      = '../test/run_dir/results'

def verify_benchmark(measure, BENCHMARK_name_list, BENCHMARK_YML) :

    run_command  = 'python3 ../src/gene_signature.py -run_directory ./run_dir -run_file ' + BENCHMARK_YML

    os.system(run_command)

    All_files_in_results_dir = os.listdir(results_dir)

    for f in All_files_in_results_dir:
        for BENCHMARK_name in BENCHMARK_name_list:

          if BENCHMARK_name in f :

              RESULT    = os.path.join(results_dir     , f                      )
              BENCHMARK = os.path.join(verification_dir, BENCHMARK_name + '.tsv')

              if filecmp.cmp(RESULT, BENCHMARK) == True:
                  print(BENCHMARK,'\t\t', '______ PASS ______' )
              else:
                  print(BENCHMARK,'\t\t', '****** FAIL ******' )

    return

def main():
    BENCHMARK = {'spearman': 
                           { 'similarity'      :[ 
                                      'BENCHMARK_1_GS_spearman.yml'
                                    , 'result_similarity_spearman'        ] 

                           ,'net_similarity'   :[  
                                      'BENCHMARK_2_GS_net_spearman.yml'
                                    , 'result_net_similarity_spearman'    ]

                           ,'cc_similarity'    :[  
                                      'BENCHMARK_3_GS_cc_spearman.yml'
                                    , 'result_cc_similarity_spearman'     ]
                           ,'cc_net_similarity':[  
                                      'BENCHMARK_4_GS_cc_net_spearman.yml'
                                    , 'result_cc_net_similarity_spearman' ]
                           }, 
                 'cosine': 
                           { 'similarity'      :[ 
                                      'BENCHMARK_1_GS_cos.yml'
                                    , 'result_similarity_cosine'        ] 

                           ,'net_similarity'   :[  
                                      'BENCHMARK_2_GS_net_cos.yml'
                                    , 'result_net_similarity_cosine'    ]

                           ,'cc_similarity'    :[  
                                      'BENCHMARK_3_GS_cc_cos.yml'
                                    , 'result_cc_similarity_cosine'     ]
                           ,'cc_net_similarity':[  
                                      'BENCHMARK_4_GS_cc_net_cos.yml'
                                    , 'result_cc_net_similarity_cosine' ]
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
