"""
sobh@illinois.edu

"""

import os
import filecmp

verification_dir = '../data/verification'
results_dir      = '../test/run_dir/results'

def verify_benchmark(BENCHMARK_name_list, BENCHMARK_YML) :

    run_command  = 'python3 ../src/gene_signature.py -run_directory ./run_dir -run_file ' + BENCHMARK_YML
    os.system(run_command)

    All_files_in_results_dir = os.listdir(results_dir)

    for f in All_files_in_results_dir:
        for BENCHMARK_name in BENCHMARK_name_list:
          if BENCHMARK_name in f :
              RESULT    = os.path.join(results_dir,      f             )
              BENCHMARK = os.path.join(verification_dir, BENCHMARK_name+'.tsv')
              if filecmp.cmp(RESULT, BENCHMARK) == True:
                  print(BENCHMARK, 'PASS' )
              else:
                  print(BENCHMARK, 'FAIL' )

def main():
    BENCHMARK = {'cos':      [ 
                                 'BENCHMARK_1_GS_cos.yml'
                               , 'result_cos'
                               ] 
               ,'net_cos':   [  
                                 'BENCHMARK_2_GS_net_cos.yml'
                               , 'result_net_cos'
                               ]
               ,'cc_cos':    [  
                                 'BENCHMARK_3_GS_cc_cos.yml'
                               , 'result_cc_cos'
                               ]
               ,'cc_net_cos':[  
                                 'BENCHMARK_4_GS_cc_net_cos.yml'
                               , 'result_cc_net_cos'
                               ]
                }

    os.system('make env_setup')
    for key in BENCHMARK.keys(): 
        BENCHMARK_list = BENCHMARK[key]
        BENCHMARK_YML  = BENCHMARK_list[0]
        # for BENCHMARK_name in BENCHMARK_list[1:] :
        verify_benchmark(BENCHMARK_list[1:], BENCHMARK_YML)
        os.system('rm ./run_dir/results/*')

if __name__ == "__main__":
    main()
