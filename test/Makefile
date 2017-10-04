SCRIPT = ../src/gene_signature.py
RUN_DIR = ./run_dir
RUN_FILES = ../data/run_files
RESULTS_DIR = $(RUN_DIR)/results
INPUT_DATA_DIR = ../data
UNIT_TEST_TEMP_DIR = ./unit/tmp

RMR = rm -rf
MKDIR_P = mkdir -p

# ----------------------------------------------------------------
# Benchmarks setup Section                                       -
# ----------------------------------------------------------------
env_setup: create_run_dir copy_run_files

create_run_dir:
	$(MKDIR_P) $(RESULTS_DIR) 

copy_run_files:
	cp $(RUN_FILES)/*.yml $(RUN_DIR) 

clean_dir_recursively:
	$(RMR) $(RUN_DIR)

# ----------------------------------------------------------------
# Large Benchmarks RUN Section                                   -
# ----------------------------------------------------------------
run_similarity:
	python3 $(SCRIPT) -run_directory $(RUN_DIR) -run_file BENCHMARK_1_GS_spearman.yml

run_net_similarity:
	python3 $(SCRIPT) -run_directory $(RUN_DIR) -run_file BENCHMARK_2_GS_net_spearman.yml

run_cc_similarity:
	python3 $(SCRIPT) -run_directory $(RUN_DIR) -run_file BENCHMARK_3_GS_cc_spearman.yml

run_cc_net_similarity:
	python3 $(SCRIPT) -run_directory $(RUN_DIR) -run_file BENCHMARK_4_GS_cc_net_spearman.yml

# ----------------------------------------------------------------
# Small Benchmarks RUN Section                                   -
# ----------------------------------------------------------------

run_cos_test:
	python3 $(SCRIPT) -run_directory $(RUN_DIR) -run_file TEST_1_GS_cos.yml

run_net_cos_test:
	python3 $(SCRIPT) -run_directory $(RUN_DIR) -run_file TEST_2_GS_net_cos.yml

run_cc_cos_test:
	python3 $(SCRIPT) -run_directory $(RUN_DIR) -run_file TEST_3_GS_cc_cos.yml

run_cc_net_cos_test:
	python3 $(SCRIPT) -run_directory $(RUN_DIR) -run_file TEST_4_GS_cc_net_cos.yml

run_cos_bench:
	python3 $(SCRIPT) -run_directory $(RUN_DIR) -run_file BENCH_1_GS_cos.yml

# ----------------------------------------------------------------
# - VERIFICATION TESTS RUN SECTION                               -
# ----------------------------------------------------------------
verification_tests:
	python3 ./integration/verify_benchmarks.py

# ----------------------------------------------------------------
# - CLEAN UP SECTION                                             -
# ----------------------------------------------------------------	
unit_test_cleanup:
	$(RMR) $(UNIT_TEST_TEMP_DIR)


