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
run_cos:
	python3 $(SCRIPT) -run_directory $(RUN_DIR) -run_file TEST_1_GS_cos.yml

# ----------------------------------------------------------------
# - CLEAN UP SECTION                                             -
# ----------------------------------------------------------------	
unit_test_cleanup:
	$(RMR) $(UNIT_TEST_TEMP_DIR)


