
LDFLAGS=-lcunit

run_test: test_time_utils test_file_utils test_string_utils test_array_utils test_fortran_utils
	./test_time_utils
	./test_file_utils 
	./test_string_utils
	./test_array_utils
	./test_fortran_utils

test_time_utils: test_time_utils.c ../time_utils.o
	$(CC) $^ -lcunit -o $@ 

test_file_utils: test_file_utils.c ../file_utils.o
	$(CC) $^ -lcunit -o $@ 

test_string_utils: test_string_utils.c ../string_utils.o
	$(CC) $^ -lcunit -o $@ 

test_array_utils: test_array_utils.c ../array_utils.o
	$(CC) $^ -lcunit -o $@ 

test_fortran_utils: test_fortran_utils.c ../fortran_utils.o
	$(CC) $^ -lcunit -o $@ 
	
trackmem_test: trackmem_test.c ../trackmem_test.o
	$(CC) $^ -lcunit -o $@ 
	
clean:
	-rm -f test_time_utils test_file_utils test_string_utils test_array_utils test_fortran_utils trackmem_test
