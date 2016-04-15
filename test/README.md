# PLL test framework

  src/    - source of test cases
  obj/    - binaries of test cases (automatically generated)
  out/    - expected output
  result/ - output of failed tests

## Update test framework

How to create a new test:

1. Place test case source file %.c into src/

  Use a descriptive but short name. 
  You could also write a long description into this file

2. Add the filename to CFILES variable in the Makefile

  e.g.,
  
  ```
  CFILES = src/alpha-cats.c \
           src/blopt-minimal.c \
           ...
           src/newtest.c
  ```
  
3. Compile with the provided Makefile

  A binary with the same name (%) should be created in obj/

4. Validate it propertly!

  Every time you place a buggy test case, god kills a kitten
  (and we do not want that)

5. Pipe the output into out/%.out

6. Proceed to next section and verify it matches the output and there are no leaks

## Run the test framework

1. Run `make` for compiling the test cases
2. Run `runtest.py` for executing the tests
3. Check the output for errors
4. If any of the tests fail, go and fix your code!
   The wrong output will be placed in result/

## Speed test

By default, runtest.py executes the validation test.
The same as 'runtest.py validation'

Use 'runtest.py speed' for running the speed tests

## Evaluating a subset

Use 'runtest.py validation|speed test1 test2 .. testN' for evaluating a
subset of test cases.

e.g., ./runtest.py speed hky alpha-cats

## Build tests for Windows

1. Build the library dll file and place them in current directory

2. Run `make -f Makefile.w64`

3. Run `./testwin.sh` for testing

4. If some dll missing files are reported, locate them and copy them
   in current or obj/ directory.
   
## Check errors

If a test fails, it will produce an error output in results directory. 
The generated files have the following format:

`testfail_ATTR_TESTNAME_DATETIME` and `testfail_ATTR_TESTNAME_DATETIME.err`

Where:
  * `ATTR` is the set of attributes used: `A` for AVX, `C` for CPU (non-vectorized version), `T` for tip vectors.
  * `TESTNAME` is the name of the test (e.g., alpha-cats)
  * `DATETIME` is the date and time when the error failed with format YYYYMMDDhhmmss
  * `.err` file is the error stream output. Most times it will be a blank file.

Tests output may be very complex if they print out, for example, the CLVs. The
best way to check what failed in the test is usually to compare the output with
the expected one:

e.g., `diff testfail_ATTR_TESTNAME_DATETIME out/TESTNAME.out`
      `vimdiff testfail_ATTR_TESTNAME_DATETIME out/TESTNAME.out`

