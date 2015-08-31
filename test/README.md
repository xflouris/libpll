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

2. Compile with the provided Makefile

  A binary with the same name (%) should be created in obj/

3. Validate it propertly!

  Every time you place a buggy test case, god kills a kitten
  (and we do not want that)

4. Pipe the output into out/%.out

5. Proceed to section B and verify it matches the output and there are no leaks

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
