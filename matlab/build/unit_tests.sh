#!/bin/sh

# Type in the path of the matlab run time library/mcr thing here
matlab_mcr=/Applications/MATLAB/MATLAB_Compiler_Runtime/v715/

# run the 9 bus test
echo '---------------------------------------------'
echo 'Running the 9 bus test';
./run_cosmic.sh $matlab_mcr case9_ps.mat unit_test_case9_events.mat

# run the 39 bus test
echo '---------------------------------------------'
echo 'Running the 39 bus test';
./run_cosmic.sh $matlab_mcr case39_ps.mat unit_test_case39_events.mat

# done
echo '---------------------------------------------'
echo 'Unit tests completed. Check the output above'
echo 'to see if the results make sense'
echo '---------------------------------------------'


