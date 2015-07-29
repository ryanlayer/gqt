#!/bin/bash
echo "+---------------------+"
echo "| functional_tests.sh |"
echo "+---------------------+"
./functional_tests.sh

echo
echo "+------------------+"
echo "| io_error_test.sh |"
echo "+------------------+"
./io_error_test.sh

echo 
echo "+------------+"
echo "| io_test.sh |"
echo "+------------+"
./io_test.sh

echo
echo "+----------------------+"
echo "| usage_error_tests.sh |"
echo "+----------------------+"
./usage_error_tests.sh
