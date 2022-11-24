make clean
make
make runlocal

# if result.txt is equal to result_100_1000_sin.txt then echo SUCCES else echo FAIL
if cmp -s "result.txt" "result_100_1000_sin.txt"; then
    echo "SUCCES"
else
    echo "FAIL - Not equal"
fi