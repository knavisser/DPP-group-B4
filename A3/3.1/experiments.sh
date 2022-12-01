module load prun
module load openmpi/gcc/64
#module load cuda11.1/toolkit/11.1.1

make clean
make

for i in {1..10}
do
    prun -v -np 1 -native "-C Titan --gres=gpu:1" ./assign3_1 1000 100
done

for i in {1..10}
do
    prun -v -np 1 -native "-C Titan --gres=gpu:1" ./assign3_1 10000 100
done

for i in {1..10}
do
    prun -v -np 1 -native "-C Titan --gres=gpu:1" ./assign3_1 100000 100
done

for i in {1..10}
do
    prun -v -np 1 -native "-C Titan --gres=gpu:1" ./assign3_1 1000000 100
done

for i in {1..10}
do
    prun -v -np 1 -native "-C Titan --gres=gpu:1" ./assign3_1 10000000 100
done