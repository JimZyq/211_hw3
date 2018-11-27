#! /bin/bash 
rm -f P3_*

echo "Loading modules..."
module purge
module load gcc/gcc-4.7.2
module load mpich-3.0.3/gcc-4.7.2
echo "Modules loaded! Now compiling..."
set -x
mpic++ -o run1 source1.cpp -std=c++11 -Ofast -march=native
mpic++ -o run2 source2.cpp -std=c++11 -Ofast -march=native
mpic++ -o run3 source3.cpp -std=c++11 -Ofast -march=native
set +x
echo "Source code compiled"

main = "P3_"
part = "P"
node = "N"
ext = ".sh"

for n in 1 2 4 8
do
    for i in 1 2 3
    do
        fn = $main$part$i$node$n$ext
        echo "#! /bin/bash" > $fn
        echo "#SBATCH --job-name=P"$i"N$n" >> $fn
        echo "#SBATCH --ouptut=result"$i"$n.txt" >> $fn
        echo "#SBATCH -N $n" >> $fn
        echo "#SBATCH -t 00:20:00" >> $fn
        echo "" >> $fn
        echo "mpirun ./run$i 10000000000" >> $fn
        printf "."
    done
done
printf "\n"
