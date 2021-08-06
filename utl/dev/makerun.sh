#!/bin/sh

job=`grep -E "::\s+job\s+=" variables.f90 | cut -d "'" -f2`
omp=`grep -E "::\s+OMPThreadNum\s+=" variables.f90 | cut -d "=" -f2`
jobname=`echo ${job}`
omp=`echo ${omp}`
FILENAME="run.sh"
export OMP_NUM_THREADS=${omp}
echo "OMP_NUM_THREADS == "${omp}

echo "make "${FILENAME}
printf "#!/bin/sh\n\n">${FILENAME}
printf "#PBS -V\n">>${FILENAME}
printf "#PBS -N "${job}"\n">>${FILENAME}
printf "#PBS -j oe\n">>${FILENAME}
printf "#PBS -o stderr.log\n">>${FILENAME}
printf "#PBS -q stix\n\n">>${FILENAME}

printf "export OMP_NUM_THREADS="${omp}"\n">>${FILENAME}
printf 'cd $PBS_O_WORKDIR\n'>>${FILENAME}
printf "mpiexec -n 2 ./main > job/"${job}"/run.log\n">>${FILENAME}
printf "python3 src/sendGmail.py "${job}"\n">>${FILENAME}
mv run.sh ../
mkdir -p ../job/${job}
cp -r ../src/ ../job/${job}/src/
