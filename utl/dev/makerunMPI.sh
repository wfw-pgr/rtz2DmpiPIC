#!/bin/sh

job=`grep -E "::\s+job\s+=" variables.f90 | cut -d "'" -f2`
omp=`grep -E "::\s+OMPNumThreads\s+=" variables.f90 | cut -d "=" -f2`
PEt=`grep -E "::\s+PEtot\s+=" variables.f90 | cut -d "=" -f2`
jobname=`echo ${job}`
omp=`echo ${omp}`
PEt=`echo ${PEt}`

FILENAME="run.sh"
export OMP_NUM_THREADS=${omp}
echo "OMP_NUMTHREADS == "${omp}

echo "make "${FILENAME}
# printf "#!/bin/sh\n\n">${FILENAME}
# printf "#PBS -V\n">>${FILENAME}
# printf "#PBS -l ncpus=6\n">>${FILENAME}
# printf "#PBS -N "${job}"\n">>${FILENAME}
# printf "#PBS -j oe\n">>${FILENAME}
# printf "#PBS -o stderr.log\n">>${FILENAME}
# printf "#PBS -q alfven\n\n">>${FILENAME}

# printf "export OMP_NUM_THREADS="${omp}"\n">>${FILENAME}
# printf 'cd $PBS_O_WORKDIR\n'>>${FILENAME}
# printf "mpiexec -n "${PEt}" ./main\n">>${FILENAME}
# printf "python3 src/sendGmail.py "${job}"\n">>${FILENAME}
printf "#!/bin/sh\n\n">${FILENAME}
printf "export OMP_NUM_THREADS="${omp}"\n">>${FILENAME}
# printf "mpiexec -machinefile hosts -n "${PEt}" ./main">>${FILENAME}
printf "mpiexec -n "${PEt}" ./main">>${FILENAME}
printf " > job/"${job}"/run.log\n">>${FILENAME}
printf "python3 src/sendGmail.py "${job}"\n">>${FILENAME}

mv run.sh ../
mkdir -p ../job/${job}
cp -r ../src/ ../job/${job}/src/
