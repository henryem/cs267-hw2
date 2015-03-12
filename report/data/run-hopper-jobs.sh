export SUBMIT_JOB=$(dirname $0)/submit-job

# serial-naive
for NUM_PARTICLES in 100 200 400 800 1600 3200; do
  $SUBMIT_JOB -m hopper -j serial-naive -n $NUM_PARTICLES
done

# serial
for NUM_PARTICLES in 100 200 400 800 1600 3200 6400 12800 250 1000 4000 16000 64000 256000; do
  $SUBMIT_JOB -m hopper -j serial -n $NUM_PARTICLES
done

# openmp (O(n) scaling)
DEFAULT_OMP_PARALLELISM = 24
for NUM_PARTICLES in 250 1000 4000 16000 64000 256000; do
  $SUBMIT_JOB -m hopper -t $DEFAULT_OMP_PARALLELISM -p 1 -j openmp -n $NUM_PARTICLES
done

# openmp (weak scaling)
for NUM_PARTICLES in 250 1000 4000 16000 64000 256000; do
  for PARALLELISM in 1 2 4 8 16 24; do
    TOTAL_PARTICLES = $((NUM_PARTICLES*PARALLELISM))
    $SUBMIT_JOB -m hopper -t $PARALLELISM -p 1 -j openmp -n $TOTAL_PARTICLES
  done
done

# openmp (strong scaling)
for NUM_PARTICLES in 250 1000 4000 16000 64000 256000; do
  for PARALLELISM in 1 2 4 8 16 24; do
    TOTAL_PARTICLES = $((NUM_PARTICLES))
    $SUBMIT_JOB -m hopper -t $PARALLELISM -p 1 -j openmp -n $TOTAL_PARTICLES
  done
done