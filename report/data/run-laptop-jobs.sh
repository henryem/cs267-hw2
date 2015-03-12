export SUBMIT_JOB=$(dirname $0)/submit-job

# serial
for NUM_PARTICLES in 250 1000 4000 16000 64000 256000; do
  $SUBMIT_JOB -m laptop -j serial -n $NUM_PARTICLES
done

# CUDA
for NUM_PARTICLES in 250 1000 4000 16000 64000 256000; do
  $SUBMIT_JOB -m laptop -j gpu -n $NUM_PARTICLES
done