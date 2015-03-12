export SUBMIT_JOB=$(dirname $0)/submit-job

# CUDA-naive
for NUM_PARTICLES in 250 1000 4000 16000 64000; do
  $SUBMIT_JOB -m stampede -j gpu-naive -n $NUM_PARTICLES
done

# CUDA
for NUM_PARTICLES in 250 1000 4000 16000 64000 256000 1024000 4096000; do
  $SUBMIT_JOB -m stampede -j gpu -n $NUM_PARTICLES
done