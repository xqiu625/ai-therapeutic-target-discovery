#!/bin/bash
# Run all notebooks sequentially with SLURM dependencies
# Usage: ./run_all_sequential.sh

# Submit jobs with dependencies
JOB1=$(sbatch --parsable scripts/hpcc/submit_data_prep.sh)
echo "Submitted data_prep job: $JOB1"

JOB2=$(sbatch --parsable --dependency=afterok:$JOB1 scripts/hpcc/submit_foundation_models.sh)
echo "Submitted foundation_models job: $JOB2 (depends on $JOB1)"

JOB3=$(sbatch --parsable --dependency=afterok:$JOB2 scripts/hpcc/submit_ml_prediction.sh)
echo "Submitted ml_prediction job: $JOB3 (depends on $JOB2)"

JOB4=$(sbatch --parsable --dependency=afterok:$JOB2 scripts/hpcc/submit_grn_scenic.sh)
echo "Submitted grn_scenic job: $JOB4 (depends on $JOB2)"

JOB5=$(sbatch --parsable --dependency=afterok:$JOB2 scripts/hpcc/submit_trajectory.sh)
echo "Submitted trajectory job: $JOB5 (depends on $JOB2)"

# Wait for ML, GRN, and trajectory to complete before prioritization
JOB6=$(sbatch --parsable --dependency=afterok:$JOB3:$JOB4:$JOB5 scripts/hpcc/submit_prioritization.sh)
echo "Submitted prioritization job: $JOB6 (depends on $JOB3, $JOB4, $JOB5)"

JOB7=$(sbatch --parsable --dependency=afterok:$JOB6 scripts/hpcc/submit_cd52.sh)
echo "Submitted cd52 job: $JOB7 (depends on $JOB6)"

echo ""
echo "Pipeline submitted! Monitor with: squeue -u \$USER"
echo "Job dependency chain:"
echo "  data_prep -> foundation_models -> [ml_prediction, grn_scenic, trajectory] -> prioritization -> cd52"
