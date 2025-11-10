#!/bin/bash -l
#SBATCH -o ./slurm.%j.out
#SBATCH -e ./slurm.%j.err
#SBATCH -D ./
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --gres=gpu:A100:4
#SBATCH --mem=236gb
#SBATCH --partition=gpu-single
#SBATCH --mail-type=ALL
#SBATCH --mail-user=youremail
#SBATCH --time=120:00:00
#SBATCH --job-name=XXX

# Load modules
module purge
module load chem/gromacs/2025.2-cuda-12.9

export OMP_NUM_THREADS=16
export GMX_ENABLE_DIRECT_GPU_COMM=1
export GMX_GPU_PME_DECOMPOSITION=1

MDRUN="gmx mdrun"
GMXBIN="gmx"

##############################################################################

input_name=$1
nvt_name=$2
npt_name=$3
pull_name=$4
phase=$5        # "nvt", "npt", or "pull"


varname="${phase}_name"
checkpoint="./${!varname}.cpt"
tpr="./${!varname}.tpr"

echo "Starting phase: $phase"
echo "Input: $input_name"
echo "Checkpoint: $checkpoint"

# Generate TPR if not already present
if [[ ! -f "$checkpoint" ]]; then
  echo "Generating $phase.tpr..."
  case $phase in
    nvt)
      $GMXBIN grompp -f nvt.mdp -c "$input_name.gro" -r "$input_name.gro" -o "$nvt_name.tpr" -p "$input_name.top" -maxwarn 2
      ;;
    npt)
      $GMXBIN grompp -f npt.mdp -c "$nvt_name.gro" -r "$nvt_name.gro" -t "$nvt_name.cpt" -o "$npt_name.tpr" -p "$input_name.top" -maxwarn 2
      ;;
    pull)
      $GMXBIN grompp -f pull.mdp -c "$npt_name.gro" -t "$npt_name.cpt" -o "$pull_name.tpr" -p "$input_name.top" -n "$input_name.ndx" -maxwarn 2
      ;;
    *)
      echo "Unknown phase: $phase"
      exit 1
      ;;
  esac
fi

PME_RANKS=1

EXTRA_FLAGS=""
if [[ "$phase" == "pull" ]]; then
  EXTRA_FLAGS="-px ${!varname}_pullx.xvg -pf ${!varname}_pullf.xvg"
fi

SECONDS=0
srun -n 1 --gpus=4 $MDRUN \
  -deffnm "${!varname}" \
  -s "$tpr" \
  -cpi "$checkpoint" \
  -ntmpi 4 \
  -ntomp 16 \
  -dlb no \
  -nb gpu \
  -pme gpu \
  -npme $PME_RANKS \
  -bonded gpu \
  -update gpu \
  -pin on \
  -maxh 119.90 \
  $EXTRA_FLAGS

# Safety check to avoid short resubmissions
runtime=$SECONDS

min_seconds=60
if (( runtime < min_seconds )); then
  echo "Runtime < ${min_seconds}s (${runtime}s). Not resubmitting to avoid loops."
  exit 42
fi

# Resubmit based on steps compared from mdp and log file
mdp_file="${phase}.mdp"
log_file="${!varname}.log"

nsteps=$(grep -i '^[[:space:]]*nsteps[[:space:]]*=' "$mdp_file" | sed 's/;.*//' | awk -F '=' '{print $2}' | tr -d ' ')
checkpoint_found=$(grep -q "^Writing checkpoint, step $nsteps" "$log_file"; echo $?)
fatal_error_found=$(grep -q "Fatal error" "$log_file"; echo $?)

if [[ $checkpoint_found -eq 0 ]]; then
  echo "Checkpoint for final step $nsteps found. Proceeding to next phase..."
  case $phase in
    nvt)
      sbatch "$0" "$input_name" "$nvt_name" "$npt_name" "$pull_name" npt
      ;;
    npt)
      sbatch "$0" "$input_name" "$nvt_name" "$npt_name" "$pull_name" pull
      ;;
    pull)
      echo "Pulling simulation complete. Not resubmitting further."
      ;;
  esac
elif [[ $fatal_error_found -eq 0 ]]; then
  echo "Fatal error detected in $log_file. Not resubmitting $phase phase."
else
  echo "Final checkpoint step $nsteps NOT found in log AND no fatal error detected. Resubmitting $phase phase..."
  sbatch "$0" "$input_name" "$nvt_name" "$npt_name" "$pull_name" "$phase"
fi
