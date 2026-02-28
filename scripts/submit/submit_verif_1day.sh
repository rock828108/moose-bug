#!/bin/bash
# Short validation run (1 day default) — align with reference before full 12-year run.
# Compare with 1-rank Dirac version (sleipner_corrected_v2_immiscible_1M.i): pressure range, co2_mass, max_sgas.
#
# DIAGNOSTIC: If segfault in Setting Up, use 20×1 first: OMP_NUM_THREADS=1 bash ...
#
# How to submit on the cluster (examples):
#   1) 20 MPI × 4 threads (default):
#        salloc -N 1 -n 20 -t 2:00:00 --mem=400G
#        cd /home/yifan/paper3/sleipner/moose_input
#        bash scripts/submit/submit_verif_1day.sh 2>&1 | tee logs/verif_1day_20x4.log
#
#   2) 20 MPI × 2 threads (less memory per rank):
#        salloc -N 1 -n 20 -t 2:00:00 --mem=400G
#        cd /home/yifan/paper3/sleipner/moose_input
#        OMP_NUM_THREADS=2 bash scripts/submit/submit_verif_1day.sh 2>&1 | tee logs/verif_1day_20x2.log
#
#   3) 10 MPI × 4 threads (requires 10-way split already generated):
#        salloc -N 1 -n 10 -t 2:00:00 --mem=400G
#        cd /home/yifan/paper3/sleipner/moose_input
#        bash scripts/submit/submit_verif_1day.sh 2>&1 | tee logs/verif_1day_10x4.log
#
#   4) 10 MPI × 2 threads:
#        salloc -N 1 -n 10 -t 2:00:00 --mem=400G
#        cd /home/yifan/paper3/sleipner/moose_input
#        OMP_NUM_THREADS=2 bash scripts/submit/submit_verif_1day.sh 2>&1 | tee logs/verif_1day_10x2.log
#
# Watching logs in another terminal:
#   tail -f logs/verif_1day_20x4.log          # or the corresponding log file name
# For line-by-line log writes when piping:
#   stdbuf -oL -eL bash scripts/submit/submit_verif_1day.sh 2>&1 | stdbuf -oL tee logs/verif_1day_20x4.log
#
# Change VERIF_END_TIME to test longer intervals:
#   VERIF_END_TIME=604800  bash ...   # 1 week
#   VERIF_END_TIME=2592000 bash ...   # 1 month
set -e
TOTAL_TASKS="${SLURM_NTASKS:-20}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
cd "$PROJECT_ROOT"

# Same setup as presplit script
module purge 2>/dev/null || true
source /etc/profile.d/modules.sh 2>/dev/null || true
module load OpenMPI/4.1.1-GCC-11.2.0 2>/dev/null || true
source ~/miniforge3/etc/profile.d/conda.sh
conda activate moose-py2
export LD_LIBRARY_PATH="$HOME/projects/moose/lib:$HOME/projects/moose/libmesh/lib:${CONDA_PREFIX}/lib:${LD_LIBRARY_PATH}"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-4}"
export OMP_PROC_BIND=close
export OMP_PLACES=cores

MOOSE_EXE="$HOME/projects/moose/modules/porous_flow/porous_flow-opt"
INPUT="inputs/sleipner_corrected_v2_immiscible_1M_20rank_sink_fast_presplit_verif.i"
# SPLIT_DIR="${PROJECT_ROOT}/meshes/sleipner_grid_1M.cpa.gz"

# Baked mesh required (no pre-split): run scripts/submit/run_mesh_bake.sh first
BAKED_MESH="${PROJECT_ROOT}/meshes/sleipner_grid_1M_with_well.e"
[ -f "$BAKED_MESH" ] || BAKED_MESH="${PROJECT_ROOT}/../meshes/sleipner_grid_1M_with_well.e"
[ ! -f "$BAKED_MESH" ] && echo "ERROR: Baked mesh not found. Run first: bash scripts/submit/run_mesh_bake.sh" && exit 1
[ ! -f "$MOOSE_EXE" ] && echo "ERROR: $MOOSE_EXE not found" && exit 1
[ ! -f "$INPUT" ] && echo "ERROR: $INPUT not found" && exit 1

mkdir -p results/exodus/verif_1day logs

# Override end_time: 86400 = 1 day; 604800 = 1 week; 2592000 = 1 month
END_TIME="${VERIF_END_TIME:-86400}"
DAYS=$((END_TIME / 86400))
echo "Short validation: end_time = ${END_TIME} s (${DAYS} days), ranks = ${TOTAL_TASKS}"
echo "Input file: ${INPUT}  (${TOTAL_TASKS}×${OMP_NUM_THREADS}, baked mesh, SMP_asm1, dt=0.01)"
echo "Compare CSV: min/max pwater, co2_mass, max_sgas vs 1-rank Dirac (sleipner_corrected_v2_immiscible_1M.i)"
echo "========================================================"

# Verif input has end_time=86400 and verif output paths; override end_time only when VERIF_END_TIME is set.
# No --use-split: FileMeshGenerator in input loads split directly. Ranks must match split (20).
EXTRA=""
[ "${END_TIME}" != "86400" ] && EXTRA="Executioner/end_time=${END_TIME}"
srun -n ${TOTAL_TASKS} --cpu-bind=cores stdbuf -oL -eL \
  "$MOOSE_EXE" -i "$INPUT" $EXTRA --n-threads=${OMP_NUM_THREADS} --distributed-mesh --timing

echo ""
echo "Results: results/exodus/verif_1day/verif_1day.csv"
echo "Sanity: total_injected_mass ≈ 28.5 * time; max_sgas < 1.1; pwater in ~MPa range"
