#!/usr/bin/env bash
# Run the FULL BAliBASE RV100 benchmark across every aligner on your PATH.
#
# Safe to run and re-run. It:
#   * RESUMES  - skips (family, aligner) pairs already in balibase.csv,
#   * won't HANG - skips any external aligner exceeding --engine-timeout,
#   * self-RETRIES - if the process dies (crash/OOM), it resumes automatically.
# If you stop it (Ctrl-C), just run it again; it picks up where it left off.
#
#   bash benchmarks/run_full_balibase.sh
# or leave it running in the background and watch the log:
#   nohup bash benchmarks/run_full_balibase.sh > balibase_run.log 2>&1 &
#   tail -f balibase_run.log
#
set -u
cd "$(dirname "$0")/.."            # repo root (holds RV100/, aligned_out/, benchmarks/)

REF='RV100/*.xml'
OUT='balibase.csv'
ALN='aligned_out'
TIMEOUT=600                        # seconds per external alignment before it is skipped
MAXMEM=10                          # GB (memory proxy): families up to this size use the accurate
                                   # consistency path and are included; larger ones are skipped.
                                   # NOTE: real peak RAM is ~2-2.5x this (so ~20-25 GB at 10), and the
                                   # consistency transform is O(N^3) in #sequences, so big families are
                                   # slow. Lower to 2-4 if the run is too slow or you see "Killed: 9".

for attempt in $(seq 1 50); do
  echo "=== attempt ${attempt}: $(date) ==="
  python benchmarks/run_benchmark.py --reference "$REF" --protein \
      --write-alignments "$ALN" --engine-timeout "$TIMEOUT" --max-mem-gb "$MAXMEM" --out "$OUT"
  code=$?
  if [ "$code" -eq 0 ]; then
    echo "=== COMPLETE: $(date) ==="
    break
  fi
  echo "=== exited (code ${code}); resuming in 5s (Ctrl-C to stop) ==="
  sleep 5
done

echo
echo "Done. Metrics: $OUT   Alignments for bali-score: $ALN/"
echo "Official scores:  python -c \"import csv,subprocess;[subprocess.run(['bali-score','-t',t,'-r',r]) for t,r,a in list(csv.reader(open('$ALN/bali_manifest.tsv'),delimiter=chr(9)))[1:]]\" | tee bali_scores.txt"
