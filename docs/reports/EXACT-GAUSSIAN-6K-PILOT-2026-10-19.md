# Exact-Gaussian 6,000-replication resource pilot (2026-10-19)

**Status:** one-task benchmark completed; no full array submitted. I6 remains FAIL, I7 remains stopped, and the pmed submission hold remains in force. This pilot measures resources and output completeness, not coverage or independent formula correctness.

## Reproducibility and result

- Hopper job `4342252_300` ran the s6, N=1,000 task: 200 outer replications and 2,000 nonparametric bootstrap refits per replication (400,000 refits). `sacct` reported `COMPLETED`, exit `0:0`, 37 seconds elapsed, 8 allocated CPUs, and `894716K` peak batch RSS under an 8 GB request. The task log reported zero failed refits.
- `analyze.R` was run from `/users/dtofighi/probmed_exact_20261019_pilot`, where it sources `code/config.R`. It verified 200 unique replication IDs, 2,000 attempted bootstrap refits each, zero failed refits, and finite ordered intervals.
- The staged source was under `/users/dtofighi/probmed_exact_20261019_pilot/code/`; `sbatch` used `--chdir` to that directory and absolute `PMED_OUT_DIR=/users/dtofighi/probmed_exact_20261019_pilot`. Local and staged SHA-256 hashes matched: `runner.R` `e38d50f7dc2855105149ba9a6f611433ed3cb12bd9b84ed69de8fab22a759e49`, `config.R` `fc31d5e4c600ad35b9f529c0913d33d3418b68e96b10b2fcfcab1880851c8862`, `submit_exact.sh` `5b62a9629869875606eab1dd907036774a56d4f60f1cf6d5d47ddc14bab9df4b`, and `analyze.R` `b5a7560cbde60caa217662cf88396c0a5e63e6e1c4f07ba3d6746b525994666c`.
- The copied pilot output's SHA-256 matched Hopper: `88dfda48110e9aebd02bc93ff4aeedf9d722e44f34b3e7e300d7d833fc7d2365` (`parts/e_300.csv`). The local copy is `/private/tmp/probmed-exact-pilot-20261019/e_300_hopper.csv`.
- The pilot runner used seed `20261020 + task`. Exclude its output from prospective coverage analysis and use a different, frozen seed family (proposed `20261021`) for every full-run task.

## Planning estimate and limits

The proposed 20-cell design at 6,000 outer replications and 2,000 resamples requires 120,000 outer datasets and 240 million bootstrap refits. With 200 outer replications per task, it would have 600 array tasks. Scaling the 37-second, largest-N pilot to every task gives **6.2 node-hours or 49.3 allocated-core-hours**, about 23 minutes of execution at 16 concurrent tasks. A 2× planning allowance is **12.3 node-hours or 98.7 core-hours**, about 46 minutes at 16 concurrent tasks. These are extrapolations, not guaranteed limits or queue-time estimates. The existing I6 Gaussian tasks took 194–271 seconds each, but they also performed the Monte Carlo candidate and other work, so they are a conservative comparison rather than a direct estimate of this exact-only runner.

At 13:54 MDT, `squeue -u dtofighi -r` showed 148 pending and 4 running tasks. A 600-task array may exceed Hopper's submitted-task cap; split ranges only after checking live capacity. No full-array submission follows from this benchmark.

Before a full run, freeze the 6,000-replication claim gate and manifest; test the formula against an independent oracle, run the I2/I6 replay and failure-injection checks in the proposal, and include all-task completeness and all-attempt coverage validation. The one-task pilot does not satisfy those scientific gates.
