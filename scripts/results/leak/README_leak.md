# Leak robustness (reviewer C2.3)

`leak_sweep.py` (in `scripts/`) re-runs, at a non-zero leak rate, exactly the draws that were
admissible in either model at gamma = 0.  It imports `fused_sweep` and reuses its cube,
`to_param`, `run_model` and `admissible`, so the parameter vectors are identical to the main
evaluation; only `gamma` differs.  `adm_idx.json` is the list of those 222 draw indices,
extracted from `results/sweep_screened.jsonl`.

Run (resumable by line count of OUT; SHARD_N/SHARD_I split the draw list across cores):

    for i in 0 1 2 3; do
      GAMMA=0.25 OUT=results/leak/leak_g025_sh$i.jsonl SHARD_N=4 SHARD_I=$i \
        python3 leak_sweep.py > results/logs/leak/g025_sh$i.log 2>&1 &
    done
    cat results/leak/leak_g025_sh*.jsonl > results/leak/leak_g025.jsonl

Records have the same shape as `sweep_screened.jsonl` minus the `raw` block.

Result, on the 148 draws admissible for Model 2 at gamma = 0 (and 164 for Model 1):

| gamma | M1 still adm. | M2 still adm. | RT coherence effect | RT duration effect | M2 all six delay slopes > 0 | M1 all six > 0 |
|-------|---------------|---------------|---------------------|--------------------|-----------------------------|----------------|
| 0.00  | 164           | 148           | 139/148             | 138/148            | 76                          | 0              |
| 0.25  | 158           | 146           | 139/148             | 141/148            | 80                          | 0              |
| 0.50  | 151           | 143           | 139/148             | 140/148            | 77                          | 0              |
| 1.00  | 144           | 131           | 141/148             | 141/148            | 75                          | 0              |

The draws that stop being admissible are those with the smallest k*c/theta (median 5.4 s^-1
against 7.9 for those that survive): a leak caps the accumulator at k*c/gamma, and once that
asymptote approaches theta the threshold stops being reached.  That ratio is the criterion for
when a leak can matter at all.
