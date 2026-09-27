# Results of the Validation Run

Test command: `python3 -m unittest discover -s tests -v`. All 10 tests passed, and the reference XOR counts were reproduced for all 11 uploaded coefficient matrices.

Evaluations of six random invertible matrices passed AddressSanitizer and UndefinedBehaviorSanitizer checks. Leak detection was not run because the current environment does not support the process inspection required by LeakSanitizer. See `sanitizer_validation.json` for the detailed configuration.

Example command:

```bash
python3 baseline.py run --dims 9 10 11 --samples 4 --iterations 2 --windows 0 --timeout 120 --out examples/validated_run
```

For each dimension, four X/Y pairs were sampled. Both the original and normalized matrices were evaluated, giving eight candidates per dimension. Each candidate completed two optimization iterations. This example validates the complete workflow; it is not a final baseline obtained through an extensive search.

| m | Base 8-bit XOR count | n=16 XOR count | n=32 XOR count | n=64 XOR count | Base submatrices exhaustively verified | Timed-out candidates |
|---:|---:|---:|---:|---:|---:|---:|
| 9 | 1673 | 3346 | 6692 | 13384 | 48619 | 0 |
| 10 | 2110 | 4220 | 8440 | 16880 | 184755 | 0 |
| 11 | 2606 | 5212 | 10424 | 20848 | 705431 | 0 |

Each base circuit and all nine lifted circuits passed exact symbolic verification, proving that each circuit agrees with its matrix for every input. Every coefficient block of the lifted matrices was also checked to be a block-diagonal replication of the corresponding base block. Together with exhaustive verification of all nonempty square minors of the base matrices, these checks establish the word-level MDS property of the lifted matrices.

The reported costs are upper bounds obtained from executable, verified in-place XOR implementations found within the specified search budget. They are not proofs of optimal s-XOR costs.
