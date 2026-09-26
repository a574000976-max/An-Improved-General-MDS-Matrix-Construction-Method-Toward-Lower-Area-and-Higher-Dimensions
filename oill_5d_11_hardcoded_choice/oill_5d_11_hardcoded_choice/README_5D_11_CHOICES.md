# 5D L-power matrices for OILL

This project keeps the original OILL C++ workflow: select a matrix with `#define CHOICE` in `matrix.h`, then compile and run.

## Source mapping from the 5D Sage result

- 4-bit uses `path_4bit`; required powers: `L^-1`, `L`, `L^3`.
- 8/16/32/64-bit use `path_higher`; required powers: `L^-1`, `L^-3`.

## CHOICE mapping

| CHOICE | Bit size | L power | Output file |
|---:|---:|---|---|
| 1 | 4 | `L^-1` | `5x5_n4_L_inv_pow_1_result.txt` |
| 2 | 4 | `L` | `5x5_n4_L_pow_1_result.txt` |
| 3 | 4 | `L^3` | `5x5_n4_L_pow_3_result.txt` |
| 4 | 8 | `L^-1` | `5x5_n8_L_inv_pow_1_result.txt` |
| 5 | 8 | `L^-3` | `5x5_n8_L_inv_pow_3_result.txt` |
| 6 | 16 | `L^-1` | `5x5_n16_L_inv_pow_1_result.txt` |
| 7 | 16 | `L^-3` | `5x5_n16_L_inv_pow_3_result.txt` |
| 8 | 32 | `L^-1` | `5x5_n32_L_inv_pow_1_result.txt` |
| 9 | 32 | `L^-3` | `5x5_n32_L_inv_pow_3_result.txt` |
| 10 | 64 | `L^-1` | `5x5_n64_L_inv_pow_1_result.txt` |
| 11 | 64 | `L^-3` | `5x5_n64_L_inv_pow_3_result.txt` |

## Run

```bash
make clean
make
timeout 5s ./evaluate
```

To switch a matrix, edit:

```cpp
#define CHOICE 1
```

in `matrix.h`, then recompile.
