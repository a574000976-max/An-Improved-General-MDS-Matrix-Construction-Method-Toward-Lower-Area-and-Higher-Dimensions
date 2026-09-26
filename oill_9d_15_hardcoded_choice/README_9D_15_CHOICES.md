# 9D L-power hardcoded choices

This project keeps only the 15 L-power matrices needed for the 9-dimensional case.
The matrices are hardcoded in `matrix.cpp`; no external txt input is read.

## Choice mapping

| CHOICE | Matrix | SIZE |
|---:|---|---:|
| 1 | GF2_16_L_inv_pow_2 | 16 |
| 2 | GF2_16_L_inv_pow_1 | 16 |
| 3 | GF2_16_L_pow_1 | 16 |
| 4 | GF2_16_L_pow_2 | 16 |
| 5 | GF2_16_L_pow_3 | 16 |
| 6 | GF2_32_L_inv_pow_2 | 32 |
| 7 | GF2_32_L_inv_pow_1 | 32 |
| 8 | GF2_32_L_pow_1 | 32 |
| 9 | GF2_32_L_pow_2 | 32 |
| 10 | GF2_32_L_pow_3 | 32 |
| 11 | GF2_64_L_inv_pow_2 | 64 |
| 12 | GF2_64_L_inv_pow_1 | 64 |
| 13 | GF2_64_L_pow_1 | 64 |
| 14 | GF2_64_L_pow_2 | 64 |
| 15 | GF2_64_L_pow_3 | 64 |

## Run

1. Edit `matrix.h`:

```cpp
#define CHOICE 1
```

2. Compile and run:

```bash
make clean
make
timeout 5s ./evaluate
```

The result is written to the file defined by `FILENAME`, for example `GF2_16_L_pow_1.txt`.
