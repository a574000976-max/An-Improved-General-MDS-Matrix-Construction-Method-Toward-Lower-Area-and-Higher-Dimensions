# 6D OILL hardcoded choices (4-bit removed)

This version removes all 4-bit matrices and keeps only 8/16/32/64-bit L-power matrices.

## CHOICE mapping

```text
1  = 6x6 n8  L^-2
2  = 6x6 n8  L^-1
3  = 6x6 n8  L

4  = 6x6 n16 L^-2
5  = 6x6 n16 L^-1
6  = 6x6 n16 L

7  = 6x6 n32 L^-2
8  = 6x6 n32 L^-1
9  = 6x6 n32 L

10 = 6x6 n64 L^-2
11 = 6x6 n64 L^-1
12 = 6x6 n64 L
```

Edit `#define CHOICE` in `matrix.h`, then run:

```bash
make clean
make
timeout 5s ./evaluate
```
