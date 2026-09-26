# 8D oill hardcoded-choice version

This project keeps the original oill `#define CHOICE` workflow.

## CHOICE mapping

 1 = 8x8 n8 L
 2 = 8x8 n8 L^3
 3 = 8x8 n8 L^-3
 4 = 8x8 n8 L^-6
 5 = 8x8 n16 L
 6 = 8x8 n16 L^2
 7 = 8x8 n16 L^-1
 8 = 8x8 n16 L^-2
 9 = 8x8 n32 L
10 = 8x8 n32 L^2
11 = 8x8 n32 L^-1
12 = 8x8 n32 L^-2
13 = 8x8 n64 L
14 = 8x8 n64 L^2
15 = 8x8 n64 L^-1
16 = 8x8 n64 L^-2

## Usage

Edit `matrix.h`:

```cpp
#define CHOICE 1
```

Then run:

```bash
make clean
make
timeout 5s ./evaluate
```
