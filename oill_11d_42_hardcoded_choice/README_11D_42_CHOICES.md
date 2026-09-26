# 11D hardcoded CHOICE version

This version keeps the original project style: select the matrix by editing `#define CHOICE` in `matrix.h`.
No command-line matrix input and no txt reading are used.

## Build and run

```bash
make clean
make
timeout 5s ./evaluate
```

## CHOICE mapping

```text
 1 = GF(2^16) L^-9
 2 = GF(2^16) L^-6
 3 = GF(2^16) L^-5
 4 = GF(2^16) L^-4
 5 = GF(2^16) L^-3
 6 = GF(2^16) L^-2
 7 = GF(2^16) L^-1
 8 = GF(2^16) L
 9 = GF(2^16) L^2
10 = GF(2^16) L^3
11 = GF(2^16) L^4
12 = GF(2^16) L^6
13 = GF(2^16) L^7
14 = GF(2^16) L^8
15 = GF(2^32) L^-9
16 = GF(2^32) L^-6
17 = GF(2^32) L^-5
18 = GF(2^32) L^-4
19 = GF(2^32) L^-3
20 = GF(2^32) L^-2
21 = GF(2^32) L^-1
22 = GF(2^32) L
23 = GF(2^32) L^2
24 = GF(2^32) L^3
25 = GF(2^32) L^4
26 = GF(2^32) L^6
27 = GF(2^32) L^7
28 = GF(2^32) L^8
29 = GF(2^64) L^-9
30 = GF(2^64) L^-6
31 = GF(2^64) L^-5
32 = GF(2^64) L^-4
33 = GF(2^64) L^-3
34 = GF(2^64) L^-2
35 = GF(2^64) L^-1
36 = GF(2^64) L
37 = GF(2^64) L^2
38 = GF(2^64) L^3
39 = GF(2^64) L^4
40 = GF(2^64) L^6
41 = GF(2^64) L^7
42 = GF(2^64) L^8
```
