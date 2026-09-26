# 10D 24 hardcoded CHOICE version

继续沿用原项目 `#define CHOICE` 方式，不读取 txt，不改命令行。

## 使用

```bash
make clean
make
timeout 5s ./evaluate
```

切换矩阵：打开 `matrix.h`，修改：

```cpp
#define CHOICE 1
```

## CHOICE 对应关系

 1 = GF(2^16) L^-4
 2 = GF(2^16) L^-3
 3 = GF(2^16) L^-2
 4 = GF(2^16) L^-1
 5 = GF(2^16) L^1
 6 = GF(2^16) L^2
 7 = GF(2^16) L^3
 8 = GF(2^16) L^4
 9 = GF(2^32) L^-4
10 = GF(2^32) L^-3
11 = GF(2^32) L^-2
12 = GF(2^32) L^-1
13 = GF(2^32) L^1
14 = GF(2^32) L^2
15 = GF(2^32) L^3
16 = GF(2^32) L^4
17 = GF(2^64) L^-4
18 = GF(2^64) L^-3
19 = GF(2^64) L^-2
20 = GF(2^64) L^-1
21 = GF(2^64) L^1
22 = GF(2^64) L^2
23 = GF(2^64) L^3
24 = GF(2^64) L^4
