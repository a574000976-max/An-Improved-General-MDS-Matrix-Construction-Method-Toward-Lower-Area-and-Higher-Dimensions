# 7D hardcoded CHOICE version

这个版本按原始项目方式运行：不读取 txt 输入文件，矩阵已经硬编码在 `matrix.cpp`。

## 使用方式

1. 打开 `matrix.h`，修改：

```cpp
#define CHOICE 1
```

2. 编译并运行：

```bash
make clean
make
timeout 5s ./evaluate
```

如果想跑久一点，把 `5s` 改成 `60s` 或更长；如果不限制时间，直接 `./evaluate`。

## CHOICE 对应关系

| CHOICE | Field | L power | SIZE | output file |
|---:|---|---|---:|---|
| 1 | GF(2^8) | L^-3 | 8 | 7x7_n8_L_inv_pow_3_result.txt |
| 2 | GF(2^8) | L^-2 | 8 | 7x7_n8_L_inv_pow_2_result.txt |
| 3 | GF(2^8) | L^-1 | 8 | 7x7_n8_L_inv_pow_1_result.txt |
| 4 | GF(2^8) | L^2 | 8 | 7x7_n8_L_pow_2_result.txt |
| 5 | GF(2^16) | L^-3 | 16 | 7x7_n16_L_inv_pow_3_result.txt |
| 6 | GF(2^16) | L^-2 | 16 | 7x7_n16_L_inv_pow_2_result.txt |
| 7 | GF(2^16) | L^-1 | 16 | 7x7_n16_L_inv_pow_1_result.txt |
| 8 | GF(2^16) | L^2 | 16 | 7x7_n16_L_pow_2_result.txt |
| 9 | GF(2^32) | L^-3 | 32 | 7x7_n32_L_inv_pow_3_result.txt |
| 10 | GF(2^32) | L^-2 | 32 | 7x7_n32_L_inv_pow_2_result.txt |
| 11 | GF(2^32) | L^-1 | 32 | 7x7_n32_L_inv_pow_1_result.txt |
| 12 | GF(2^32) | L^2 | 32 | 7x7_n32_L_pow_2_result.txt |
| 13 | GF(2^64) | L^-3 | 64 | 7x7_n64_L_inv_pow_3_result.txt |
| 14 | GF(2^64) | L^-2 | 64 | 7x7_n64_L_inv_pow_2_result.txt |
| 15 | GF(2^64) | L^-1 | 64 | 7x7_n64_L_inv_pow_1_result.txt |
| 16 | GF(2^64) | L^2 | 64 | 7x7_n64_L_pow_2_result.txt |

注意：`FILENAME` 是输出结果文件名，不是输入矩阵文件名。
