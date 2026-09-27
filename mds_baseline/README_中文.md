# 9、10、11 维 Cauchy MDS 比较基线：生成、OILL 评估与验证

本项目把用户提供的 `oill_5d_11_hardcoded_choice.zip` 改为可以读取任意可逆方形二进制矩阵的有限次评估器，并实现完整的 Cauchy 基线流程。

## 快速运行

需要 Python 3.10 或更新版本，以及支持 C++17 的 g++。没有第三方 Python 依赖，也不需要 SageMath。Windows 建议使用 WSL/Linux；在项目目录运行：

```bash
python3 -m unittest discover -s tests -v
python3 baseline.py run --dims 9 10 11 --samples 4 --iterations 2 --out results_demo
```

首次运行会自动编译 C++ 程序，后续按源文件摘要复用。输出目录必须是新的或为空，以免把不同配置的结果混在一起。

扩大实验规模的示例：

```bash
python3 baseline.py run --dims 9 10 11 --samples 1000 --iterations 8 --windows 20 --window-max 24 --timeout 120 --seed 20260927 --out results_1000
```

这条命令是实验配置示例，并不是运行时间承诺。每个维数进行 1000 次 X、Y 抽样，每次同时评估原始和归一化矩阵，即每个维数 2000 个候选。完整 OILL 简化和局部窗口重综合可能耗时较长。先用小规模试验估算本机开销，再确定预算。

`--timeout` 是每个候选的优化时间上限。C++ 会先保存一个已验证的高斯消元实现，再在改善时原子更新文件。超时后只接收已经完整保存且再次通过 Python 验证的电路；日志标记 `timed_out=true`，不把它当作完成全部迭代。固定种子且迭代全部完成时，结果可复现；按墙钟超时截断的结果可能受机器负载影响。

## 计算流程

1. 在 GF(2^8) 中选择互不相交、各有 m 个不同元素的 X、Y。
2. 生成 A[i,j] = 1 / (X[i] + Y[j])，这里加法为 XOR。
3. 同时生成 A'[i,j] = A[i,j] A[0,0] / (A[i,0] A[0,j])。A' 的第一行和第一列均为 1。非零行列倍乘保持 MDS 性质。
4. 根据指定的二进制基，把 A 或 A' 展开为 8m×8m 二进制矩阵。
5. 对整个二进制矩阵运行用户提供的 OILL 核心算法，输出原地 XOR 电路。
6. 在候选中保留已找到的最小门数，独立检查最终矩阵的全部非空方形子矩阵。
7. 对 n=16、32、64，将每个 8×8 系数块替换为其 2、4、8 个副本的分块对角矩阵，并复制对应电路。
8. 精确验证每个扩展矩阵的分块结构及完整电路，输出 CSV、JSON 和可读 XOR 序列。

默认有限域：f(x)=x^8+x^7+x^2+x+1，十六进制编码 0x187。

默认基：`(1,g,g^2+g,g^3+g^2,...,g^7+g^6)`，与上传论文的八比特表示一致。`g` 在这里表示 `x mod f(x)`；本项目使用域的四则运算，不依赖 g 是乘法群生成元。

也可以使用多项式基：

```bash
python3 baseline.py run --dims 9 10 11 --samples 4 --basis polynomial --polynomial 0x11b --out results_polynomial_basis
```

不同域表示和不同基会改变电路成本，论文比较中必须说明选用了哪些表示。默认试验只覆盖指定的一种表示，没有搜索所有不可约多项式或所有基。

## 正确性与计数口径

### Cauchy 矩阵的 MDS 性质

当 X、Y 内部元素互异且互不相交时，Cauchy 行列式公式保证每个方形子矩阵的行列式非零。程序额外使用独立 C++ 有限域实现逐一验证最终实例的所有非空方形子矩阵：

| m | 检查数量 C(2m,m)-1 |
|---:|---:|
| 9 | 48619 |
| 10 | 184755 |
| 11 | 705431 |

只对最终选中的矩阵穷举全部子式。其余候选依据 Cauchy 构造和非零行列倍乘保证 MDS，不重复穷举所有子式。

### 扩展矩阵仍然 MDS

扩展后任意 k×k 分块子矩阵，在重新排列比特后都是相应八比特子矩阵的若干独立副本。基础子矩阵可逆，就保证扩展后的子矩阵可逆。程序精确检查每个系数块确实具有这一结构，因此无需重新穷举数十万次最高 704×704 的二进制秩运算。

这种构造保持按 n 比特字计算的 MDS 性质；它内部有相互独立的八比特通道，不宣称拥有超出字级 MDS 性质的其他扩散特征。

### 电路验证

电路 JSON 中：

- `input_permutation[i]` 表示初始 `z[i] = x[input_permutation[i]]`，所有寄存器同时初始化。
- `gates` 中每个 `[dst,src]` 表示 `z[dst] ^= z[src]`，按列出顺序执行。
- 输出为 `y[i] = z[i]`，没有省略的输出变换。
- 每条更新计一个二输入 XOR；比特排列不计 XOR。

Python 验证器将每个输入当作独立形式变量，完整传播每个寄存器的线性表达式，并逐行比较目标矩阵。这相当于验证全部单位基输入，能够证明所有输入上的电路等价性，不只是随机测试。

如果基础电路含 b 个 XOR，复制到 n 比特字的指定实现恰好使用 `(n/8)*b` 个 XOR。这是已构造实现的门数，也是相应原地 s-XOR 的一个上界，不是证明全局最优。

`comparison.csv` 中的 `paper_xor` 来自用户上传论文的结果，仅作为参考列；本程序没有重新验证论文的高维构造及其成本。`reduction_percent=(baseline_xor-paper_xor)/baseline_xor*100`，负值表示本文参考计数更高。不能把小规模示例的优势比例表述为相对于整个 Cauchy 矩阵族的最优优势。

## 文件输出

```text
results_demo/
  config.json                  实际配置、编译器信息与计数口径
  comparison.csv               三种维数与三个字长的对照表
  comparison.md                便于阅读的结果表
  m9/                          m10、m11 结构相同
    search_log.json             每个候选的 X/Y、种子、成本、时间、超时状态
    candidate_00000_raw/
      matrix.txt               完整二进制矩阵
      optimizer.json           已验证电路
    best_checkpoint.json       搜索过程中实时保存的最佳基础实例
    best.json                  最终基础矩阵、X/Y、基、完整电路及扩展文件索引
    best_n16.json              完整 144×144 矩阵及电路
    best_n32.json
    best_n64.json
    circuit_n8.txt             初始化、所有 XOR 操作与输出规则
    circuit_n16.txt
    circuit_n32.txt
    circuit_n64.txt
    verification.json          独立验证结果
```

`examples/validated_run` 是本次实际运行并验证的示例。其规模是每维 4 次抽样、每次两个版本、每候选 2 次优化迭代。它用于证明完整流程可运行，不作为已经充分搜索的最终论文基线。

独立重新验证已保存结果：

```bash
python3 baseline.py verify examples/validated_run/m9/best.json
python3 baseline.py verify examples/validated_run/m10/best.json
python3 baseline.py verify examples/validated_run/m11/best.json
```

## 评估自己的二进制矩阵

输入格式：首行是矩阵尺寸 N，随后 N 行，每行 N 个 0/1，左边第一位是第零列。例如：

```text
4
1100
0010
0001
1000
```

运行：

```bash
python3 baseline.py evaluate examples/L4_inverse.txt --out my_evaluation --iterations 20 --windows 20
```

得到 `optimizer.json` 和 `circuit.txt`。这个接口支持任意可逆方形二进制矩阵，按尺寸自动编译。奇异矩阵会明确拒绝。内部 Python 整数矩阵的第零列是最低有效位，外部文本的第零列是最左位；转换函数已统一处理，不要直接用 `int(row,2)` 读取外部行。

## 与上传 OILL 程序的关系

保留了 `strategy.cpp`、`reduce.cpp` 的贪心行列消元、局部重写及等价子序列重综合算法。没有把它换成 Paar 或 Boyar–Peralta 算法。

修改包括：

1. 用文件输入替代 `CHOICE` 中硬编码的 11 个矩阵，增加矩阵尺寸和可逆性检查。
2. 修正 `strgy1` 的 pivot 搜索：必须先检查 r 没有超过行数，再访问 `mark[r]`。
3. 修正 `reduce3` 删除元素后仍按旧序列长度遍历的问题。
4. 默认初始化 XOR 操作结构，避免未初始化标志。
5. 使用固定种子和单线程，消除共享 `mt19937` 的并发访问；不声称与原多线程搜索在同一时间预算下得到相同结果。
6. 将无限循环改为指定迭代次数；原有窗口搜索替换为可配置的有限次随机窗口尝试。
7. 在 C++ 保存结果前检查行操作确实将目标矩阵化为置换矩阵；再由 Python 独立验证最终正向电路。
8. 导出完整初始化排列与输出规则，不依赖原文本中零散的 `y[i]` 注释。

`oill/UPSTREAM_SHA256.json` 记录了上传 C++ 源文件摘要。原上传材料未附许可证，本包保留其来源说明，不额外宣称对该部分拥有授权或改变其许可。

## 验证入口

`tests/test_baseline.py` 包括：

- 独立的无进位乘法与多项式长除，对比全部 65536 个 GF(256) 乘积；
- 全部 256 个系数、每个 256 个输入的二进制乘法表示检查；
- 对照上传文件中 11 个系数矩阵的参考 XOR 数量；
- 小维数全部子式的有限域与二进制秩交叉验证；
- 奇异输入、损坏电路、错误门数和损坏扩展块的拒绝测试；
- 随机可逆矩阵与全部输入的电路检查；
- 固定随机种子的复现检查。

## 方法出处

Cauchy 构造、完整线性电路优化和子域扩展可参考：Thorsten Kranz, Gregor Leander, Ko Stoffelen, Friedrich Wiemer, *Shorter Linear Straight-Line Programs for MDS Matrices*, ToSC 2017(4)。其中 Definition 8 和 Lemma 1 给出了子域构造及保持 MDS、复制 XOR 代价的性质。

作者公开版本：https://ko.stoffelen.nl/papers/tosc2017i4-slpmds.pdf

本包的 OILL 优化实现来源于本次用户上传的代码，而非上述论文的优化器。
