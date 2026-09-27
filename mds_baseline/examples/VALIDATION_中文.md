# 本次实际验证结果

测试命令：`python3 -m unittest discover -s tests -v`，10 项测试通过；上传的 11 个系数矩阵全部复现参考 XOR 数。

AddressSanitizer 和 UndefinedBehaviorSanitizer 对 6 个随机可逆矩阵的评估通过。由于当前环境不支持 LeakSanitizer 的进程检查，泄漏检测未运行；详细配置见 sanitizer_validation.json。

示例命令：

```bash
python3 baseline.py run --dims 9 10 11 --samples 4 --iterations 2 --windows 0 --timeout 120 --out examples/validated_run
```

每个维数抽样 4 组 X/Y，同时评估原矩阵和归一化矩阵，合计 8 个候选。每个候选完成 2 次优化迭代。该示例用于验证完整流程，不是充分搜索后的最终基线。

| m | 基础8-bit XOR | n=16 XOR | n=32 XOR | n=64 XOR | 穷举验证的基础子矩阵数 | 超时候选 |
|---:|---:|---:|---:|---:|---:|---:|
| 9 | 1673 | 3346 | 6692 | 13384 | 48619 | 0 |
| 10 | 2110 | 4220 | 8440 | 16880 | 184755 | 0 |
| 11 | 2606 | 5212 | 10424 | 20848 | 705431 | 0 |

每个基础电路和九个扩展电路均通过精确符号验证，证明电路与矩阵对全部输入一致。扩展矩阵的全部系数块也已逐项检查为对应基础块的分块对角副本；结合基础矩阵的全子式验证，可保证其字级 MDS 性质。

代价为当前预算内找到的、可执行且已验证的原地 XOR 实现上界，不是最优 s-XOR 证明。
