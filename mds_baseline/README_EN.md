# Cauchy MDS Baselines for Dimensions 9, 10, and 11: Generation, OILL Evaluation, and Verification

This project adapts the user-provided `oill_5d_11_hardcoded_choice.zip` into an evaluator that accepts any invertible square binary matrix and runs for a finite number of iterations. It also implements a complete workflow for constructing and evaluating Cauchy baselines.

## Quick Start

Requirements: Python 3.10 or later and g++ with C++17 support. No third-party Python packages or SageMath installation are required. On Windows, WSL/Linux is recommended. Run the following commands from the project directory:

```bash
python3 -m unittest discover -s tests -v
python3 baseline.py run --dims 9 10 11 --samples 4 --iterations 2 --out results_demo
```

The C++ programs are compiled automatically on the first run. Subsequent runs reuse the compiled programs based on source-file hashes. The output directory must be new or empty to avoid mixing results from different configurations.

Example of a larger experiment:

```bash
python3 baseline.py run --dims 9 10 11 --samples 1000 --iterations 8 --windows 20 --window-max 24 --timeout 120 --seed 20260927 --out results_1000
```

This command illustrates an experimental configuration; it does not imply a guaranteed runtime. For each dimension, it samples X and Y 1,000 times and evaluates both the original and normalized matrices for every sample, yielding 2,000 candidates per dimension. Full OILL simplification and local window resynthesis can be time-consuming. Start with a small experiment to estimate the cost on your machine, then choose a search budget.

`--timeout` limits the optimization time for each candidate. The C++ program first saves a verified Gaussian-elimination implementation, then atomically updates the file whenever it finds an improvement. After a timeout, only a fully saved circuit that passes another round of Python verification is accepted. The log records `timed_out=true`; such a run is not treated as having completed all iterations. Results are reproducible with a fixed seed when all iterations finish. Results truncated by a wall-clock timeout may depend on machine load.

## Workflow

1. Select disjoint sets X and Y in GF(2^8), each containing m distinct elements.
2. Construct A[i,j] = 1 / (X[i] + Y[j]), where addition is XOR.
3. Also construct A'[i,j] = A[i,j] A[0,0] / (A[i,0] A[0,j]). The first row and first column of A' consist entirely of ones. Scaling rows and columns by nonzero elements preserves the MDS property.
4. Expand A or A' into an 8m × 8m binary matrix using the specified binary basis.
5. Run the core OILL algorithm from the user-provided code on the entire binary matrix and export an in-place XOR circuit.
6. Retain the candidate with the lowest gate count found and independently check every nonempty square submatrix of the final matrix.
7. For n = 16, 32, and 64, replace each 8 × 8 coefficient block with a block-diagonal matrix containing 2, 4, or 8 copies of that block, respectively, and replicate the corresponding circuit.
8. Verify the block structure and complete circuit of each lifted matrix exactly, then export CSV files, JSON files, and readable XOR sequences.

Default finite field: f(x) = x^8 + x^7 + x^2 + x + 1, encoded as hexadecimal 0x187.

Default basis: `(1,g,g^2+g,g^3+g^2,...,g^7+g^6)`, matching the eight-bit representation in the uploaded paper. Here, `g` denotes `x mod f(x)`. This project uses finite-field arithmetic and does not require g to generate the multiplicative group.

A polynomial basis can also be used:

```bash
python3 baseline.py run --dims 9 10 11 --samples 4 --basis polynomial --polynomial 0x11b --out results_polynomial_basis
```

Different field representations and bases can change circuit costs. Comparisons in a paper must specify which representations were used. The default experiment covers only the specified representation; it does not search all irreducible polynomials or all bases.

## Correctness and Cost Model

### MDS Property of Cauchy Matrices

When the elements within X and Y are distinct and the two sets are disjoint, the Cauchy determinant formula guarantees that every square submatrix has a nonzero determinant. The program additionally checks every nonempty square submatrix of the final instance using an independent C++ finite-field implementation:

| m | Number of checks: C(2m,m) − 1 |
|---:|---:|
| 9 | 48619 |
| 10 | 184755 |
| 11 | 705431 |

Exhaustive minor verification is performed only for the final selected matrix. For the other candidates, the MDS property follows from the Cauchy construction and nonzero row and column scalings, so exhaustive verification is not repeated.

### The Lifted Matrices Remain MDS

After a bit permutation, every k × k block submatrix of a lifted matrix consists of independent copies of the corresponding eight-bit submatrix. Invertibility of the base submatrix therefore guarantees invertibility of the lifted submatrix. The program checks that every coefficient block has exactly this structure, avoiding hundreds of thousands of additional binary rank computations on matrices as large as 704 × 704.

This construction preserves the MDS property measured in n-bit words. Internally, it contains independent eight-bit channels; no additional diffusion properties beyond the word-level MDS property are claimed.

### Circuit Verification

In the circuit JSON:

- `input_permutation[i]` specifies the initialization `z[i] = x[input_permutation[i]]`. All registers are initialized simultaneously.
- Each `[dst,src]` entry in `gates` specifies `z[dst] ^= z[src]`. Gates are executed in the listed order.
- The output is `y[i] = z[i]`, with no omitted output transformation.
- Each update counts as one two-input XOR gate. Bit permutations have zero XOR cost.

The Python verifier treats each input as an independent formal variable, propagates the complete linear expression in every register, and compares the result with the target matrix row by row. This is equivalent to checking all standard basis inputs and proves circuit equivalence for every input, rather than merely testing random inputs.

If the base circuit contains b XOR gates, the specified implementation replicated for n-bit words uses exactly `(n/8)*b` XOR gates. This is the gate count of a constructed implementation and an upper bound on the corresponding in-place s-XOR cost; it is not a proof of global optimality.

The `paper_xor` column in `comparison.csv` contains results from the user-uploaded paper for reference only. This program does not independently verify the paper's high-dimensional constructions or their costs. The reported percentage is `reduction_percent=(baseline_xor-paper_xor)/baseline_xor*100`; a negative value means that the reference count from the paper is higher. Improvements observed in a small example must not be described as improvements over the optimum of the entire family of Cauchy matrices.

## Output Files

```text
results_demo/
  config.json                  Actual configuration, compiler information, and cost model
  comparison.csv               Comparison across three dimensions and three word sizes
  comparison.md                Readable results table
  m9/                          m10 and m11 have the same structure
    search_log.json            X/Y, seed, cost, runtime, and timeout status for each candidate
    candidate_00000_raw/
      matrix.txt               Complete binary matrix
      optimizer.json           Verified circuit
    best_checkpoint.json       Best base instance saved during the search
    best.json                  Final base matrix, X/Y, basis, complete circuit, and index of lifted files
    best_n16.json              Complete 144 × 144 matrix and circuit
    best_n32.json
    best_n64.json
    circuit_n8.txt             Initialization, all XOR operations, and output rules
    circuit_n16.txt
    circuit_n32.txt
    circuit_n64.txt
    verification.json          Independent verification results
```

`examples/validated_run` contains an example that was actually run and verified. It uses four samples per dimension, two versions per sample, and two optimization iterations per candidate. It demonstrates that the complete workflow runs successfully; it is not a final paper baseline obtained through an extensive search.

To independently reverify the saved results:

```bash
python3 baseline.py verify examples/validated_run/m9/best.json
python3 baseline.py verify examples/validated_run/m10/best.json
python3 baseline.py verify examples/validated_run/m11/best.json
```

## Evaluate Your Own Binary Matrix

Input format: the first line contains the matrix dimension N, followed by N lines, each containing N binary digits. The leftmost digit is column zero. For example:

```text
4
1100
0010
0001
1000
```

Run:

```bash
python3 baseline.py evaluate examples/L4_inverse.txt --out my_evaluation --iterations 20 --windows 20
```

This produces `optimizer.json` and `circuit.txt`. The interface supports any invertible square binary matrix and compiles the evaluator automatically for its size. Singular matrices are explicitly rejected. In the internal Python integer representation, column zero is the least significant bit; in the external text format, it is the leftmost digit. The conversion functions handle this distinction consistently. Do not read external rows directly with `int(row,2)`.

## Relationship to the Uploaded OILL Program

The greedy row and column elimination, local rewriting, and equivalent-subsequence resynthesis algorithms from `strategy.cpp` and `reduce.cpp` are retained. They have not been replaced with the Paar or Boyar–Peralta algorithms.

Changes include:

1. Replacing the 11 matrices hardcoded under `CHOICE` with file input, and adding matrix-size and invertibility checks.
2. Fixing the pivot search in `strgy1`: the row index r must be checked against the row count before accessing `mark[r]`.
3. Fixing `reduce3`, which continued iterating using the old sequence length after deleting elements.
4. Default-initializing XOR operation structures to avoid uninitialized flags.
5. Using a fixed seed and a single thread to eliminate concurrent access to a shared `mt19937`. No claim is made that this produces the same results as the original multithreaded search under an equal time budget.
6. Replacing the infinite loop with a specified number of iterations, and replacing the original window search with a configurable, finite number of random window attempts.
7. Checking in C++ that the row operations transform the target matrix into a permutation matrix before saving a result, then independently verifying the final forward circuit in Python.
8. Exporting the complete initialization permutation and output rules, without relying on scattered `y[i]` comments in the original text output.

`oill/UPSTREAM_SHA256.json` records hashes of the uploaded C++ source files. The original upload did not include a license. This package retains its provenance information and makes no additional claim of authorization or change to the licensing of that code.

## Verification Tests

`tests/test_baseline.py` includes:

- Comparison of all 65536 GF(256) products against an independent implementation using carryless multiplication and polynomial long division.
- Checks of the binary multiplication representation for all 256 coefficients and all 256 inputs per coefficient.
- Comparison with the reference XOR counts for the 11 coefficient matrices in the uploaded files.
- Cross-checks between finite-field minors and binary ranks for all square submatrices at small dimensions.
- Rejection tests for singular inputs, corrupted circuits, incorrect gate counts, and corrupted lifted blocks.
- Circuit checks on random invertible matrices using every possible input.
- Reproducibility checks with a fixed random seed.

## Methodological Reference

For Cauchy constructions, complete linear circuit optimization, and subfield lifting, see Thorsten Kranz, Gregor Leander, Ko Stoffelen, and Friedrich Wiemer, *Shorter Linear Straight-Line Programs for MDS Matrices*, ToSC 2017(4). Definition 8 and Lemma 1 describe the subfield construction, preservation of the MDS property, and replication of XOR costs.

Author's public version: https://ko.stoffelen.nl/papers/tosc2017i4-slpmds.pdf

The OILL optimizer in this package is derived from the code uploaded by the user in this session, not from the optimizer used in that paper.
