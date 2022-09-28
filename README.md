# Multi-Language kD-tree N-body Benchmarks

This repository has small benchmark programs for N-body simulations using 
kD-trees in a variety of languages. This has grown out of an effort to test the
performance of Rust for planetary ring simulations into something of an
updated version of the n-body benchmark that is part of Benchmark Games
(https://benchmarksgame-team.pages.debian.net/benchmarksgame/performance/nbody.html)

## Language Status

| Language   | Basic Status     | Note                                                              |
| ---------- | ---------------- | ----------------------------------------------------------------- |
| Rust       | Complete         | All versions currently use SIMD and no dynamic memory after init. |
| C++        | Compiling        | Needs more verification, but kD-tree looks good.                  |
| Scala      | Stubbed          |                                                                   |
| Java       | Stubbed          |                                                                   |
| Python     | Needs validation | Tests pass & KD-tree looks sensible.                              |
| JavaScript | Needs validation | Tests pass & KD-tree looks sensible.                              |
| GoLang     | Not Started      |                                                                   |
| C          | Needs validation | Tests pass & KD-tree looks sensible.                              |
