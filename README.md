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
| C++        | Complete         | Worked. Benchmarked. Faster with Clang than gcc.                  |
| Scala      | Complete         | This works, but for some reason it is slow (not as slow as Python) Not taking the time right now to investigate. |
| Java       | Complete         | Benchmarked with both array and OO versions. |
| Python     | Complete         | Benchmarked and HORRIBLY slow.              |
| JavaScript | Complete         | Benchmarked. Not nearly as slow as Python.   |
| GoLang     | Complete         | Benchmarking in progress.           |
| C          | Complete         | Benchmarked. Faster with Clang than gcc.    |
