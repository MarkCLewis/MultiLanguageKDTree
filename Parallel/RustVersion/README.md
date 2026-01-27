## Benchmaark Results

The kD-tree build was run on two different machines with both 10 million and
100 million particles. The first machine is
a dual Intel(R) Xeon(R) CPU E5-2680 v3 with 64 GB of RAM. Here are the results
for 100 million particles.
```
Sequential Runtime = 180.993414175
Par1 Runtime = 1453.826940036
Par1 Chunk Runtime = 228.038510697
Par2 Runtime = 330.221600026
Par3 Runtime = 240.49525577
Par3 Chunk Runtime = 37.146947319
Par4 Runtime = 23.656383378
```

The second machine is a dual AMD EPYC 9655 96-Core with 786 GB of RAM. That
machine had the following results for 100 million particles.
```
Sequential Runtime = 68.001499262
Par1 Runtime = 5899.986586855
Par1 Chunk Runtime = 51.514659845
Par2 Runtime = 93.765809181
Par3 Runtime = 868.118120671
Par3 Chunk Runtime = 7.85350348
Par4 Runtime = 10.16924106
```

The only present the 100 million particle results because the two are
basically scaled. Some people might find it intersting that some of the
parallel versions are slower than the sequential version. This is a common
issue because adding parallelism adds overhead. So bad parallelism, where
jobs are too small, can often slow things down.

The most interesting aspect of this is that the optimal approach is different
on the two machines. Par3 Chunk and Par4 are the fastest two on both processors,
but which one actually wins out differs, and by a reasonable amount. Both
are significantly faster than sequential on both machines. Given that the
AMD EPYC is the newer hardware, we will probably go with the Par3 Chunk
approach for real implementations.
