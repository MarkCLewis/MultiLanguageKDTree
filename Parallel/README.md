# Parallel Versions

The version in this directory support parallel operations. The first cut
was just to parallelize the force calculations. These are embarassingly
parallel, so any language support for parallel for loops works well.

## Parallel kD-Tree Construction

The real challenge is building the kD-tree in parallel. While this is less of
the total runtime than the force calculations, it takes a significant amount
of time for large simulations and becomes a limiting factor how how much
speed we can gain from parallelism.

The tree construction is O(n) at each level and is recursively splits in two.
Once enough splits have happened, we can have each one use a thread, but we
need to handle the early levels and do the quick-stat function in parallel.
If we can do it all with functions that are available in Rayon it would be
avantageous for research simulations.

The bottom levels could use the Rayon join calls on the recursive calls. That
only needs to happen at one level when there are enough tasks to distribute
well.

The bigger challenge is the top levels. Many of the operations can be easily
done with parallel iterators. The hard one is the quick-stat. Partitioning
requires moving memory around in ways that aren't inherently thread-safe.
This could be done with a mpsc channel, but that is probably more overhead
than the work involved. The alternative is to do something akin to a
counting sort. Break the array into chuncks and in each chunk count how many
elements go before (after is num-before). This can be done in parallel. Then
there is an O(P) calculation for the indices each thread would copy into.
An O(P) number of splits makes he groups that can then be handed off to
do the actual moving. Unfortunately, this doesn't happen in place. So a second
array is needed. That could be tricky for ownership.
