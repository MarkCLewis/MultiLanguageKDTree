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

## Results

### Quickstat

The first thing I did was refactor out the quickstat so I could easily switch
between sequential and parallel versions. I started with the sequential code.
Then I wrote a parallel version that uses a "counting sort" approach. I then
wrote a version that uses an MPSC channel. Both of the paralle version require
a second array to copy into.

Here is a sample output. The exact times vary dramatically, but the counting
sort parallel version tends to be 2x faster than the sequential version.
However, the MPSC approach was terribly slow. Honestly, I'm surprised the
counting sort approach got speed ups. The partitioning is a very simple
operation with very little work per element. So it is hard to get a speed boost
by trying to break it up.

```
Sequential runtime = 2.250117283
Parallel runtime = 0.536187488
MPSC runtime = 162.455712889
```

### Tree Building

I tried a number of different approaches to building the tree in parallel.
There are the main aspects to the tree building algorithm.

- Find the bounds, mass, and center of mass.
- Do quickstat on the correct dimension.
- Call on the two partitions to build the children.

I tested doing each of these in parallel. Getting the aggregate information
is much faster with chunks because we are calculating a lot of values and
the normal fold spends all the time passing from one iteration to the next.
As discussed above, I considered two parallel options for quickstat. One
was faster than sequential in standalone performacne tests. Lastly, for the
recursive calls, we can use `rayon::join` to do those in parallel.

We tried these in various configurations. The `rayon::join` is the only one
that actually produced faster results. Adding parallel versions of either
of the other two features actually slowed things down.

```
Sequential Runtime = 11.200045432
Par1 Runtime = 195.531032224
Par1 Chunk Runtime = 14.420209785
Par2 Runtime = 32.775782068
Par3 Runtime = 26.056673998
Par3 Chunk Runtime = 2.344795801
Par4 Runtime = 2.217027904
```

What is very odd about these results is that while the parallel quickstat is
faster than the sequential one, the fastest version of the tree building
didn't use parallel quick stat. Indeed, the fastest version is very simple. It
uses the `rayon::join` function to build children in parallel until the
number of branches exceeds the number of virtual cores on the machine.

### Simple Sim

Running with a million particles for 100 steps. Used the sequential tree
construction as well as `build_tree_par4`. The sequential version had a runtime
of 206.4 seconds. Watching the load bar, the steps were obvious because there
was a point where the load when down significantly for each step.

For the parallel version, the run time was 170.2 seconds and the length of
time the load was lower was much smaller and didn't fall so much. So it was
consistently using all the threads the machine had available.

This shows a 15% decrease in total runtime. That is a significant improvement
for just altering the tree building.
