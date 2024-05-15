particle_counts=(100000 1000000)
thread_counts=(2 4 6 8 12 24 48)

rm times.txt
for parts in "${particle_counts[@]}"
do
	for threads in "${thread_counts[@]}"
	do
		export RAYON_NUM_THREADS=$threads
		echo $parts $threads >> times.txt
		for cnt in {1..7}
		do
			{ time ./target/release/rust_kdtree_nbody --steps 10 --number $parts ; } 2>> times.txt
		done
	done
done
