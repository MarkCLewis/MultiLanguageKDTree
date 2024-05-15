particle_counts=(1000000)
thread_counts=(48 24 12 8 6 4 2)

for parts in "${particle_counts[@]}"
do
	for threads in "${thread_counts[@]}"
	do
		echo $parts $threads >> times.txt
		for cnt in {1..7}
		do
			{ time julia --threads $threads kdtreeparallel.jl 10 $parts ; } 2>> times.txt
		done
	done
done
