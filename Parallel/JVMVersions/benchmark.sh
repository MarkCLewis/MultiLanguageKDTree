particle_counts=(100000 1000000)
thread_counts=(1 3 5 7 11 23 47)

rm times.txt
for parts in "${particle_counts[@]}"
do
	for threads in "${thread_counts[@]}"
	do
		echo $parts $threads >> times.txt
		for cnt in {1..7}
		do
			{ time java -Djava.util.concurrent.ForkJoinPool.common.parallelism=$threads -cp target/scala-3.1.3/jvmversions_3-0.1.0-SNAPSHOT.jar JavaMain 10 $parts ; } 2>> times.txt
		done
	done
done
