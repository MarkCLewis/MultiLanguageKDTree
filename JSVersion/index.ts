import { simple_sim } from "./kdtree";
import { circular_orbits } from "./particle";
import { big_solar, big_solar_with_steps, single_node, two_leaves } from "./test";


let args = process.argv.map(e => e)
// args is node index.js args...
args.shift() // pop 'node'

if (args[1] == 'test') {
  single_node()
  two_leaves()
  big_solar()
  big_solar_with_steps()

  console.log('Tests done')
}
else {
  let n: number = parseInt(args[1])
  let steps: number = parseInt(args[2])

  if (Number.isNaN(n) || Number.isNaN(steps)) {
    console.error(`Usage: ${args[0]} particle_count steps`)
    console.error(args)
    process.exit(1)
  }

  let dt = 1e-3

  simple_sim(circular_orbits(n), dt, steps)
}


