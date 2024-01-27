

export class Vector3 {
  x: number
  y: number
  z: number

  constructor(x: number, y: number, z: number) {
    this.x = x;
    this.y = y;
    this.z = z;
  }

  get(i: number): number {
    switch (i) {
      case 0:
        return this.x
      case 1:
        return this.y
      case 2:
        return this.z
      default:
        throw RangeError(`Vector3 only has 3 dimensions, asking for dim ${i}`);
    }
  }

  set(i: number, val: number) {
    switch (i) {
      case 0:
        this.x = val
        break
      case 1:
        this.y = val
        break
      case 2:
        this.z = val
        break
      default:
        throw RangeError(`Vector3 only has 3 dimensions, asking for dim ${i}`);
    }
  }

  copy(): Vector3 {
    return new Vector3(this.x, this.y, this.z)
  }

  minUpdate(o: Vector3) {
    this.x = this.x < o.x ? this.x : o.x
    this.y = this.y < o.y ? this.y : o.y
    this.z = this.z < o.z ? this.z : o.z
  }

  maxUpdate(o: Vector3) {
    this.x = this.x > o.x ? this.x : o.x
    this.y = this.y > o.y ? this.y : o.y
    this.z = this.z > o.z ? this.z : o.z
  }

  addUpdate(o: Vector3) {
    this.x += o.x
    this.y += o.y
    this.z += o.z
  }

  add(o: Vector3): Vector3 {
    return new Vector3(this.x + o.x, this.y + o.y, this.z + o.z)
  }

  subUpdate(o: Vector3) {
    this.x -= o.x
    this.y -= o.y
    this.z -= o.z
  }

  sub(o: Vector3) {
    return new Vector3(this.x - o.x, this.y - o.y, this.z - o.z)
  }

  scaleUpdate(s: number) {
    this.x *= s
    this.y *= s
    this.z *= s
  }

  scale(s: number): Vector3 {
    return new Vector3(this.x * s, this.y * s, this.z * s)
  }

  dot(o: Vector3) {
    return this.x * o.x + this.y * o.y + this.z * o.z
  }

  magSquared() {
    return this.dot(this)
  }

  mag() {
    return Math.sqrt(this.dot(this))
  }

  distanceSquared(o: Vector3) {
    return (this.x - o.x) * (this.x - o.x) + (this.y - o.y) * (this.y - o.y) + (this.z - o.z) * (this.z - o.z)
  }

  distance(o: Vector3) {
    return Math.sqrt(this.distanceSquared(o))
  }
}

const zero = () => new Vector3(0, 0, 0)

export class Particle {
  p: Vector3
  v: Vector3
  r: number
  m: number
  constructor(pos: Vector3, vel: Vector3, r: number, m: number) {
    this.p = pos
    this.v = vel
    this.r = r
    this.m = m
  }

  toString() {
    return `@(${this.p.x}, ${this.p.y}, ${this.p.z})`
  }
}


export function two_bodies(): Particle[] {
  let bodies: Particle[] = []
  bodies.push(new Particle(zero(), zero(), 1.0, 1.0))
  bodies.push(new Particle(new Vector3(1.0, 0.0, 0.0,), new Vector3(0.0, 1.0, 0.0), 1e-4, 1e-20))
  return bodies
}

export function circular_orbits(n: number): Particle[] {
  let particle_buf = []
  particle_buf.push(new Particle(zero(), zero(), 0.00465047, 1.0))

  for (let i = 0; i < n; i++) {
    let d = 0.1 + (i * 5.0 / n)
    let v = Math.sqrt(1.0 / d)
    let theta = Math.random() * 6.28
    let x = d * Math.cos(theta)
    let y = d * Math.sin(theta)
    let vx = -v * Math.sin(theta)
    let vy = v * Math.cos(theta)
    particle_buf.push(new Particle(
      new Vector3(x, y, 0.0),
      new Vector3(vx, vy, 0.0),
      1e-14,
      1e-7,
    ))
  }
  return particle_buf
}

// function calc_accel(i: number, j: number, pi: Particle, pj: Particle, acc: Vector3[]) {
//   let dp = pi.p.sub(pj.p)
//   const dist = dp.mag()
//   let magi = -pj.m / (dist * dist * dist)
//   acc[i].addUpdate(dp.scale(magi))
//   let magj = pi.m / (dist * dist * dist)
//   acc[j].addUpdate(dp.scale(magj))
// }

export function calc_pp_accel(pi: Particle, pj: Particle): Vector3 {
  let dp = pi.p.sub(pj.p);
  let dist = dp.mag();
  let magi = -pj.m / (dist * dist * dist);
  //   println!("magi={}", magi[0]);
  dp.scaleUpdate(magi);
  return dp
}
