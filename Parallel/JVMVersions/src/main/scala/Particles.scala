object Particles {
  opaque type Particle = Int

  object Particle {
    def apply(index: Int): Particle = index
  }

  extension (part: Particle)
    inline def p(dim: Int)(using system: SSystem) = system.p(part, dim)
    inline def v(dim: Int)(using system: SSystem) = system.v(part, dim)
    inline def r(using system: SSystem) = system.r(part)
    inline def m(using system: SSystem) = system.m(part)

    inline def incP(dim: Int, delta: Double)(using system: SSystem) = system.incP(part, dim, delta)
    inline def incV(dim: Int, delta: Double)(using system: SSystem) = system.incV(part, dim, delta)
}
