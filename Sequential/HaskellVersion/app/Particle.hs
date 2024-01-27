{-# LANGUAGE StrictData #-}
module Particle where

import qualified System.Random.Stateful as Random

data F64x3 = F64x3 {
  x :: !Double,
  y :: !Double,
  z :: !Double}
  deriving (Show, Eq, Ord)

zero :: F64x3
zero = F64x3 0 0 0

-- Vector ops
(@+) :: F64x3 -> F64x3 -> F64x3
(F64x3 x1 y1 z1) @+ (F64x3 x2 y2 z2) = F64x3 (x1 + x2) (y1 + y2) (z1 + z2)
infixl 6 @+

-- dot product
(@*) :: F64x3 -> F64x3 -> Double
(F64x3 x1 y1 z1) @* (F64x3 x2 y2 z2) = (x1 * x2) + (y1 * y2) + (z1 * z2)
infixl 7 @*

(@-) :: F64x3 -> F64x3 -> F64x3
(F64x3 x1 y1 z1) @- (F64x3 x2 y2 z2) = F64x3 (x1 - x2) (y1 - y2) (z1 - z2)
infixl 6 @-

negate :: F64x3 -> F64x3
negate (F64x3 x2 y2 z2) = F64x3 (-x2) (-y2) (-z2)


data Particle = Particle {
  p :: !F64x3,
  v :: !F64x3,
  r :: !Double,
  m :: !Double
} deriving (Show, Eq, Ord)


twoBodies :: [Particle]
twoBodies = [Particle zero zero 1 1, Particle (F64x3 1 0 0) (F64x3 0 1 0) 1e-4 1e-20]


circularOrbits :: Random.StatefulGen g m => g -> Int -> m [Particle]
circularOrbits g n = (first:) <$> traverse (makeOne g) [0..n-1]
  where
    first = Particle zero zero 0.00465047 1

    makeOne :: Random.StatefulGen g m => g -> Int -> m Particle
    makeOne g i = do
      rand <- Random.uniformRM (0 :: Double, 1 :: Double) g
      return $
        let
          d = 0.1 + (fromIntegral i * 5.0 / fromIntegral n)
          v = sqrt $ 1.0 / d
          theta = rand * 6.28
          x = d * cos theta
          y = d * sin theta
          vx = -v * sin theta
          vy = v * cos theta
        in Particle (F64x3 x y 0) (F64x3 vx vy 0) 1e-14 1e-7




