module Main where

import Data.Packed.Matrix
import qualified Data.Packed.Vector as V
import Numeric.LinearAlgebra.Algorithms
import Graphics.Plot (imshow)
import Numeric.Container
import Criterion.Main
import Criterion.Config
import Data.Monoid

h :: Double
h = 1.0

-- Initial condition, box with heat
x = 3
y = 3
width = 2
height = 2

laplace :: Int -> Matrix Double
laplace n = let fun (i,j)
                  | i==j = -2.0
                  | i==j-1 || i==j+1 = 1.0
                  | otherwise = 0.0
      in buildMatrix n n fun

initial :: Int -> Matrix Double
initial n = let fun (i,j) = if (x < i) && (i <= (x+width)) && (y < j) && (j <= (y+height))
                         then 1.0
                         else 0.0
      in buildMatrix n n fun

lambda n j = 4/h^2 * (sin (pi*fromIntegral j/(2*(fromIntegral n+1))))^2

eigVectors :: Int -> (Int, Int) -> Double
eigVectors n (i,j) = sqrt(2/(fromIntegral n+1)) * (sin (((fromIntegral i+1)*(fromIntegral j+1)*pi)/(fromIntegral n+1)))

calcSol n =
  let tMat = laplace n
      (dVec,qMat) = eigSH tMat
      --dVec = V.fromList (map (lambda n) [1..n])
      --qMat = buildMatrix n n (eigVectors n)
      gMat = scale (h^2) (initial n)
      gMatM = (trans qMat) `multiply` (gMat `multiply` qMat)
      uFun (i,j) = (gMatM @@> (i,j)) / ((dVec @> i) + (dVec @> j))
      uMatM = buildMatrix n n uFun
  in (qMat `multiply` uMatM) `multiply` (trans qMat)

--myMain = defaultMainWith (defaultConfig { cfgSamples = ljust 2 }) (return ())

--main = myMain [bench "sol" (nf calcSol 1000)]
main = saveMatrix "out.txt" "%f" (calcSol 10)
