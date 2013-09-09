module Main where

import Data.Packed.Matrix
import qualified Data.Packed.Vector as V
import Numeric.LinearAlgebra.Algorithms
import Graphics.Plot (imshow)
import Numeric.Container
import Criterion.Main
import Criterion.Config
import Data.Monoid

import Control.Parallel.MPI.Simple (mpiWorld, commWorld, unitTag, send, recv)

h :: Double
h = 1.0

-- Initial condition, box with heat
x = 3
y = 3
width = 2
height = 2

initial :: Int -> Matrix Double
initial n = let fun (i,j) = if (x < i) && (i <= (x+width)) && (y < j) && (j <= (y+height))
                         then 1.0
                         else 0.0
      in buildMatrix n n fun

lambda n j = 4/h^2 * (sin (pi*fromIntegral j/(2*(fromIntegral n+1))))^2

eigVectors :: Int -> (Int, Int) -> Double
eigVectors n (i,j) = sqrt(2/(fromIntegral n+1)) * (sin (((fromIntegral i+1)*(fromIntegral j+1)*pi)/(fromIntegral n+1)))

syncMatrix :: Quad a -> IO (Matrix a)

performMult :: Quad a -> 

calcSol n =
  let dVec = V.fromList (map (lambda n) [1..n])
      qMat = buildMatrix n n (eigVectors n)
      gMat = scale (h^2) (initial n)
      tmp1 = gMatM `multiply` qMat
      gMatM = (trans qMat) `multiply` tmp1 -- (gMat `multiply` qMat)
      uFun (i,j) = (gMatM @@> (i,j)) / ((dVec @> i) + (dVec @> j))
      uMatM = buildMatrix n n uFun
      tmp2 = qMat `multiply` uMatM
  in tmp2 `multiply` (trans qMat)
  --in (qMat `multiply` uMatM) `multiply` (trans qMat)

main = saveMatrix "out.txt" "%f" (calcSol 10)

--main :: IO ()
--main = mpiWorld $ \size rank ->
--   if size < 2
--      then putStrLn "At least two processes are needed"
--      else case rank of
--         0 -> do (msg, _status) <- recv commWorld 1 unitTag
--                 putStrLn msg
--         _ -> send commWorld 0 unitTag "Hello World"
--         _ -> return ()
