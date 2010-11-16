{-# LANGUAGE NoMonomorphismRestriction, BangPatterns, PatternGuards #-}
{-# OPTIONS_GHC -W #-}

-- Compile with -threaded

import Data.List
import System.Process
import System.IO
import System.Exit
import Text.Printf
import Control.Applicative
import Control.Monad
import Control.Concurrent
import Control.Exception
import Control.Concurrent.Spawn

-- mass 1 - 50
-- radius 0.05 - 1.0
-- orbit time 1 - 5
-- sim time 1 - 5

bin = "./milkyway_nbody"
outFile = "runtime_results"

inFile = "arst.js"
histogram = "histogram"

nthreads = 2

--nbodySets = [ 1024, 2048, 3072, 4096, 8192, 10000, 15000 ]
nbodySets = [ 100 ]

simpleArguments = arguments 1234 inFile histogram


data FitParams = FitParams { mass        :: Double,
                             radius      :: Double,
                             reverseTime :: Double,
                             forwardTime :: Double
                           } deriving (Eq, Show)

data Workunit = Workunit { nbody  :: Int,
                           fitParams :: FitParams
                         } deriving (Eq, Show)

high = FitParams { mass        = 50.0,
                   radius      = 1.0,
                   reverseTime = 5.0,
                   forwardTime = 5.0
                 }

low = FitParams { mass        = 1.0,
                  radius      = 0.05,
                  reverseTime = 1.0,
                  forwardTime = 1.0
                }

steps = FitParams { mass        = 5.0,
                    radius      = 0.1,
                    reverseTime = 1.0,
                    forwardTime = 1.0
                  }
fpRange :: (FitParams -> Double) -> [Double]
fpRange f = [f low, f low + f steps .. f high]


mapM' :: Monad m => (a-> m b) -> [a] -> m [b]
mapM' _ []     = return []
mapM' f (x:xs) = do y  <- f x
                    ys <- y `seq` mapM' f xs
                    return (y:ys)

findWorkunit :: Int -> [Workunit]
findWorkunit n = map (Workunit n) fps
  where fps = [ FitParams m r 4.0 4.0 | m <- fpRange mass,
                                        r <- fpRange radius ]

findWorkunits :: [Int] -> [Workunit]
findWorkunits ns = concatMap findWorkunit ns

arguments :: Int -> FilePath -> FilePath -> Workunit -> [String]
arguments seed file histogram wu = [ "-t",
                                     "-f", file,
                                     "-h", histogram,
                                     "-e", show seed,
                                     "-u", show (nbody wu),
                                     "-np", show (length params),
                                     "-p" ] ++ params
  where params = map show [ mass fp, radius fp, reverseTime fp, forwardTime fp ]
        fp = fitParams wu

-- assumes only thing written to stdout is "<run_time> %g </run_time>\n"
-- worst function ever
readTimeTag :: String -> Double
readTimeTag str | Just rest <- stripPrefix begin str = let numLen = length rest - endlen
                                                       in if end `isSuffixOf` rest
                                                            then read $ take numLen rest
                                                            else err
                | otherwise = err
  where
    endlen = length end
    begin  = "<run_time> "
    end    = " </run_time>\n"
    err    = (-1.0)

readResult :: Handle -> IO Double
readResult h = do !x <- hGetContents h
                  return (readTimeTag x)

wuString :: Workunit -> Double -> String
wuString wu val = let fp = fitParams wu
                  in printf "%d, %g, %g, %g, %g, %g\n"
                            (nbody wu)
                            (mass fp)
                            (radius fp)
                            (reverseTime fp)
                            (forwardTime fp)
                            val

runWorkunit :: Chan String -> Workunit -> IO ()
runWorkunit results wu = do
  (_, pout, _, h) <- runInteractiveProcess bin (simpleArguments wu) Nothing Nothing
  rc <- waitForProcess h
  rtime <- readResult pout
  hClose pout
  let !resString = if rc /= ExitSuccess
                     then printf "Failed to run process: %s\n" (show rc)
                     else wuString wu rtime
  hPrintf stdout "Completed workunit in %g: %s\n" rtime (show wu)
  writeChan results resString


runWorkunits results runPool wus = do
  mapM' (spawn . runPool . runWorkunit results) wus

readChanN :: Int -> Chan a -> IO [a]
readChanN n chan = take n <$> getChanContents chan

main = do
  results <- newChan :: IO (Chan String)
  runPool <- pool nthreads
  let wus = findWorkunits nbodySets
      n = length wus
  hPrintf stdout "Running %d workunits\n" n
  runWorkunits results runPool wus
  out <- concat <$> readChanN n results
  putStrLn out
  writeFile outFile out

