{-# LANGUAGE NoMonomorphismRestriction #-}

import Data.Char
import Text.Printf
import System.Directory
import System.FilePath
import System.Process

data Criterion = NewCriterion | Exact | BH86 | SW93
               deriving (Eq, Ord, Enum, Bounded)

instance Show Criterion where
  show NewCriterion = "New-Criterion"
  show Exact        = "Exact"
  show BH86         = "BH86"
  show SW93         = "SW93"

data Disk = Exponential | MiyamotoNagai
          deriving (Eq, Ord, Enum, Bounded)

instance Show Disk where
  show Exponential   = "Exponential"
  show MiyamotoNagai = "Miyamoto-Nagai"

data Halo = Logarithmic | NFW | Triaxial
          deriving (Eq, Show, Ord, Enum, Bounded)

data NBodyPerm = NBodyPerm { useQuad   :: Bool,
                             criterion :: Criterion,
                             disk      :: Disk,
                             halo      :: Halo
                           }
               deriving (Eq, Ord)

instance Show NBodyPerm where
  show (NBodyPerm q c d h) = printf "disk=%s_halo=%s_quad=%s_criterion=%s"
                                     (lowshow d) (lowshow h) (lowshow q) (lowshow c)

allTests = [NBodyPerm q c d h | q <- es, c <- es, d <- es, h <- es]
  where es = [minBound .. maxBound]

lower = map toLower
lowshow = lower . show

printTemplate :: String -> NBodyPerm -> String
printTemplate ts np@(NBodyPerm q c d h) = printf ts (lowshow np) (lowshow c) (lowshow q) (lowshow d) (lowshow h)

runTest dir bin name str = do
  let fname = dir </> (name ++ "_results")
  rawSystem bin [ "-s", str, "-o", fname, "-h", histFile ]

makeRunTest dir bin ts np = do let name  = show np
                                   str   = printTemplate ts np
                                   testf = dir </> name <.> ".js"
                               writeFile testf str
                               runTest dir bin name str

generateTests template testDir bin = do
  createDirectoryIfMissing True testDir
  ts <- readFile template
  mapM_ (makeRunTest testDir bin ts) allTests

histFile = "orphan_test_histogram"
templateFile = "orphan_test_template.js"

--main = generateTests "test_set2_template.js" "test_set2" "../bin/milkyway_nbody"
main = generateTests templateFile "test_set1" "../bin/milkyway_nbody"

