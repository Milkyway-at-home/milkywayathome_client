#!/usr/bin/env runhaskell

import System.IO
import Text.Printf
import Control.Applicative

makeCFile :: String -> String
makeCFile = printf "#include \"kernelsource.h\"\n\nconst char* kernelSrc = %s;\n\n"

main = makeCFile . show <$> readFile "arst.cl" >>= writeFile "kernelsource.c"


