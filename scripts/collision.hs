{-# LANGUAGE OverloadedStrings #-}

import           Control.Applicative
import           Control.DeepSeq
import           Control.Monad
import           Control.Parallel.Strategies
import           Data.Aeson
import           Data.Aeson.Types
import qualified Data.ByteString.Lazy.Char8  as B
import           Data.Maybe
import qualified Data.Packed.Vector          as V
import qualified Data.Vector                 as M
import           Data.Vector.Strategies
import           Numeric.Container
import           System.Environment          (getArgs)
import           System.IO.Unsafe            (unsafeInterleaveIO)

data Particle = Particle
    { updated  :: {-# UNPACK #-} !Double
    , charge   :: {-# UNPACK #-} !Int
    , element  :: {-# UNPACK #-} !Int
    , position :: !(V.Vector Double)
    , velocity :: !(V.Vector Double)
    } deriving (Show)

pattern :: String
pattern = "SavePoint*.txt"

instance FromJSON Particle where
    parseJSON (Object v) = Particle
                       <$> v .: "Updated"
                       <*> v .: "Charge"
                       <*> v .: "Type"
                       <*> v .: "Position"
                       <*> v .: "Velocity"
    parseJSON _ = mzero

particleRadii :: M.Vector Double
particleRadii = M.fromList [1.09e-10, 1.30e-10, 1.83e-10, 2.46e-10]

intersecting :: Particle -> Particle -> Bool
intersecting p1 p2 = norm2 diff < (r1 + r2)
    where r1   = particleRadii M.! element p1
          r2   = particleRadii M.! element p2
          diff = position p1 `sub` position p2

getParticles :: FilePath -> IO (M.Vector Particle)
getParticles path = do
    json <- decode `fmap` B.readFile path
    return $! case json of
                   Just a -> M.fromList $ concat $ maybeToList
                           $ parseMaybe (.: "Particles") a
                   _      -> M.empty

{-# INLINE findIntersecting #-}
findIntersecting :: Particle -> M.Vector Particle -> M.Vector Particle
findIntersecting p v = M.filter (intersecting p) v

{-# INLINE countIntersecting #-}
countIntersecting :: M.Vector Particle -> Int
countIntersecting v = M.sum (counts `using` parVector 16)
    where counts = M.map (M.length . intersectingFirst . flip M.drop v)
                 $ M.enumFromStepN 0 1 (M.length v - 1)

          intersectingFirst v | M.null v  = M.empty
                              | otherwise = findIntersecting p (M.drop 1 v)
                              where p = M.head v

main = do
    files <- mapM (unsafeInterleaveIO . getParticles) =<< getArgs
    let counts = map countIntersecting $ filter (not . M.null) files
    mapM_ print counts
