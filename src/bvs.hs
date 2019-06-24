{-=BasicVariantSelector (BVS): A Haskell-based solution for=-}
{-=selecting variants from a variant calling format=-}
{-=(vcf) file from a mgibed file.=-}
{-=Author: Matthew Mosior=-}
{-=Version: 1.0=-}
{-=Synopsis: This haskell script will take a variant calling=-}
{-=format (vcf) file and a mgi variant=-}
{-=file and will output the filtered vcf=-}
{-=file based on the mgi variant file.=-}

{-Language Extension.-}

{-# LANGUAGE MultiWayIf #-}

{-------------------}

{-Imports-}

import Data.List as DL
import Data.List.Extra as DLE
import Data.List.Split as DLS
import Data.Ord as DO
import Data.Tuple as DT
import System.Console.GetOpt as SCG
import System.Directory as SD
import System.Environment as SE
import System.Exit as SX
import System.IO as SIO
import System.IO.Temp as SIOT
import System.Process as SP

{---------} 

{-Custom CML Option Datatype.-}

data Flag
    = Verbose               -- -v
    | Version               -- -V -?
    | InFormat String       -- -I
    | MgibedFile FilePath   -- -m
    | OutputFile String     -- -o 
    | GzipIn                -- -g
    | GzipOut               -- -G
    | Help                  -- --help
    deriving (Eq,Ord,Show)

{-----------------------------}

{-Option Description function relating to datatype above.-}

--options -> This function will
--describe flags.
options :: [OptDescr Flag]
options =
    [ Option ['v']     ["verbose"]             (NoArg Verbose)                "Output on stderr.",
      Option ['V','?'] ["version"]             (NoArg Version)                "Show version number.",
      Option ['I']     ["InFormat"]            (ReqArg InFormat "IN")         "The format of the input file.",
      Option ['m']     ["MgibedFile"]          (ReqArg MgibedFile "MGIBED")   "The mgibed file from which to select variants.",
      Option ['o']     ["OutputFile"]          (ReqArg OutputFile "OUTFILE")  "The output file.", 
      Option ['g']     ["GzipIn"]              (NoArg GzipIn)                 "Gzipped input file?",
      Option ['G']     ["GzipOut"]             (NoArg GzipOut)                "Gzipped output file?",
      Option []        ["help"]                (NoArg Help)                   "Print this help message."
    ]

{---------------------------------------------------------}

{-Custom bool functions for Flag Datatype.-}

--isOutputFile -> This function will
--test for OutputFile flag.
isOutputFile :: Flag -> Bool
isOutputFile (OutputFile _) = True
isOutputFile _              = False

--isMgibedFile -> This function will
--test for the Mgibedfile flag.
isMgibedFile :: Flag -> Bool
isMgibedFile (MgibedFile _) = True
isMgibedFile _              = False

--isInFormat -> This function will
--test for InFormat flag.
isInFormat :: Flag -> Bool
isInFormat (InFormat _) = True
isInFormat _            = False

--isGzipOut -> This function will
--test for the GzipOut flag.
isGzipOut :: Flag -> Bool
isGzipOut GzipOut = True
isGzipOut _       = False

{------------------------------------------}

{-Custom extraction functions for Flag Datatype.-}

--extractOutputFile -> This function will
--extract the string associated with 
--OutputFile.
extractOutputFile :: Flag -> String
extractOutputFile (OutputFile x) = x

--extractMgibedFile -> This function will
--extract the filepath associated with
--MgibedFile.
extractMgibedFile :: Flag -> FilePath
extractMgibedFile (MgibedFile x) = x

--extractInFormat -> This function will
--extract the string associated with
--InFormat.
extractInFormat :: Flag -> String
extractInFormat (InFormat x) = x

{------------------------------------------------}

{-compilerOpts-related functions-}

--checkInFormat -> This function will
--check the format of IN string.
checkInFormat :: String -> Bool
checkInFormat [] = False
checkInFormat xs = if xs == "vcf" || xs == "tvcf"
                       then True
                       else False

{--------------------}

{-Function to correctly parse the flags.-}

--compilerOpts -> This function will
--parse incoming command line arguments.
compilerOpts :: [String] -> IO ([Flag],String)
compilerOpts argv =
    case getOpt Permute options argv of
        (args,file,[]) ->
            if DL.elem Help args
                then do hPutStrLn stderr (greeting ++ SCG.usageInfo header options)
                        SX.exitWith SX.ExitSuccess
                else if DL.elem Version args
                    then do hPutStrLn stderr (version ++ SCG.usageInfo header options)
                            SX.exitWith SX.ExitSuccess
                    else if (DL.length (DL.filter (isInFormat) args) < 1)
                        then do hPutStrLn stderr (inerror ++ formats ++ SCG.usageInfo header options)
                                SX.exitWith (SX.ExitFailure 1)
                        else if (DL.length (DL.filter (isInFormat) args) > 0) &&
                                (not (checkInFormat (extractInFormat (DL.head (DL.filter (isInFormat) args)))))
                            then do hPutStrLn stderr (inferror ++ formats ++ SCG.usageInfo header options)
                                    SX.exitWith (SX.ExitFailure 1)
                            else if (DL.length (DL.filter (isGzipOut) args) > 0) &&
                                    (DL.length (DL.filter (isOutputFile) args) < 1) 
                                then do hPutStrLn stderr (gziperror ++ SCG.usageInfo header options)
                                        SX.exitWith (ExitFailure 1)
                                else if (DL.length (DL.filter (isMgibedFile) args) < 1)
                                    then do hPutStrLn stderr (mgibedferror ++ SCG.usageInfo header options)
                                            SX.exitWith (ExitFailure 1)
                                    else if DL.length file > 1
                                        then do hPutStrLn stderr (flerror ++ greeting ++ github ++ SCG.usageInfo header options)
                                                SX.exitWith (SX.ExitFailure 1)
                                            else return (DL.nub args,DL.concat file)
        (_,_,errors) -> do
            hPutStrLn stderr (DL.concat errors ++ SCG.usageInfo header options)
            SX.exitWith (SX.ExitFailure 1)
        where
            greeting       = "Basic Variant Selector, Copyright (c) 2019 Matthew Mosior.\n"
            header         = "Usage: vs [-vV?IOmogG] [file]"
            version        = "Basic Variant Selector (VS), Version 1.0.\n"
            github         = "Please see https://github.com/Matthew-Mosior/Variant-to-bam-readcount.\n"
            flerror        = "Incorrect number of input files:  Please provide a single input file.\n" 
            inerror        = "Please provide an input format (-I).\n"
            inferror       = "Input format not recognized.\n" 
            mgibedferror   = "Please provide a single mgibed file from which to select variants (-m)."
            gziperror      = "OutputFile argument (-o) necessary to use GzipOut argument (-G).\n"
            formats        = "Accepted input formats are vcf and tvcf.\nThe accepted output format is mgibed.\n"
           
{----------------------------------------}

{-Vcf or Tvcf pipeline function.-}

--pipeLine -> This function will
--help decide which pipeline the script
--will follow.
pipeLine :: [Flag] -> String
pipeLine []      = []
pipeLine options = if instring == "vcf"
                             then "vcf"
                             else "tvcf"
    where
        --Local definitions.--
        instring  = extractInFormat (DL.head (DL.filter (isInFormat) options))
        ----------------------  

{-------------------------------}

{-Shared General Utility Functions.-}

--lineFeed -> This function will
--read the file in and split on
--whitespace, returning a list
--of lists.
lineFeed :: String -> [[String]]
lineFeed [] = []
lineFeed xs = DL.map DL.words (DL.lines xs)

--mapNotLast -> This function will
--work like the traditional map
--function in Data.List, but not
--map to the last element of a list.
mapNotLast :: (a -> a) -> [a] -> [a]
mapNotLast fn []     = []
mapNotLast fn [x]    = [x]
mapNotLast fn (x:xs) = fn x : mapNotLast fn xs

--tuplifyTwo -> This function will
--turn a list of two elements into
--a two-tuple.
tuplifyTwo :: [a] -> (a,a)
tuplifyTwo [x,y] = (x,y)

--tuplifyThree -> This function will
--turn a list of three elements into
--a triplet.
tuplifyThree :: [a] -> (a,a,a)
tuplifyThree [x,y,z] = (x,y,z)

--singleunnest -> This function will
--unnest a list.
singleunnest :: [a] -> a
singleunnest [a] = a

{-----------------------------------}

{-General utility functions (vcf/tvcf).-}

--onlyInfoBool -> This function will
--return True for only lines that 
--contain the "##INFO" fields.
onlyInfoBool :: [String] -> Bool
onlyInfoBool xs = DL.isInfixOf "##INFO" (DL.head xs) 

--onlyInfoGrabber -> This function will
--grab only lines of file that contain
--the initial header lines.
onlyInfoGrabber :: [[String]] -> [[String]]
onlyInfoGrabber [] = []
onlyInfoGrabber xs = DL.filter (onlyInfoBool) xs

--onlyMetadataVcfBool -> This function will
--return true for only lines that 
--metadata.
onlyMetadataVcfBool :: [String] -> Bool
onlyMetadataVcfBool xs = DL.any (\x -> DL.isPrefixOf "##" x) xs

--onlyMetadataVcfGrabber -> This function will
--grab only lines of the file that contain
--all header lines.
onlyMetadataVcfGrabber :: [[String]] -> [[String]]
onlyMetadataVcfGrabber [] = []
onlyMetadataVcfGrabber xs = DL.filter (onlyMetadataVcfBool) xs

--onlyDataVcfBool -> This function will
--return True for only lines of the
--file that contain tab-delimited
--information.
onlyDataVcfBool :: [String] -> Bool
onlyDataVcfBool xs = not (DL.any (\x -> DL.isPrefixOf "##" x) xs)

--onlyDataVcfGrabber -> This function will
--grab only lines of the file that
--contain tab-delimited information.
onlyDataVcfGrabber :: [[String]] -> [[String]]
onlyDataVcfGrabber [] = []
onlyDataVcfGrabber xs = DL.filter (onlyDataVcfBool) xs

--variantSelectorVcf -> This function will
--filter out variants from a vcf tuple that
--isn't seen in a mgibed file tuple.
variantSelectorVcf :: Int -> Int -> [(String,String)] -> [[String]] -> [[String]]
variantSelectorVcf _          _        _          [] = []
variantSelectorVcf _          _        []         _  = []
variantSelectorVcf chromindex posindex ((a,b):xs) ys = (DL.filter (\c -> (a == (map (\y -> dropPrefix "Chr" y) (map (\x -> dropPrefix "chr" x) c) DL.!! chromindex)) && (b == (c DL.!! posindex))) ys) ++ (variantSelectorVcf chromindex posindex xs ys)
 
{----------------------------------}

{-Printing functions.-}

--tempFileCreation -> This function will
--print the file to stdout using
--readProcess of the unix tool cat.
catFile :: [[String]] -> IO ()
catFile [] = return ()
catFile xs = do
    --Open a temporary file.
    (tempfile,temph) <- SIOT.openTempFile "." "temp.txt"
    --mapNotLast newlines in xs.
    let intercalatedxs = DL.concat (DL.concat (mapNotLast (++ ["\n"]) xs))
    --Add intercalatedxs to temp.txt.
    hPutStrLn temph intercalatedxs
    --Close the temporary file's handle.
    hClose temph
    --Print out the contents of tempfile to the screen using cat unix tool.
    (_,_,_,ph) <- SP.createProcess (SP.proc "cat" [tempfile])
    ec <- SP.waitForProcess ph
    case ec of
        SX.ExitSuccess   -> do SP.readProcess "rm" [tempfile] []
                               return ()
        SX.ExitFailure _ -> do error "Could not cat file."
                               SP.readProcess "rm" [tempfile] []
                               return ()

--noGzipPrintFile -> This function will
--print the file, not gzipped.
noGzipPrintFile :: [Flag] -> [[String]] -> IO ()
noGzipPrintFile [] [] = return ()
noGzipPrintFile [] _  = return ()
noGzipPrintFile _  [] = return ()
noGzipPrintFile opts xs = do
    --Grab just "OUTFILE".
    let outfile = DL.head (DL.filter (isOutputFile) opts)
    --Extract the string from FilterFields.
    let outfilestring = extractOutputFile outfile
    --mapNotLast newlines in xs.
    let intercalatedxs = DL.concat (DL.concat (mapNotLast (++ ["\n"]) xs))
    --Write the output to the user-specified filename.
    SIO.writeFile outfilestring intercalatedxs
               
--gzipPrintFile -> This function will
--will print the file, but gzipped.
gzipPrintFile :: [Flag] -> [[String]] -> IO String
gzipPrintFile [] [] = return []
gzipPrintFile [] _  = return []
gzipPrintFile _  [] = return []
gzipPrintFile opts xs = do
    --Grab just "OUTFILE".
    let outfile = DL.head (DL.filter (isOutputFile) opts)
    --Extract the string from FilterFields.
    let outfilestring = extractOutputFile outfile
    --mapNotLast newlines in xs.
    let intercalatedxs = DL.concat (DL.concat (mapNotLast (++ ["\n"]) xs))
    --Write the output to the user-specified filename.
    SIO.writeFile outfilestring intercalatedxs
    --Gzip outfile.
    SP.readProcess "gzip" [outfilestring] []

{---------------------}

{-BVS Specific Functions.-}
 
--processArgsAndFilesVcfMgibed -> This function will
--walk through each of the command-line
--arguments and files provided by the user.
processArgsAndFilesVcfMgibed :: ([Flag],String) -> IO ()
processArgsAndFilesVcfMgibed ([],[]) = return ()
processArgsAndFilesVcfMgibed (options,inputfile) = do
    --Check to see if inputfile is gzip compressed.
    if DL.elem GzipIn options
        then do --Process vcf file.--
                --Decompress the file.
                SP.readProcess "gunzip" [inputfile] []
                --Read the decompressed vcf file.
                gunzippedfile <- SIO.readFile (DLE.dropSuffix ".gz" inputfile)
                --Apply lineFeed to gunzippedfile.
                let gvcffile = lineFeed gunzippedfile
                --Grab only the metadata lines in gvcffile.
                let gvcfmetadataonly = onlyMetadataVcfGrabber gvcffile
                --Grab only the data lines in gprocessedfile.
                let gvcfdataonly = onlyDataVcfGrabber gvcffile
                --Grab indices of #CHROM of gvcfdataonly.
                let gvcfchromindex = singleunnest (DL.elemIndices "#CHROM" (DL.head gvcfdataonly))
                --Grab indices of POS of gvcfdataonly.
                let gvcfposindex = singleunnest (DL.elemIndices "POS" (DL.head gvcfdataonly)) 
                --Grab the first element of gvcfdataonly.
                let gvcfheadonly = DL.head gvcfdataonly
                --Grab all but the first elements of gvcfdataonly.
                let gtruefinalvcf = DL.tail gvcfdataonly
                ---------------------   
                --Process the mgibed variant file.--
                --Grab just "MGIBED".
                let gmgibedfile = DL.head (DL.filter (isMgibedFile) options)
                --Extract the filepath from MgibedFile.
                let gmgibedfilestring = extractMgibedFile gmgibedfile
                --Read the mgibed file.
                gmgibedfilefinal <- SIO.readFile gmgibedfilestring
                --Apply lineFeed to gmgibedfilefinal.
                let gmgifile = lineFeed gmgibedfilefinal
                --Transpose gmgiprocessedfile.
                let gtransposedmgifile = DL.transpose gmgifile
                --Grab only columns that start with chromosome and start.
                let gfinalmgi = DL.filter (\x -> DL.isInfixOf "Chr" (DL.head x) || DL.isInfixOf "Chromosome" (DL.head x) || DL.isInfixOf "Start" (DL.head x) || DL.isInfixOf "start" (DL.head x)) gtransposedmgifile
                --Grab all but first elements of gfinalmgi.
                let gtruefinalmgi = DL.map (DL.tail) gfinalmgi
                --Combine head and last of gtruefinalmgi into a list of tuples.
                let gtupledmgi = DL.zip (DL.head gtruefinalmgi) (DL.last gtruefinalmgi)                
                --Remove "chr" from beginning of the first element.
                let gnochrtupledmgi = DL.map (\(a,b) -> if
                                                  | DL.isPrefixOf "chr" a -> (DLE.dropPrefix "chr" a,b)
                                                  | DL.isPrefixOf "Chr" a -> (DLE.dropPrefix "Chr" a,b)
                                                  | otherwise             -> (a,b)) gtupledmgi 
                ------------------------------------
                --Filter out variants from gnochrtupledvcf that aren't seen in gnochrtupledmgi.--
                let gvcfvariants = variantSelectorVcf gvcfchromindex gvcfposindex gnochrtupledmgi gtruefinalvcf
                --Add header back onto gfinalvcfvariants.
                let gheadervcfvariants = gvcfheadonly : gvcfvariants
                --mapNotLast tabs to gheadervcfvariants.
                let gheaderreadyvcfvariants = DL.map (mapNotLast (++ "\t")) gheadervcfvariants
                --Add metadata back onto gheadervcfvariants.
                let gprintvariantselector = gvcfmetadataonly ++ gheaderreadyvcfvariants
                ---------------------------------------------------------------------------------
                --Print the file to stdout (cat) or to a file.
                if DL.length (DL.filter (isOutputFile) options) > 0
                    --Check to see if outfile is to be gzipped.
                    then if DL.elem GzipOut options
                        then do
                            _ <- gzipPrintFile options gprintvariantselector
                            return ()
                        else noGzipPrintFile options gprintvariantselector
                else catFile gprintvariantselector
                ----------------------------------------------
        else do --Process vcf file.--
                --Read in the file.
                readinputfile <- SIO.readFile inputfile
                --Apply lineFeed to readinputfile.
                let vcffile = lineFeed readinputfile
                --Grab only the metadata lines in vcffile.
                let vcfmetadataonly = onlyMetadataVcfGrabber vcffile
                --Grab only the data lines in vcffile.
                let vcfdataonly = onlyDataVcfGrabber vcffile
                --Grab indices of #CHROM of vcfdataonly.
                let vcfchromindex = singleunnest (DL.elemIndices "#CHROM" (DL.head vcfdataonly))
                --Grab indices of POS of vcfdataonly.
                let vcfposindex = singleunnest (DL.elemIndices "POS" (DL.head vcfdataonly))
                --Grab the first element of vcfdataonly.
                let vcfheadonly = DL.head vcfdataonly
                --Grab all but the first elements of vcfdataonly.
                let truefinalvcf = DL.tail vcfdataonly
                ---------------------
                --Process the mgibed variant file.--
                --Grab just "MGIBED".
                let mgibedfile = DL.head (DL.filter (isMgibedFile) options)
                --Extract the filepath from MgibedFile.
                let mgibedfilestring = extractMgibedFile mgibedfile
                --Read the mgibed file.
                mgibedfilefinal <- SIO.readFile mgibedfilestring
                --Apply lineFeed to mgibedfilefinal.
                let mgifile = lineFeed mgibedfilefinal
                --Transpose mgiprocessedfile.
                let transposedmgifile = DL.transpose mgifile
                --Grab only columns that start with chromosome and start.
                let finalmgi = DL.filter (\x -> DL.isInfixOf "Chr" (DL.head x) || DL.isInfixOf "Chromosome" (DL.head x) || DL.isInfixOf "Start" (DL.head x) || DL.isInfixOf "start" (DL.head x)) transposedmgifile
                --Grab all but first elements of finalmgi.
                let truefinalmgi = DL.map (DL.tail) finalmgi
                --Combine head and last of truefinalmgi into a list of tuples.
                let tupledmgi = DL.zip (DL.head truefinalmgi) (DL.last truefinalmgi)
                --Remove "chr" from beginning of the first element.
                let nochrtupledmgi = DL.map (\(a,b) -> if
                                                  | DL.isPrefixOf "chr" a -> (DLE.dropPrefix "chr" a,b)
                                                  | DL.isPrefixOf "Chr" a -> (DLE.dropPrefix "Chr" a,b)
                                                  | otherwise             -> (a,b)) tupledmgi
                ------------------------------------
                --Filter out variants from nochrtupledvcf that aren't seen in nochrtupledmgi.--
                let vcfvariants = variantSelectorVcf vcfchromindex vcfposindex nochrtupledmgi truefinalvcf
                --Add header back onto finalvcfvariants.
                let headervcfvariants = vcfheadonly : vcfvariants
                --mapNotLast tabs to headervcfvariants.
                let headerreadyvcfvariants = DL.map (mapNotLast (++ "\t")) headervcfvariants
                --Add metadata back onto headervcfvariants.
                let printvariantselector = vcfmetadataonly ++ headerreadyvcfvariants
                ---------------------------------------------------------------------------------
                --Print the file to stdout (cat) or to a file.
                if DL.length (DL.filter (isOutputFile) options) > 0
                    --Check to see if outfile is to be gzipped.
                    then if DL.elem GzipOut options
                       then do
                            _ <- gzipPrintFile options printvariantselector
                            return ()
                       else noGzipPrintFile options printvariantselector
                else catFile printvariantselector
                ----------------------------------------------

--processArgsAndContentsVcfMgibed -> This function will
--walk through each of the command-line
--arguments and files provided by the user.
processArgsAndContentsVcfMgibed :: ([Flag],String) -> IO ()
processArgsAndContentsVcfMgibed ([],[]) = return ()
processArgsAndContentsVcfMgibed (options,content) = do
    --Process vcf file.--
    --Apply lineFeed to readinputfile.
    let vcffile = lineFeed content
    --Grab only the metadata lines in vcffile.
    let vcfmetadataonly = onlyMetadataVcfGrabber vcffile
    --Grab only the data lines in vcffile.
    let vcfdataonly = onlyDataVcfGrabber vcffile
    --Grab indices of #CHROM of vcfdataonly.
    let vcfchromindex = singleunnest (DL.elemIndices "#CHROM" (DL.head vcfdataonly))
    --Grab indices of POS of vcfdataonly.
    let vcfposindex = singleunnest (DL.elemIndices "POS" (DL.head vcfdataonly))
    --Grab the first element of vcfdataonly.
    let vcfheadonly = DL.head vcfdataonly
    --Grab all but the first elements of vcfdataonly.
    let truefinalvcf = DL.tail vcfdataonly
    ---------------------
    --Process the mgibed variant file.--
    --Grab just "MGIBED".
    let mgibedfile = DL.head (DL.filter (isMgibedFile) options)
    --Extract the filepath from MgibedFile.
    let mgibedfilestring = extractMgibedFile mgibedfile
    --Read the mgibed file.
    mgibedfilefinal <- SIO.readFile mgibedfilestring
    --Apply lineFeed to mgibedfilefinal.
    let mgifile = lineFeed mgibedfilefinal
    --Transpose mgiprocessedfile.
    let transposedmgifile = DL.transpose mgifile
    --Grab only columns that start with chromosome and start.
    let finalmgi = DL.filter (\x -> DL.isInfixOf "Chr" (DL.head x) || DL.isInfixOf "Chromosome" (DL.head x) || DL.isInfixOf "Start" (DL.head x) || DL.isInfixOf "start" (DL.head x)) transposedmgifile
    --Grab all but first elements of finalmgi.
    let truefinalmgi = DL.map (DL.tail) finalmgi
    --Combine head and last of truefinalmgi into a list of tuples.
    let tupledmgi = DL.zip (DL.head truefinalmgi) (DL.last truefinalmgi)
    --Remove "chr" from beginning of the first element.
    let nochrtupledmgi = DL.map (\(a,b) -> if
                                     | DL.isPrefixOf "chr" a -> (DLE.dropPrefix "chr" a,b)
                                     | DL.isPrefixOf "Chr" a -> (DLE.dropPrefix "Chr" a,b)
                                     | otherwise             -> (a,b)) tupledmgi
    ------------------------------------
    --Filter out variants from nochrtupledvcf that aren't seen in nochrtupledmgi.--
    let vcfvariants = variantSelectorVcf vcfchromindex vcfposindex nochrtupledmgi truefinalvcf
    --Add header back onto finalvcfvariants.
    let headervcfvariants = vcfheadonly : vcfvariants
    --mapNotLast tabs to headervcfvariants.
    let headerreadyvcfvariants = DL.map (mapNotLast (++ "\t")) headervcfvariants
    --Add metadata back onto headervcfvariants.
    let printvariantselector = vcfmetadataonly ++ headerreadyvcfvariants
    ---------------------------------------------------------------------------------
    --Print the file to stdout (cat) or to a file.
    if DL.length (DL.filter (isOutputFile) options) > 0
        --Check to see if outfile is to be gzipped.
        then if DL.elem GzipOut options
            then do
                _ <- gzipPrintFile options printvariantselector
                return ()
            else noGzipPrintFile options printvariantselector
    else catFile printvariantselector
    ----------------------------------------------

--processArgsAndFilesTvcfMgibed -> This function will
--walk through each of the command-line
--arguments and files provided by the user.
processArgsAndFilesTvcfMgibed :: ([Flag],String) -> IO ()
processArgsAndFilesTvcfMgibed ([],[]) = return ()
processArgsAndFilesTvcfMgibed (options,inputfile) = do
    --Check to see if inputfile is gzip compressed.
    if DL.elem GzipIn options
        then do --Process tvcf file.--
                --Decompress the file.
                SP.readProcess "gunzip" [inputfile] []
                --Read the decompressed vcf file.
                gunzippedfile <- SIO.readFile (DLE.dropSuffix ".gz" inputfile)
                --Apply lineFeed to gunzippedfile.
                let gtvcffile = lineFeed gunzippedfile
                --Grab only the metadata lines in gtvcffile.
                let gtvcfmetadataonly = onlyMetadataVcfGrabber gtvcffile
                --Grab only the data lines in gvcffile.
                let gtvcfdataonly = onlyDataVcfGrabber gtvcffile
                --Grab indices of #CHROM of gtvcfdataonly.
                let gtvcfchromindex = singleunnest (DL.elemIndices "#CHROM" (DL.head gtvcfdataonly))
                --Grab indices of POS of gtvcfdataonly.
                let gtvcfposindex = singleunnest (DL.elemIndices "POS" (DL.head gtvcfdataonly))
                --Grab the first element of gtvcfdataonly.
                let gtvcfheadonly = DL.head gtvcfdataonly
                --Grab all but the first elements of gtvcfdataonly.
                let gtruefinaltvcf = DL.tail gtvcfdataonly
                ---------------------
                --Process the mgibed variant file.--
                --Grab just "MGIBED".
                let gmgibedfile = DL.head (DL.filter (isMgibedFile) options)
                --Extract the filepath from MgibedFile.
                let gmgibedfilestring = extractMgibedFile gmgibedfile
                --Read the mgibed file.
                gmgibedfilefinal <- SIO.readFile gmgibedfilestring
                --Apply lineFeed to gmgibedfilefinal.
                let gmgifile = lineFeed gmgibedfilefinal
                --Transpose gmgiprocessedfile.
                let gtransposedmgifile = DL.transpose gmgifile
                --Grab only columns that start with chromosome and start.
                let gfinalmgi = DL.filter (\x -> DL.isInfixOf "Chr" (DL.head x) || DL.isInfixOf "Chromosome" (DL.head x) || DL.isInfixOf "Start" (DL.head x) || DL.isInfixOf "start" (DL.head x)) gtransposedmgifile
                --Grab all but first elements of gfinalmgi.
                let gtruefinalmgi = DL.map (DL.tail) gfinalmgi
                --Combine head and last of gtruefinalmgi into a list of tuples.
                let gtupledmgi = DL.zip (DL.head gtruefinalmgi) (DL.last gtruefinalmgi)
                --Remove "chr" from beginning of the first element.
                let gnochrtupledmgi = DL.map (\(a,b) -> if
                                                  | DL.isPrefixOf "chr" a -> (DLE.dropPrefix "chr" a,b)
                                                  | DL.isPrefixOf "Chr" a -> (DLE.dropPrefix "Chr" a,b)
                                                  | otherwise             -> (a,b)) gtupledmgi
                ------------------------------------
                --Filter out variants from gnochrtupledvcf that aren't seen in gnochrtupledmgi.--
                let gtvcfvariants = variantSelectorVcf gtvcfchromindex gtvcfposindex gnochrtupledmgi gtruefinaltvcf
                --Add header back onto gfinalvcfvariants.
                let gheadertvcfvariants = gtvcfheadonly : gtvcfvariants
                --mapNotLast tabs to gheadertvcfvariants.
                let gheaderreadytvcfvariants = DL.map (mapNotLast (++ "\t")) gheadertvcfvariants
                --Add metadata back onto gheadertvcfvariants.
                let gprintvariantselector = gtvcfmetadataonly ++ gheaderreadytvcfvariants
                ---------------------------------------------------------------------------------
                --Print the file to stdout (cat) or to a file.
                if DL.length (DL.filter (isOutputFile) options) > 0
                    --Check to see if outfile is to be gzipped.
                    then if DL.elem GzipOut options
                        then do
                            _ <- gzipPrintFile options gprintvariantselector
                            return ()
                        else noGzipPrintFile options gprintvariantselector
                else catFile gprintvariantselector
                ----------------------------------------------
        else do --Process tvcf file.--
                --Read in the file.
                readinputfile <- SIO.readFile inputfile
                --Apply lineFeed to readinputfile.
                let tvcffile = lineFeed readinputfile
                --Grab only the metadata lines in tvcffile.
                let tvcfmetadataonly = onlyMetadataVcfGrabber tvcffile
                --Grab only the data lines in tvcffile.
                let tvcfdataonly = onlyDataVcfGrabber tvcffile
                --Grab indices of #CHROM of tvcfdataonly.
                let tvcfchromindex = singleunnest (DL.elemIndices "#CHROM" (DL.head tvcfdataonly))
                --Grab indices of POS of tvcfdataonly.
                let tvcfposindex = singleunnest (DL.elemIndices "POS" (DL.head tvcfdataonly))
                --Grab the first element of tvcfdataonly.
                let tvcfheadonly = DL.head tvcfdataonly
                --Grab all but the first elements of tvcfdataonly.
                let truefinaltvcf = DL.tail tvcfdataonly
                ---------------------
                --Process the mgibed variant file.--
                --Grab just "MGIBED".
                let mgibedfile = DL.head (DL.filter (isMgibedFile) options)
                --Extract the filepath from MgibedFile.
                let mgibedfilestring = extractMgibedFile mgibedfile
                --Read the mgibed file.
                mgibedfilefinal <- SIO.readFile mgibedfilestring
                --Apply lineFeed to mgibedfilefinal.
                let mgifile = lineFeed mgibedfilefinal
                --Transpose mgiprocessedfile.
                let transposedmgifile = DL.transpose mgifile
                --Grab only columns that start with chromosome and start.
                let finalmgi = DL.filter (\x -> DL.isInfixOf "Chr" (DL.head x) || DL.isInfixOf "Chromosome" (DL.head x) || DL.isInfixOf "Start" (DL.head x) || DL.isInfixOf "start" (DL.head x)) transposedmgifile
                --Grab all but first elements of finalmgi.
                let truefinalmgi = DL.map (DL.tail) finalmgi
                --Combine head and last of truefinalmgi into a list of tuples.
                let tupledmgi = DL.zip (DL.head truefinalmgi) (DL.last truefinalmgi)
                --Remove "chr" from beginning of the first element.
                let nochrtupledmgi = DL.map (\(a,b) -> if
                                                  | DL.isPrefixOf "chr" a -> (DLE.dropPrefix "chr" a,b)
                                                  | DL.isPrefixOf "Chr" a -> (DLE.dropPrefix "Chr" a,b)
                                                  | otherwise             -> (a,b)) tupledmgi
                ------------------------------------
                --Filter out variants from nochrtupledtvcf that aren't seen in nochrtupledmgi.--
                let tvcfvariants = variantSelectorVcf tvcfchromindex tvcfposindex nochrtupledmgi truefinaltvcf
                --Add header back onto finaltvcfvariants.
                let headertvcfvariants = tvcfheadonly : tvcfvariants
                --mapNotLast tabs to headertvcfvariants.
                let headerreadytvcfvariants = DL.map (mapNotLast (++ "\t")) headertvcfvariants
                --Add metadata back onto headertvcfvariants.
                let printvariantselector = tvcfmetadataonly ++ headerreadytvcfvariants
                ---------------------------------------------------------------------------------
                --Print the file to stdout (cat) or to a file.
                if DL.length (DL.filter (isOutputFile) options) > 0
                    --Check to see if outfile is to be gzipped.
                    then if DL.elem GzipOut options
                        then do
                            _ <- gzipPrintFile options printvariantselector
                            return ()
                        else noGzipPrintFile options printvariantselector
                else catFile printvariantselector
                ----------------------------------------------         

--processArgsAndContentsTvcfMgibed -> This function will
--walk through each of the command-line
--arguments and files provided by the user.
processArgsAndContentsTvcfMgibed :: ([Flag],String) -> IO ()
processArgsAndContentsTvcfMgibed ([],[]) = return ()
processArgsAndContentsTvcfMgibed (options,content) = do
    --Process tvcf file.--
    --Apply lineFeed to readinputfile.
    let tvcffile = lineFeed content
    --Grab only the metadata lines in tvcffile.
    let tvcfmetadataonly = onlyMetadataVcfGrabber tvcffile
    --Grab only the data lines in tvcffile.
    let tvcfdataonly = onlyDataVcfGrabber tvcffile
    --Grab indices of #CHROM of tvcfdataonly.
    let tvcfchromindex = singleunnest (DL.elemIndices "#CHROM" (DL.head tvcfdataonly))
    --Grab indices of POS of tvcfdataonly.
    let tvcfposindex = singleunnest (DL.elemIndices "POS" (DL.head tvcfdataonly))
    --Grab the first element of tvcfdataonly.
    let tvcfheadonly = DL.head tvcfdataonly
    --Grab all but the first elements of tvcfdataonly.
    let truefinaltvcf = DL.tail tvcfdataonly
    ---------------------
    --Process the mgibed variant file.--
    --Grab just "MGIBED".
    let mgibedfile = DL.head (DL.filter (isMgibedFile) options)
    --Extract the filepath from MgibedFile.
    let mgibedfilestring = extractMgibedFile mgibedfile
    --Read the mgibed file.
    mgibedfilefinal <- SIO.readFile mgibedfilestring
    --Apply lineFeed to mgibedfilefinal.
    let mgifile = lineFeed mgibedfilefinal
    --Transpose mgiprocessedfile.
    let transposedmgifile = DL.transpose mgifile
    --Grab only columns that start with chromosome and start.
    let finalmgi = DL.filter (\x -> DL.isInfixOf "Chr" (DL.head x) || DL.isInfixOf "Chromosome" (DL.head x) || DL.isInfixOf "Start" (DL.head x) || DL.isInfixOf "start" (DL.head x)) transposedmgifile
    --Grab all but first elements of finalmgi.
    let truefinalmgi = DL.map (DL.tail) finalmgi
    --Combine head and last of truefinalmgi into a list of tuples.
    let tupledmgi = DL.zip (DL.head truefinalmgi) (DL.last truefinalmgi)
    --Remove "chr" from beginning of the first element.
    let nochrtupledmgi = DL.map (\(a,b) -> if
                                     | DL.isPrefixOf "chr" a -> (DLE.dropPrefix "chr" a,b)
                                     | DL.isPrefixOf "Chr" a -> (DLE.dropPrefix "Chr" a,b)
                                     | otherwise             -> (a,b)) tupledmgi
    ------------------------------------
    --Filter out variants from nochrtupledvcf that aren't seen in nochrtupledmgi.--
    let tvcfvariants = variantSelectorVcf tvcfchromindex tvcfposindex nochrtupledmgi truefinaltvcf
    --Add header back onto finaltvcfvariants.
    let headertvcfvariants = tvcfheadonly : tvcfvariants
    --mapNotLast tabs to headertvcfvariants.
    let headerreadytvcfvariants = DL.map (mapNotLast (++ "\t")) headertvcfvariants
    --Add metadata back onto headertvcfvariants.
    let printvariantselector = tvcfmetadataonly ++ headerreadytvcfvariants
    ---------------------------------------------------------------------------------
    --Print the file to stdout (cat) or to a file.
    if DL.length (DL.filter (isOutputFile) options) > 0
        --Check to see if outfile is to be gzipped.
        then if DL.elem GzipOut options
            then do
                _ <- gzipPrintFile options printvariantselector
                return ()
            else noGzipPrintFile options printvariantselector
    else catFile printvariantselector
    ----------------------------------------------

{-------------------------}

{-Main function.-}

main :: IO ()
main = do
    --Get command line arguments.
    (args,files) <- SE.getArgs >>= compilerOpts
    --See if files is null.
    if null files
        then do --Get stdin.
                contents <- SIO.getContents
                --mgibed parsing pipeline?
                if (pipeLine args) == "vcf"
                    then do --Run args and contents through processArgsAndContentsVcfTvcf.
                            processArgsAndContentsVcfMgibed (args,contents)
                    else do --Run args and contents through processArgsAndContentsTvcfVcf.
                            processArgsAndContentsTvcfMgibed (args,contents) 

        else do --mgibed parsing pipeline?
                if (pipeLine args) == "vcf" 
                    then do --Run args and contents through processArgsAndFilesVcfTvcf.
                            processArgsAndFilesVcfMgibed (args,files)
                    else do --Run args and contents through processArgsAndFilesTvcfVcf.
                            processArgsAndFilesTvcfMgibed (args,files)

{----------------}
