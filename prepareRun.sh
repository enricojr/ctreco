#!/bin/bash

function ExtractMatrixSize () {
    dirRaw=$1
    cd $dirRaw
    oneFile=`ls $dirRaw | tail -n 1`
    size=`cat $oneFile | wc -l`
    echo $size
}

function ExtractBaseFileName () {
    dirRaw=$1
    cd $dirRaw
    oneFile=`ls $dirRaw | head -n 1`
    idx=`expr index $oneFile '[0-9]'`
    # for large numbers this is needed
    if [ $idx -eq 0 ] ; then
	idx=`expr index $oneFile '[0-100000000]'`
    fi
    # back one
    #let idx=$idx-1
    #baseFileName=${oneFile:0:${idx}}
    # prefix
    prefix=${oneFile:0:${idx}-1}
    # Put on the tail
    let restLength=${#oneFile}-${idx}+1
    #baseFileName=${prefix}${baseFileName}"%d"${oneFile:${idx}+1:$restLength}
    baseFileName=${prefix}"%d"${oneFile:${idx}:$restLength}
    # Now get rid of other numbers after "%d"
    ##baseFileName=${baseFileName//[0-9]/}
    echo ${baseFileName}
 }

# check parameters
if [ ${#@} -lt 2 ] ; then
    echo "use: $0  data_directory(absolute_path)  work_directory(absolute_path)"
    exit 1;
fi

originalDataDir=$1
workDataDir=$2

echo ""
echo "[INFO] Creating the work directory $workDataDir"
mkdir -p $workDataDir/raw
echo "[INFO] Working on $originalDataDir"

# Check the number of files
cd $originalDataDir
nFiles=`ls | wc -l`
echo "[INFO] there are $nFiles in the directory. Here the first 10 files ..."
echo "       ---------------------------------------------------------------"
ls | head -n 10
echo "       -------------------------------------------------- scanning ..."
# Find the frame files in the directory
frameFileFound=0

for a in *
do
    numericOnly=0
    # 0 means that it contains alpha characters.
    # Allow spaces, dots, commas, and exponents, i.e. e+06 e-01
    res=`cat $a | tr "\t" "1" | tr " " "1" | tr "." "1" | tr "," "1" | tr "-" "1" | tr "+" "1" | tr "e" "1" | tr "E" "1" | awk ' { if ( $0 ~ /^[0-9\ ]*$/ ) printf "" ; else printf "0" }' | grep 0`
    #echo $res
    # and empty result means a complete numeric file
    if [ ${#res} -eq 0 ] ; then
	cp $a $workDataDir/raw
    fi
done


cd $workDataDir/raw
nFiles=`ls | wc -l`
echo "[INFO] The following seem to be the frame files. Showing the first 10 out of $nFiles."
echo "       ---------------------------------------------------------------"
ls -1v | head -n 10
echo " ... and the last 10 files ..."
ls -1v | tail -n 10
echo "       ---------------------------------------------------------------"
echo "[INFO] These files have been copied to: $workDataDir/raw"
echo "[-->?] Does this look good ? (Hit Enter to continue, Ctrl+C to quit)"
read res
echo "[INFO] D'accord !.  Creating the configuration files ..."


echo "[INFO] Creating the configuration directory $workDataDir/config"
mkdir -p $workDataDir/config
cd $workDataDir/config

#######################################################################
## interpolateBadPixels
echo ""
echo " ***** 1) interpolateBadPixels"
if [ -e interpolateBadPixels.txt ] ; then
    rm -f interpolateBadPixels.txt
fi
touch interpolateBadPixels.txt
if [ ! -e $workDataDir/DB.txt ]; then
    echo "[-->?] WARNING ! This is VERY important. Do not continue without satisfying this requirement"
    echo "       Enter the full path+filename location of your flat field file."
    read dbFile
    cp $dbFile $workDataDir/DB.txt
fi
echo "[INFO] Creating: $workDataDir/config/interpolateBadPixels.txt"
# Extract the base input file name
genericFilename=$(ExtractBaseFileName $workDataDir/raw/)
echo "$workDataDir/raw/$genericFilename" >> interpolateBadPixels.txt
# DB filename
echo $workDataDir/DB.txt >> interpolateBadPixels.txt
# start angle
echo "[-->?] Enter start angle. This doesn't have to be the real angle, just follow the file numbering convention. (and Hit Enter)"
read startAngle
echo $startAngle >> interpolateBadPixels.txt
# finale angle
echo "[-->?] Enter final angle(inclusive). This doesn't have to be the real angle, just follow the file numbering convention. (and Hit Enter)"
read finalAngle
echo $finalAngle >> interpolateBadPixels.txt
# angle step
echo "[-->?] Enter angle step. This doesn't have to be the real angle step, just follow the file numbering convention. (Hit Enter)"
read stepAngle
echo $stepAngle >> interpolateBadPixels.txt
# Matrix size
msize=$(ExtractMatrixSize $workDataDir/raw/)
echo $msize >> interpolateBadPixels.txt
echo $msize >> interpolateBadPixels.txt


#######################################################################
## normalizeToOB
echo ""
echo " ***** 2) normalizeToOB"
if [ -e normalizeToOB.txt ] ; then
    rm -f normalizeToOB.txt
fi
touch normalizeToOB.txt
echo "[INFO] Creating: $workDataDir/config/normalizeToOB.txt"
# input filename
echo "$workDataDir/raw/$genericFilename" >> normalizeToOB.txt
# DB filename
echo $workDataDir/DB.txt >> normalizeToOB.txt
# start angle
echo $startAngle >> normalizeToOB.txt
# finale angle
echo $finalAngle >> normalizeToOB.txt
# angle step
echo $stepAngle >> normalizeToOB.txt
# Matrix size
echo $msize >> normalizeToOB.txt
echo $msize >> normalizeToOB.txt
# normalized frames go to
echo "[INFO] Creating the normalized data directory $workDataDir/normalized"
mkdir -p $workDataDir/normalized
echo $workDataDir/normalized/$genericFilename >> normalizeToOB.txt

#######################################################################
## makeSinogram
echo ""
echo " ***** 3) makeSinogram"
if [ -e makeSinogram.txt ] ; then
    rm -f makeSinogram.txt
fi
touch makeSinogram.txt
echo "[INFO] Creating: $workDataDir/config/makeSinogram.txt"
# input filename
echo "$workDataDir/normalized/$genericFilename" >> makeSinogram.txt
# Matrix size
echo $msize >> makeSinogram.txt
echo $msize >> makeSinogram.txt
# start angle
echo $startAngle >> makeSinogram.txt
# finale angle
echo $finalAngle >> makeSinogram.txt
# angle step
echo $stepAngle >> makeSinogram.txt
# starting slice
echo "[-->?] Enter starting slice (ex: 0) (and Hit Enter)"
read startSlice
echo $startSlice >> makeSinogram.txt
# ending slice
echo "[-->?] Enter ending slice (ex: 255 for 256*256 matrix) (and Hit Enter)"
read endSlice
echo $endSlice >> makeSinogram.txt
# slice direction
echo "[-->?] Enter slice direction: 0=x, 1=y (and Hit Enter)"
read sliceDirection
echo $sliceDirection >> makeSinogram.txt
# format of the output sinogram file names
echo "[INFO] Creating the sinograms directory $workDataDir/sinograms"
mkdir -p $workDataDir/sinograms
echo "$workDataDir/sinograms/sinogram_%d.txt" >> makeSinogram.txt
# data range: minimum pixel row (column)
echo $startSlice >> makeSinogram.txt
# data range: maximum pixel row (column)
echo $endSlice >> makeSinogram.txt
# angle rescale factor
echo "[-->?] Enter the rescale factor.  Each angle will be divided by this factor to obtain the angle in degrees."
read angleRescale
echo $angleRescale >> makeSinogram.txt

#######################################################################
## findAxis
echo ""
echo " ***** 4) findAxis"
if [ -e findAxis.txt ] ; then
    rm -f findAxis.txt
fi
touch findAxis.txt
echo "[INFO] Creating: $workDataDir/config/findAxis.txt"
# input sinograms
echo "$workDataDir/sinograms/sinogram_%d.txt" >> findAxis.txt
# projection size
echo $msize >> findAxis.txt
# format of the reconstruction file name
echo "[INFO] Creating the reconstruction directory $workDataDir/reconstruction/OSEM"
mkdir -p $workDataDir/reconstruction/OSEM
echo "$workDataDir/reconstruction/OSEM/frame_%d.txt" >> findAxis.txt
# OSEM subset size (a subset of 10 yields a resonable computing time to reach the cuts below)
echo "10" >> findAxis.txt
# format of the chi2 file name
echo "$workDataDir/reconstruction/OSEM/chi2_%d.txt" >> findAxis.txt
# chi2 cut on subset iterations
echo "0.01" >> findAxis.txt
# chi2 cut on full set iteration
echo "0.01" >> findAxis.txt
# format of the file names for the chi2 of the rotation axis finding routine
echo "$workDataDir/reconstruction/OSEM/findAxisChi2_%d.txt" >> findAxis.txt
# range around the center where to look for the rotation axis
echo "[-->?] Enter the maximum shift search value. Typically within 5mm (i.e. enter a value of 90)"
read maxShift
echo $maxShift >> findAxis.txt
# format of the sinogram after shift of the rotation axis
echo "$workDataDir/sinograms/sinogram_centeredAxis_%d.txt" >> findAxis.txt
# boolean: search axis --> 1 = do axis search
echo "1" >> findAxis.txt
# manual axis shift
echo "0" >> findAxis.txt
# name of the offset file
echo "$workDataDir/reconstruction/OSEM/offset.txt" >>  findAxis.txt

#######################################################################
## OSEM
echo ""
echo " ***** 5) OSEM"
if [ -e OSEM.txt ] ; then
    rm -f OSEM.txt
fi
touch OSEM.txt
echo "[INFO] Creating: $workDataDir/config/OSEM.txt"
# input sinograms
echo "$workDataDir/sinograms/sinogram_%d.txt" >> OSEM.txt
# projection size
echo $msize >> OSEM.txt
# format of the reconstruction file name
echo "$workDataDir/reconstruction/OSEM/frame_%d.txt" >> OSEM.txt
# OSEM subset size (a subset of 10 yields a resonable computing time to reach the cuts below)
echo "10" >> OSEM.txt
# format of the chi2 file name
echo "$workDataDir/reconstruction/OSEM/chi2_%d.txt" >> OSEM.txt
# chi2 cut on subset iterations
echo "0.01" >> OSEM.txt
# chi2 cut on full set iteration
echo "0.01" >> OSEM.txt
# format of the file names for the chi2 of the rotation axis finding routine
echo "$workDataDir/reconstruction/OSEM/findAxisChi2_%d.txt" >> OSEM.txt
# range around the center where to look for the rotation axis
echo "30" >> OSEM.txt
# format of the sinogram after shift of the rotation axis
echo "$workDataDir/sinograms/sinogram_centeredAxis_%d.txt" >> OSEM.txt
# boolean: search axis --> 1 = do not perform axis search, this is a OSEM reco job
echo "0" >> OSEM.txt
# manual axis shift
echo "-1" >> OSEM.txt
# name of the offset file
echo "$workDataDir/reconstruction/OSEM/offset.txt" >> OSEM.txt

cd ~/analysis/ctreco


# Launch stuff
echo "[INFO] You can run the whole chain following these steps:"
echo "./interpolateBadPixels $workDataDir/config/interpolateBadPixels.txt"
echo "./normalizeToOB $workDataDir/config/normalizeToOB.txt"
echo "./makeSinogram $workDataDir/config/makeSinogram.txt"
echo "./reconstruct_OSEM $workDataDir/config/findAxis.txt $workDataDir/config/makeSinogram.txt"
echo "./reconstruct_OSEM $workDataDir/config/OSEM.txt $workDataDir/config/makeSinogram.txt"

# ask
#echo "[RUN?] Enter to run, Crtl+C to finish"
#read rr1 
#./interpolateBadPixels $workDataDir/config/interpolateBadPixels.txt ; ./normalizeToOB $workDataDir/config/normalizeToOB.txt ; ./makeSinogram $workDataDir/config/makeSinogram.txt ; ./reconstruct_OSEM $workDataDir/config/findAxis.txt $workDataDir/config/makeSinogram.txt ; ./reconstruct_OSEM $workDataDir/config/OSEM.txt $workDataDir/config/makeSinogram.txt
echo "[DONE]"

echo "[-->?] Want to auto run all subsequent steps? (Hit Enter to continue, Ctrl+C to quit)"
read res
echo "[INFO] D'accord !.  Lets go ..."

./interpolateBadPixels $workDataDir/config/interpolateBadPixels.txt; ./normalizeToOB $workDataDir/config/normalizeToOB.txt; ./makeSinogram $workDataDir/config/makeSinogram.txt; ./reconstruct_OSEM $workDataDir/config/findAxis.txt $workDataDir/config/makeSinogram.txt; ./reconstruct_OSEM $workDataDir/config/OSEM.txt $workDataDir/config/makeSinogram.txt
echo "[DONE]"


