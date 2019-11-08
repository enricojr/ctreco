******************************************
*******
******* Enrico Jr. Schioppa & Michele Doni
******* Nikhef detector R&D, April 2014
******* 
******************************************

Hello, this is John.

******* This folder contains all the software you need to reconstruct tomographic data using a Filtered Backprojection (FBP) algorithm or an Ordered Subsets Expectation Maximization (OSEM) algorithm slice by slice, plus additional functions.

******* The reconstruction itself is written in C++, fftw libraries for Fourier transform calculations are required, visualization of data requires ROOT. C++ is not at all the most efficient way to handle this type of data, but if time is not a particular concern than you might want to use this stuff.

******* Raw data must be 1 ASCII matrix frame per angle. Remember that both algorithms apply to -log(I/I0) data or to signal-to-thickness calibrated data!

******* The rotation axis of the images should be aligned with the matrix rows or columns. If it's not the case, the frames must first be rotated. No function here is yet implemented that finds the optimal rotation angle in these cases. One has to determine it "by hand". Best of all, align the sample when you measure the data!

******* Extension convention:
	- .cc -> C++
	- .hh -> C++ headers
	- .C  -> ROOT
	- none -> bash script

******* WARNING! File names ordered by numbers are usually supposed to have integers (%d)

******* WARNING! Even if sometimes the frame dimensions are requested in terms of nx and ny, still it always has to be nx=ny

******* The libraries/ folder contains the classes that are called by the main functions. These classes are overdefined: they contain more stuff than it's actually needed

******* The demo/ folder contains a running example

******* to compile all the functions (except the root functions, extension .C), run the compile script. You might need to explicitly set it into an executable:
	- chmod +x compile

******* plotFrame.C draws a frame + spectrum + cumulative spectrum. Compile in root and execute:
	[] .L plotFrame.C++
	[] plot(const char *fileName, int nx, int ny, bool logscale=false, double min=0., double max=0., const char *ZaxisTitle = "CT number")

******* makeSinogram reads the raw data and makes slice by slice sinograms. It takes a configuration file as an argument, which has to be formatted the following way:

   file name format of input raw frames
   number of pixels nx
   number of pixels ny
   minimum angle position
   maximum angle position
   angle step
   starting slice
   ending slice
   slice direction: 0=x, 1=y
   format of the output sinogram file names
   data range: minimum pixel row (column)
   data range: maximum pixel row (column)
   angle rescale factor
   correction of quad cross artifact

       - the correction of the quad cross artifact starts at the position indicated in the configuration file. If this value is 0, no correction is performed.

******* reconstruct_OSEM serves to both do the OSEM reconstruction, but also to find the rotation axis shift. It takes as second argument the makeSinogram configuration file and as first argument a new configuration file, that has to be formatted as follows:

   format of the input sinogram file name
   projection size
   format of the reconstruction file name
   OSEM subset size
   format of the chi2 file name
   chi2 cut on subset iterations
   chi2 cut on full set iterations
   format of the file names for the chi2 of the rotation axis finding routine
   range around the center where to look for the rotation axis
   format of the sinogram after shift of the rotation axis
   boolean: search axis
   manual axis shift
   name of the offset file

	- If the boolean variable for search axis is set to true, the offset file is created containing the axis shift for each slice. Otherwise, it runs OSEM.
	- WARNING! After fixing the rotation axis, the new frames have a different size! The new size is printed out at the end of the reconstruction, but one can check it later by simply counting words (using the wc command) in one frame.

******* reconstruct_FBP runs the FBP algorithm on sinograms with already centered axis (see reconstruct_OSEM). It takes as second argument the makeSinogram configuration file and as first argument a new configuration file, that has to be formatted as follows:

   file name format of the (centered axis) sinogram
   projection size
   file name format of the reconstruction

	- WARNING! Projection size must be the one after the rotation axis shift

******* drawImages.C creates images out of frames, be them reconstructions or raw frames from the CT scan. Image format can be chosen by just specifying the file name extension (e.g. .png, .jpeg, .gif, etc.). Compile in ROOT and execute:
	[] .L drawImages.C++
	[] int drawImages(int sliceMin, int sliceMax, int sliceStep, int size, double min, double max, const char *fileNameFormat, const char*imageNameFormat, bool logscale = false)
		or
	[] int drawImagesCT(int sliceMin, int sliceMax, int sliceStep, int size, double min, double max, const char *fileNameFormat, const char*imageNameFormat, bool logscale = false)

	WARNING! Open ROOT by calling root -b
	SUGGESTION! Use %0Nd in the image file name format, in order to have N-C zeros preceding a C cyphers number. For example, with %03d, images will be numbered as 001, 002, 003, ... 999. This is useful if you want to make movies out of these images.

******* To make a movie out of a collection of images, use whatever you prefer. I suggest avconv (see demo for an example).

******* interpolateBadPixels.cc looks for pixels with zero counts and interpolates with the neighboring ones. It is meant for CT scan data before flat field correction, so it loops over the raw frames and the OB frame. It takes a configuration file as argument, which has to be formatted as follows:

   format of the input frames files
   format of the input open beam files
   minimum angle position
   maximum angle position
   angle step
   number of pixels nx
   number of pixels ny

	WARNING! Frames are overwritten with the interpolated ones.

******* normalizeToOB.cc calculates -log(I/I0) from the raw data. It takes a configuration file as argument, which has to be formatted as follows:

   format of the input frames files
   format of the input open beam files
   minimum angle position
   maximum angle position
   angle step
   number of pixels nx
   number of pixels ny
   format of the output frames files

******* sumFramesWithFormat.cc sums all frames at a given angle in a CT scan. It takes a configuration file as argument, which has to be formatted as follows:

   format of the input frames files
   minimum angle position
   maximum angle position
   angle step
   minimum frame number
   maximum frame number
   frame number step
   number of pixels nx
   number of pixels ny
   format of the output frames files
