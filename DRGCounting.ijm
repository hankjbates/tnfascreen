// This ImageJ Macro is meant to expedite the process of cell counting and position identification for DRGs of zebrafish 
// Run in Fiji/ImageJ, language is IJ1 Macro
// Henry Bates, Smith Lab, University of Notre Dame
// July 2021

macro "DRG Counting" {

#@ String (visibility="MESSAGE", value="This is an analysis script designed to count fluorescent zebrafish DRG cells", required=false) msg
#@ File (label = "Input Directory", style = "directory") input 
#@ File (label = "Output Directory", style = "directory") output
#@ File (label = "Error Files", style = "directory") manualCount
#@ String (label = "File suffix", value = ".tif") suffix 
#@ String (label = "Output File Name") tableName

	setBatchMode(true);
	resultsArray = newArray(0);
	list = getFileList(input);
	list = Array.sort(list);
	coordArray = newArray(0);
	locationArray = newArray(0);
	duplicateArray = newArray(0);
	wholeArray = newArray(0);
	Table.create("Position of Fluorescent DRGs");
	processFolder(input);

	function processFolder(input) { // Sorts through input folder to find files with .tif suffix
		list = getFileList(input); 
		list = Array.sort(list); 
		for (i = 0; i < list.length; i++) { 
			if (File.isDirectory(input + File.separator + list[i])) {  
				processFolder(input + File.separator + list[i]); 
			}
			if (endsWith(list[i], suffix)) {
				run("Bio-Formats Windowless Importer", "open=["+input+"/"+list[i]+"]");
				setAutoThreshold("Moments dark");
				setOption("BlackBackground", true);
				getThreshold(lower, upper);
				fileName = getInfo("image.filename");
				if (lower < 150) {
					File.rename(input+"/"+fileName, manualCount+"/"+fileName);
				}
				run("Analyze Particles...", "size=1000-30000 add"); //Removes any large and fluorescent noisemakers in the image
				bigROIs = roiManager("count");
				if (bigROIs>0) {
					roiManager("Select", 0);
					getSelectionBounds(x, y, width, height);
					y2 = (y + height + 75);
					makeLine(600, (y - 75), 600, y2, 1200);
					run("Line to Area");
					run("Clear", "slice");
					roiManager("Show All without labels");
					roiManager("Delete");
				}
				setAutoThreshold("Triangle dark"); 
				setOption("BlackBackground", true);
				getThreshold(lower, upper);
				lowthresh = (lower*0.9);
				if (lowthresh > 150) {
					setThreshold(lowthresh, upper); 
					setOption("BlackBackground", true);
				}
				processFile(input, output, list[i]); 
			}
		}
	}
	
	function processFile(input, output, file) { // Opens file and sets conditions to plot ROI and analyze particles  
		getVoxelSize(voxWidth, voxHeight, voxDepth, unit);
		fileName = getInfo("image.filename");
		ROIWidth = newArray(65,25,40); 
		ROIOffset = newArray(20,0,32.5);
		toUnscaled(ROIWidth, ROIOffset);
		IncX=100*voxWidth;
		ImYt=getHeight(); 
		ImXt=getWidth(); 
		xPoints = newArray(0);
		yPoints = newArray(0);
		xLow = newArray(0);
		yLow = newArray(0);
		InitX = 0; 
		LineNum = floor(((ImXt - InitX)/IncX));
		for (i = 0; i <= LineNum; i++) {
			profile = 0;
			SetX=InitX+(IncX*i);
			makeLine(SetX, 0, SetX, (ImYt*1/2)); 
			run("Line Width...", "line="+IncX+"");
			profile = getProfile();
			yMaxima = Array.findMaxima(profile, 15); 
			MaxYArray = Array.trim(yMaxima, 1); 
			if (MaxYArray.length==0) {
				MaxYArray = newArray(1);
				MaxYArray[0] = NaN;
			}
			SetXArray = newArray(1);
			SetXArray[0] = SetX;
			xPoints = Array.concat(xPoints, SetXArray);
			yPoints = Array.concat(yPoints, (MaxYArray[0]+80));
		}
		run("Rotate... ", "angle=180 grid=1 interpolation=Bilinear");
		for (i = 0; i <= LineNum; i++) {
			lowProfile = 0;
			lowSetX=InitX+(IncX*i);
			makeLine(lowSetX, 0, lowSetX, (ImYt*1/2));
			run("Line Width...", "line="+IncX+"");
			lowprofile= getProfile();
			lowYMaxima = Array.findMaxima(lowprofile, 15);
			lowMaxYArray = Array.trim(lowYMaxima, 1);
			if (lowMaxYArray.length==0) {
				lowMaxYArray = newArray(1);
				lowMaxYArray[0] = NaN;
			}
			lowSetXArray = newArray(1);
			lowSetXArray[0] = lowSetX;
			xLow = Array.concat(xLow,lowSetXArray);
			yLow = Array.concat(yLow,(lowMaxYArray[0]+80));
		}
		Fit.doFit("Straight Line", xPoints, yPoints);
		highRSq = Fit.rSquared;
		Fit.doFit("Straight Line", xLow, yLow);
		lowRSq = Fit.rSquared;
		if (highRSq > lowRSq) { // images can come in two different orientations, so we have to differentiate between the two
			yPoints = Array.sort(yPoints);
			run("Rotate... ", "angle=180 grid=1 interpolation=Bilinear");
			if (yPoints[0] > 560) {
				run("Translate...", "x=0 y=-120 interpolation=None");
				for (i = 0; i < yPoints.length; i++) {
					subTransArray = yPoints;
					transArray = newArray(1);
					transArray[0] = (yPoints[i] - 120);
					yPoints = Array.deleteIndex(subTransArray, i);
					subTransArray = yPoints;
					yPoints = Array.concat(subTransArray,transArray);
					Array.sort(yPoints);
				}
			}
			plotROI(input, ROIWidth, ROIOffset, voxWidth, file); 
		}
		if (lowRSq > highRSq) {
			yPoints = yLow;
			xPoints = xLow;
			plotROI(input, ROIWidth, ROIOffset, voxWidth, file);
		}
	}
	
	function plotROI(ROIFile, ROIWidth, ROIOffset, scale, ID) { // Sets ROI for each file	
		for (i = 0; i < yPoints.length; i++) {
			if (yPoints[i] != yPoints[i]) {
				yPoints1 = yPoints;
				yPoints = Array.deleteIndex(yPoints1, i);
				xFit = Array.deleteIndex(xPoints, i);
				nanY = newArray(1);
				Fit.doFit("Straight Line", xFit, yPoints);
				nanY[0] = Fit.f(xPoints[i]);
				yPoints = Array.concat(yPoints,nanY);
				yPoints = Array.sort(yPoints);
			}
		}
		if (yPoints[0] > 560) {
			File.rename(input+"/"+fileName, manualCount+"/"+fileName);
		}
		Array.sort(yPoints);
		ynew = Array.slice(yPoints, 0, (yPoints.length - 4));
		xnew = Array.slice(xPoints, 0, (xPoints.length - 4));
		lastarray = newArray(1);
		Fit.doFit("Straight Line", xnew, ynew);
		lastarray[0] = Fit.f(1200);
		yPoints = Array.concat(ynew, lastarray);
		xPoints = Array.concat(xnew, 1200);
		yTemp = Array.deleteIndex(yPoints, 0);
		xTemp = Array.deleteIndex(xPoints, 0);
		firstpt = newArray(1);
		Fit.doFit("Straight Line", xTemp, yTemp);
		firstpt[0] = Fit.f(0);
		firstx = newArray(1);
		firstx[0] = 0;
		yPoints = Array.concat(yTemp, firstpt);
		xPoints = Array.concat(xTemp, firstx);
		Fit.doFit("Straight Line", xPoints, yPoints);
		highRSq = Fit.rSquared;
		if (highRSq < 0.9) {
			File.rename(input+"/"+fileName, manualCount+"/"+fileName);
		} 
		
		run("Line Width...", "line="+60+"");
		makeSelection("polyline", xPoints, yPoints);
		run("Line to Area");
		ypts = yPoints;
		ypts = Array.sort(ypts);
		if (abs(ypts[ypts.length-1] - ypts[0]) > 350) {
			File.rename(input+"/"+fileName, manualCount+"/"+fileName);
		}
		selectvalue = getValue("selection.size");
		if (selectvalue == 0) {
			File.rename(input+"/"+fileName, manualCount+"/"+fileName);
		}
		else {
			setBackgroundColor(0, 0, 0);
			run("Clear Outside");
			setAutoThreshold("Triangle dark");
			getThreshold(lower, upper);
			highthresh = upper;
			newlow = lower*1.1;
			if (lower < 115) {
				setThreshold(newlow, highthresh);
				setOption("BlackBackground", true);	
			}
			run("Analyze Particles...", "size=20-300 circularity=0.10-1.00 add"); 
			roiManager("deselect"); 
			numResults = roiManager("count"); 
			if (numResults > 0) {
				for (i = 0; i < numResults; i++) {
					roiManager("select", i);
					Roi.getCoordinates(x, y); 
					xFirst = newArray(1);
					xFirst[0] = x[0];
					coordArray = Array.concat(coordArray,xFirst);
					coordArray = Array.sort(coordArray);
					xFirst = newArray(1);
				}
				for (i = 1; i < numResults; i++) {
					if ((abs(coordArray[i] - coordArray[i-1])) < 75) {
						deletion = newArray(1);
						deletion[0] = coordArray[i-1];
						duplicateArray = Array.concat(duplicateArray,deletion);
					}	
				}
				for (i = 0; i < duplicateArray.length; i++) {
					coordArray = Array.deleteValue(coordArray, duplicateArray[i]);
				}
				if (coordArray.length <= 1) {
					File.rename(input+"/"+fileName, manualCount+"/"+fileName);
				}
				if (coordArray.length > 0) {
					coordArray = Array.deleteValue(coordArray, 0);
					roundArray = newArray(1);
					for (i = 0; i < coordArray.length; i++) {
						roundArray[0] = Math.ceil(coordArray[i]/(1200/8)) * (643.4/8);
						wholeArray = Array.concat(wholeArray,roundArray);
					}
					wholeArray = Array.deleteValue(wholeArray, 0);
					for (i = 1; i < wholeArray.length; i++) {
						if (wholeArray[i] == wholeArray[i-1]) {
							wholeArray = Array.deleteIndex(wholeArray, i);
						}
					}
		 			if (wholeArray.length < 1) {
						File.rename(input+"/"+fileName, manualCount+"/"+fileName);
					}
				}
				if (fileName.contains("Capture")) {
					newName = substring(fileName, 0, 24);
					fileName = newName;
				}
			Table.setColumn(""+fileName, wholeArray);
			roiManager("deselect");
			roiManager("delete");
			}
			else {
				File.rename(input+"/"+fileName, manualCount+"/"+fileName);
			}
		}
	} 
	numberImages = newArray(1);
	numberImages[0] = nImages;
	Table.setColumn("Number of Images", numberImages);
	Table.save(output+"/"+tableName+"results.csv"); 
}