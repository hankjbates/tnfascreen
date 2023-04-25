macro "VAST Image Organizer" {
	
#@ File (label = "Input Directory", style = "directory") input 
#@ File (label = "Output Directory", style = "directory") output 
#@ String (label = "Output File Name") tableName

	list = getFileList(input);
	daList = newArray(0);
	for (i = 0; i < list.length; i++) {
		subList = newArray(1);
		subList[0] = substring(list[i], 16, 19);
		daList = Array.concat(daList,subList);
	}
	daList = Array.sort(daList);
	Table.create("VAST Images");
	Table.setColumn("Well", daList);
	Table.save(output+"/"+tableName+".csv");
}