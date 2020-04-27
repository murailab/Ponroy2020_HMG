macro "Benny is still trying" {

		Dialog.create("Select folders:");
		Dialog.addMessage(" input \n output  \n MACRO WILL CLOSE ALL OPEN IMAGES");
		Dialog.show();
		input = getDirectory("Select Directory"); 
		output_tif = getDirectory("Select Directory"); 
	//	output_max = getDirectory("Select Directory"); 

		run("Close All");
		
		Dialog.create("filetype: select file ending");
		Dialog.addString("filetype", ".oif");
		Dialog.show();
		filetype=Dialog.getString();

		list=getFileList(input);
			for(i=0;i<list.length;i++){
				name=list[i];
				image=input+name;
				if (endsWith(image, filetype)) {
					run("Bio-Formats Importer", "open=image color_mode=Default split_channels view=Hyperstack stack_order=XYCZT series_1");
					title=File.nameWithoutExtension;
							open_list = getList("image.titles");
									for (j=0; j<open_list.length; j++) {
										saveAs("tiff", output_tif+title+"_ch_"+j+".tif");
									}
				run("Close All");	
					}
				}
}
