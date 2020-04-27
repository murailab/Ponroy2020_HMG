macro "cell ROI time analysis" {

		Dialog.create("Select folders:");
		Dialog.addMessage("input then results");
		Dialog.show();
		input = getDirectory("Select Directory"); 
		results = getDirectory("Select Directory"); 
		// saves results as .csv --> can be read in R
		run("Input/Output...", "jpeg=100 gif=-1 file=.csv copy_column copy_row save_column save_row");

			
 			do {
 				Dialog.create("Radio Buttons");
				items = newArray("ROI", "analyze", "finished");
				Dialog.addRadioButtonGroup("Gee Brain, what are we gonna do tonight?", items, 3, 1, "ROI");
				Dialog.addMessage("The same thing we do every night, Pinky. \n Try to take over the world.");
 				Dialog.show;
				choice=Dialog.getRadioButton;
					
								if (choice=="ROI"){
										Dialog.create("Select SD folders:");
										Dialog.addMessage("select SD folder");
										Dialog.show();
										SD_tif = getDirectory("Select Directory"); 
										setKeyDown("none");
										
										list=getFileList(input);
											for(i=0;i<list.length;i++){
												name=list[i];
												image=input+name;
												if (endsWith(image, ".tif")) {
													roiManager("reset")
													open(image);
													title=File.nameWithoutExtension;
													SD=SD_tif+title+"_SD_.tif";
													open(SD);
													run("Tile");
													waitForUser("add ROIs \n also this is the time to make any necessary corrections \n CHANGES NEED TO BE SAVED MANUALLY \n hold down shift to skip image");
													if (isKeyDown("shift")) {
														run("Close All");
														Dialog.create("thanks for playing");
														Dialog.addMessage("What did you mess up now?");
														Dialog.show();	
														
															
													}
													else {
														roiManager("Save", input+title+"_ROI.zip");
														roiManager("show all with labels");
														run("Flatten");
														saveAs("Tiff", results+title+"_SD_with_ROIs.tif");
														run("Close All");	
													}
													}
												}
											Dialog.create("DONE!!!");
											Dialog.addMessage("I think so, Brain, but there's still a bug stuck in here from last time.");
											Dialog.show();
										}
								if (choice=="analyze"){
									list=getFileList(input);
									for(i=0;i<list.length;i++){
												name=list[i];
												image=input+name;
												roiManager("reset")
												if (endsWith(image, ".tif")) {
													open(image);
													
													title=File.nameWithoutExtension;
													roiManager("Open", input+title+"_ROI.zip");
													selectWindow(name);
														run("Duplicate...", "title=first_pulse duplicate range=40-140");
														first_pulse=getTitle();
															run("Z Project...", "start=1 stop=28 projection=[Average Intensity]");
															first_F0=getTitle();
															imageCalculator("Divide create 32-bit stack", first_pulse, first_F0);
															run("Set Measurements...", "mean redirect=None decimal=3");
															roiManager("Multi Measure");
															saveAs("Results", results+title+"_first_pulse_Results.csv");
															run("Clear Results");
														selectWindow(name);
														run("Duplicate...", "title=second_pulse duplicate range=130-230");
														second_pulse=getTitle();
															run("Z Project...", "start=1 stop=28 projection=[Average Intensity]");
															second_F0=getTitle();
															imageCalculator("Divide create 32-bit stack", second_pulse, second_F0);
															roiManager("Multi Measure");
															saveAs("Results", results+title+"second_pulse_Results.csv");
															run("Clear Results");
														selectWindow(name);
														run("Duplicate...", "title=third_pulse duplicate range=220-320");
														third_pulse=getTitle();
															run("Z Project...", "start=1 stop=28 projection=[Average Intensity]");
															third_F0=getTitle();
															imageCalculator("Divide create 32-bit stack", third_pulse, third_F0);
															roiManager("Multi Measure");
															saveAs("Results", results+title+"_third_pulse_Results.csv");
															run("Clear Results");
												run("Close All");	
													}
												}

											Dialog.create("DONE!!!");
											Dialog.addMessage("Um, I think so, Brainie, but why would anyone want to Pierce Brosnan?");
											Dialog.show();
								
											}
								if (choice=="finished"){
									exit("I think so, Brain, but me and Pippi Longstocking -- I mean, what would the children look like?");
								}
					} while(choice!="finished");
}
