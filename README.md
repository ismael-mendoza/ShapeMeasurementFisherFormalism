# Research-NoiseBias
Repository that contains code for my Weak Lensing summer research project in Prof. Burchat Research Group at Stanford University.
##Interface files.
###generate.py 
Allows the user to create galaxies and save them into a .csv file that can be read by display.py or fitting.py
###display.py
Allows the user to create plots and other visualizations of the statistical results produced by fisher.py and the 
plots produced by draw.py
##Analysis files
###fisher.py
###runfits.py
###readfits.py


##Plotting and Results, 
###draw.py
###readfits.py
###info.py

##Utilities
###defaults.py
###names.py

##ToDo
*galfun organize? it is okay now. 
*change omit -> omit_fits for models.py? 
*readfits.py can be moved to an ipython notebook, as results can all be nicely done only in ipython. 
*generate.py integrate omit_fit? right now is manual passing a dictionary to GParameters, but don't know how necessary this could be. 
*make a demo for producing images and looking at the different parts of GParameters which is the most confusing package. 
*integrate GParameters, fisher and models better?
*images from fisher only in ipython do not save as pdf? 
*can be more efficient if fisher does not have to produce ALL results at once. 

