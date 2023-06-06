library(av)

# Set the directory where your images are located
image_dir <- "D:/PSU Thesis/data/"

regex_folder <- list.files(image_dir,"GaugeCorr_")

for(i in regex_folder){
  # Get a list of image file names in the directory
  image_files <- list.files(paste0(image_dir,i), pattern = "QPE", full.names = TRUE)
 
   # Create a video writer
  output_file <- paste0(i,".mp4")
  
  av::av_encode_video(image_files, paste0(image_dir,output_file), framerate = 6)

}


#utils::browseURL(output_file)