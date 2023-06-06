library(av)

# Set the directory where your images are located
image_dir <- "D:/PSU Thesis/scripts/NCFR_Thesis/NCFR_Thesis/"

# Get a list of image file names in the directory
image_files <- list.files(image_dir, pattern = "QPEt", full.names = TRUE)

# Create a video writer
output_file <- "output_video.mp4"

av::av_encode_video(image_files, paste0(image_dir,output_file), framerate = 6)
#utils::browseURL(output_file)