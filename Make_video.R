library(av)
library(gtools)
# Set the directory where your images are located
# image_dir <- c("D:\\PSU Thesis\\data\\kuki_GaugeCorr_QPE_01H_2017162_20171132",
#                "D:\\PSU Thesis\\data\\kuki_GaugeCorr_QPE_01H_20192811_201921711",
#                "D:\\PSU Thesis\\data\\kove_GaugeCorr_QPE_01H_2017260_20172100",
#                "D:\\PSU Thesis\\data\\kral_GaugeCorr_QPE_01H_20191142_20191182")


image_dir <- c("G:\\NCFR Thesis\\NCFR_Thesis\\IVT_kove_20172615_20172715\\")
for (j in image_dir){
  regex_folder <- rev(mixedsort(list.files(j,"*")))


  #for(i in regex_folder){
    # Get a list of image file names in the directory
    image_files <- paste0(j,"\\",regex_folder)#list.files(paste0(j,"\\",i), pattern = "QPE", full.names = TRUE)
   
     # Create a video writer
    output_file <- paste0("KOVE_MRMS_Pulse.mov")
    
    av::av_encode_video(image_files, paste0(j,"\\",output_file), framerate = 6)
  
 # }
}

#utils::browseURL(output_file)