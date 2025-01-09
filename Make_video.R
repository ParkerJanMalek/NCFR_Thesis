library(av)
library(gtools)
# Set the directory where your images are located
# image_dir <- c("D:\\PSU Thesis\\data\\kuki_GaugeCorr_QPE_01H_2017162_20171132",
#                "D:\\PSU Thesis\\data\\kuki_GaugeCorr_QPE_01H_20192811_201921711",
#                "D:\\PSU Thesis\\data\\kove_GaugeCorr_QPE_01H_2017260_20172100",
#                "D:\\PSU Thesis\\data\\kral_GaugeCorr_QPE_01H_20191142_20191182")


#radar_dirs <- c("G:\\NCFR Thesis\\NCFR_Thesis\\combined_kove_20172615_20172715\\")
dirs <- list.dirs("G:/NCFR Thesis/NCFR_Thesis/")
# ivt_dirs <- grep("IVT_*", dirs, value = TRUE)
# radar_dirs <- grep("Radar_*", dirs, value = TRUE)
combined_dirs <- grep("G:/NCFR Thesis/NCFR_Thesis/combined_*", dirs, value = TRUE)
# for (j in ivt_dirs){
#   regex_folder <- rev(mixedsort(list.files(j,"^IVT.*\\.png$")))
# 
#   #for(i in regex_folder){
#     # Get a list of image file names in the directory
#     image_files <- paste0(j,"\\",regex_folder)#list.files(paste0(j,"\\",i), pattern = "QPE", full.names = TRUE)
#    
#      # Create a video writer
#     output_file <- paste0("KOVE_IVT_Pulse.gif")
#     
#     av::av_encode_video(image_files, paste0(j,"\\",output_file), framerate = 6)
#   
#  # }
# }
# 
# for (j in ivt_dirs){
#   regex_folder <- rev(mixedsort(list.files(j,"^850.*\\.png$")))
#   
#   #for(i in regex_folder){
#   # Get a list of image file names in the directory
#   image_files <- paste0(j,"\\",regex_folder)#list.files(paste0(j,"\\",i), pattern = "QPE", full.names = TRUE)
#   
#   # Create a video writer
#   output_file <- paste0("KOVE_Tadv_Pulse.gif")
#   
#   av::av_encode_video(image_files, paste0(j,"\\",output_file), framerate = 6)
#   
#   # }
# }
# 
# for (j in radar_dirs){
#   regex_folder <- mixedsort(list.files(j,"^KBBX.*\\.png$"))
#   
#   #for(i in regex_folder){
#   # Get a list of image file names in the directory
#   image_files <- paste0(j,"\\",regex_folder)#list.files(paste0(j,"\\",i), pattern = "QPE", full.names = TRUE)
#   
#   # Create a video writer
#   output_file <- paste0("KOVE_Radar_Pulse.gif")
#   
#   av::av_encode_video(image_files, paste0(j,"\\",output_file), framerate = 6)
#   
#   # }
# }

#IVT combinedG:\\NCFR Thesis\\NCFR_Thesis\\combined_kove_201721718_201721818
k<- c("G:\\NCFR Thesis\\NCFR_Thesis\\combined_kove_201721718_201721818")#,"G:\\NCFR Thesis\\NCFR_Thesis\\combined_kove_2017220_2017240")
for (j in k){
  print(j)
  try({
  regex_folder <- mixedsort(list.files(j,"^KBBX.*\\presentation.png$"))
  
  #for(i in regex_folder){
  # Get a list of image file names in the directory
  image_files <- paste0(j,"\\",regex_folder)#list.files(paste0(j,"\\",i), pattern = "QPE", full.names = TRUE)
  
  # Create a video writer
  output_file <- paste0("KOVE_combined_Pulse.mov")
  
  av::av_encode_video(image_files, paste0(j,"\\",output_file), framerate = 7)
})
  
  # }
}
# #Tadv combined
# for (j in combined_dirs){
#   regex_folder <- mixedsort(list.files(j,"^850.*\\.png$"))
#   
#   #for(i in regex_folder){
#   # Get a list of image file names in the directory
#   image_files <- paste0(j,"\\",regex_folder)#list.files(paste0(j,"\\",i), pattern = "QPE", full.names = TRUE)
#   
#   # Create a video writer
#   output_file <- paste0("KOVE_combined_Tadv_Pulse.mov")
#   
#   av::av_encode_video(image_files, paste0(j,"\\",output_file), framerate = 6)
#   
#   # }
# }


#utils::browseURL(output_file)