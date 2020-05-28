#script that creates the sysdata file
#Masterfile <- readr::read_delim("~/MasterThesis/Annotation_data_sets/Own_AnnotationFile/Annotation_files/Methylation_Masterfile_hg19_270220.txt",
#                            "\t", escape_double = FALSE, comment = "#", trim_ws = TRUE)

#names(Masterfile)[1] = "Markername"
#test = Masterfile$chr
#test = sub("...", "", test)
#head(test)
#test[test == "X"] = 23
#test[test == "Y"] = 24
#test = as.numeric(test)
#summary(test)
#Masterfile$CHR = test
#Masterfile = Masterfile[!is.na(Masterfile$CHR),]
#usethis::use_data(Masterfile, internal = T, overwrite = T)
