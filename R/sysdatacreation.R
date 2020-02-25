#script that creates the sysdata file
Masterfile <- read_delim("~/MasterThesis/Annotation_data_sets/Own_AnnotationFile/Annotation_files/Methylation_Masterfile_hg19_260319.txt",
                            "\t", escape_double = FALSE, comment = "#", trim_ws = TRUE)

names(Masterfile)[3] = "CHR"
test = Masterfile$CHR
test = sub("...", "", test)
head(test)
test[test == "X"] = 23
test[test == "Y"] = 24
test = as.numeric(test)
summary(test)
Masterfile$CHR = test
usethis::use_data(Masterfile, internal = T, overwrite = T)
