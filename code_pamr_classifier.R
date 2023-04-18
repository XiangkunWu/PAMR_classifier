###pamr_classifier
#To load a custom function defined in an R script file:
source("function_pamr_classifier.R")
###pamr_classifier
##Please prepare your TME score file. 
#If you do not have a TME score file yet, please prepare a gene expression profile and use code_calculate_TME_signature_score.R code to calculate TME scores.
pamr_classifier(train_file = "GSVA_score_CRCAFFYcohort.csv",#pamr_classifier file,No modifications required. Please ensure that this file is in the same folder as your preparation file.
                test_file = "GSVA_score_CRCRNAseqcohort.csv",#Please rename your file to follow the format of "GSVA_score_CRCRNAseqcohort.csv" for the file format.
                subtype_file = "CCCRC_CRCAFFYcohort.csv",#pamr_classifier file,No modifications required. Please ensure that this file is in the same folder as your preparation file.
                output_file = "CCCRC_pamr.txt")#Please rename your output file
