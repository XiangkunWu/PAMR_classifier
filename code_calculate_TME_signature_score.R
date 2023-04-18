# To load a custom function defined in an R script file:
source("function_calculate_TME_signature_score.R")

# To calculate TME signature scores by calling a function:
#The gene expression profiles processing approach for input can follow the instructions of the GSVA package.
calculate_TME_signature_score("test_data_expression.csv", "TME_signature_score_gsva_testdata.csv")
