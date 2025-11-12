
#' Encode DNA Sequences into One-Hot Format
#'
#' This function converts a vector of DNA strings into a one-hot encoded
#' data frame. Each position in the DNA sequence becomes a set of columns
#' representing the nucleotide (A, T, C, G). This is an internal helper
#' function.
#'
#' @param dna_strings A character vector of DNA sequences. All strings must
#'   have the same length.
#' @return A data frame where each row corresponds to an input DNA string,
#'   and columns represent the one-hot encoded nucleotides at each position.
#' @keywords internal

dna_encoding <- function(dna_strings){
  nn <- nchar( dna_strings[1] )
  seq_m <- matrix( unlist( strsplit(dna_strings, "") ), ncol = nn, byrow = TRUE)
  colnames(seq_m) <- paste0("nt_pos", 1:nn)
  seq_df <- as.data.frame(seq_m)
  seq_df[] <- lapply(seq_df, factor, levels = c("A", "T", "C", "G"))
  return(seq_df)
}

#' Predict m6A Modification Status for Multiple Sites
#'
#' This function takes a data frame of input features and predicts the m6A
#' modification probability and status for each site using a pre-trained
#' random forest model.
#'
#' @param ml_fit A fitted random forest model object, typically loaded from an
#'   .rds file.
#' @param feature_df A data frame containing the input features for prediction.
#'   Must include columns: `gc_content`, `RNA_type`, `RNA_region`,
#'   `exon_length`, `distance_to_junction`, `evolutionary_conservation`, and
#'   `DNA_5mer`.
#' @param positive_threshold A numeric value between 0 and 1. The probability
#'   threshold above which a site is classified as "Positive". Defaults to 0.5.
#'
#' @return The input `feature_df` augmented with two new columns:
#'   `predicted_m6A_prob` (the probability of being a positive m6A site) and
#'   `predicted_m6A_status` ("Positive" or "Negative").
#'
#' @import randomForest
#' @importFrom stats predict
#' @export
#'
#' @examples
#' # Load the pre-trained model from the package's external data
#' rf_model <- readRDS(
#'   system.file("extdata", "rf_fit.rds", package = "m6APrediction")
#' )
#'
#' # Load the example input data from the package
#' example_df <- read.csv(
#'   system.file("extdata", "m6A_input_example.csv", package = "m6APrediction")
#' )
#'
#' # Run the prediction function with the loaded data
#' predictions <- prediction_multiple(rf_model, example_df, positive_threshold = 0.6)
#'
#' # Print the first few rows of the results
#' head(predictions)
#'

prediction_multiple <- function(ml_fit, feature_df, positive_threshold = 0.5){
  stopifnot(all(c("gc_content", "RNA_type", "RNA_region", "exon_length", "distance_to_junction", "evolutionary_conservation", "DNA_5mer") %in% colnames(feature_df))) #Check errors if incorrect column names of input data.frame

  feature_df <- cbind(feature_df, dna_encoding(feature_df$DNA_5mer))
  feature_df$RNA_type <- factor(feature_df$RNA_type, levels = c("mRNA", "lincRNA", "lncRNA", "pseudogene"))
  feature_df$RNA_region <- factor(feature_df$RNA_region, levels = c("CDS", "intron", "3'UTR", "5'UTR"))

  probs <- predict(ml_fit, newdata = feature_df, type = "prob")[, "Positive"]
  status <- ifelse(probs > positive_threshold, "Positive", "Negative")

  feature_df$predicted_m6A_prob <- probs
  feature_df$predicted_m6A_status <- factor(status, levels = c("Negative", "Positive"))

  return(feature_df) #return a data.frame with supplied columns of predicted m6A prob and predicted m6A status
}


#' Predict m6A Modification Status for a Single Site
#'
#' This function predicts the m6A modification probability and status for a
#' single site based on manually provided feature values.
#'
#' @param ml_fit A fitted random forest model object.
#' @param gc_content Numeric value for GC content (0-1).
#' @param RNA_type Character string for RNA type (e.g., "mRNA").
#' @param RNA_region Character string for RNA region (e.g., "CDS").
#' @param exon_length Numeric value for exon length.
#' @param distance_to_junction Numeric value for distance to junction.
#' @param evolutionary_conservation Numeric value for evolutionary conservation (0-1).
#' @param DNA_5mer A 5-character string representing the DNA sequence.
#' @param positive_threshold A numeric value (0-1) used as the classification
#'   threshold. Defaults to 0.5.
#'
#' @return A named character vector of length 2, containing the
#'   `predicted_m6A_prob` and `predicted_m6A_status`.
#'
#' @export
#'
#' @examples
#' # Load the pre-trained model from the package's external data
#' rf_model <- readRDS(
#'   system.file("extdata", "rf_fit.rds", package = "m6APrediction")
#' )
#'
#' # Run the prediction for a single data point
#' single_prediction <- prediction_single(
#'   ml_fit = rf_model,
#'   gc_content = 0.5,
#'   RNA_type = "mRNA",
#'   RNA_region = "CDS",
#'   exon_length = 10,
#'   distance_to_junction = 8,
#'   evolutionary_conservation = 0.5,
#'   DNA_5mer = "GGACA"
#' )
#'
#' # Print the result
#' print(single_prediction)
#'

prediction_single <- function(ml_fit, gc_content, RNA_type, RNA_region, exon_length, distance_to_junction, evolutionary_conservation, DNA_5mer, positive_threshold = 0.5){

  single_df <- data.frame(gc_content = gc_content,
                          RNA_type = factor(RNA_type, levels = c("mRNA", "lincRNA", "lncRNA", "pseudogene")),
                          RNA_region = factor(RNA_region, levels = c("CDS", "intron", "3'UTR", "5'UTR")),
                          exon_length = exon_length,
                          distance_to_junction = distance_to_junction,
                          evolutionary_conservation = evolutionary_conservation,
                          DNA_5mer = DNA_5mer,
                          stringsAsFactors = FALSE)

  single_df <- cbind(single_df, dna_encoding(single_df$DNA_5mer))
  pred_df <- prediction_multiple(ml_fit, single_df, positive_threshold)

  returned_vector <- c(predicted_m6A_prob = pred_df$predicted_m6A_prob[1],
                       predicted_m6A_status = as.character(pred_df$predicted_m6A_status[1]))

  return(returned_vector) #return a named vector with values for predicted m6A prob and predicted m6A status
}
