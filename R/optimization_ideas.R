library(dplyr)

# ideas on weight from chatgpt
#
# Inverse Rank Weighting:
# Assign a weight to each item based on the inverse of its rank.
# Formula: weight = 1 / rank.
# This method heavily emphasizes the top-ranked items, as their weights will be significantly higher than those of lower-ranked items.
#
# Linearly Decreasing Weights:
#   Assign weights that linearly decrease with increasing rank.
# Formula: weight = max_rank - rank + 1 (where max_rank is the highest rank in your data).
# This method offers a more balanced approach, where the difference in weights between successive ranks is constant.
#
# Exponential Decay Weighting:
#   Assign weights that exponentially decrease with increasing rank.
# Formula: weight = exp(-alpha * rank) (where alpha is a decay rate parameter).
# This method can be tuned to heavily favor top-ranked items (with a high decay rate) or to provide a more gradual decay.
#
# Normalized Rank Weighting:
#   Normalize the ranks and then take the inverse.
# Formula: weight = 1 / (rank / max_rank).
# This method normalizes the ranks to a [0, 1] scale before inverting, ensuring that weights are scaled relative to the size of the data set.
#
# Quadratic or Higher-Order Polynomial Weighting:
#   Assign weights based on a polynomial function of the rank.
# Formula: weight = (max_rank - rank + 1)^n (where n is the degree of the polynomial).
# This allows for more flexible weighting schemes, where higher-order polynomials can provide a sharper distinction between high and low ranks.
#
# Logarithmic Weighting:
#   Use the logarithm of the rank to assign weights.
# Formula: weight = 1 / log(rank + 1) (adding 1 to avoid log(0)).
# This method gives a slower decay of weights and is less aggressive than the inverse rank weighting.
#
# Custom Weight Functions:
#   Develop a custom weight function based on the specific needs of your analysis.
# This could involve combining different methods or introducing specific cutoffs where the weighting scheme changes.


# Label Genes as Responsive or Non-responsive: Use your label_genes function to categorize genes based on expr_effect_size and expr_p_value thresholds.
#
# Add Weight to Each Gene: Use your add_weight function to assign weights to genes based on binding_effect_size and a decay factor alpha.
#
# Fit a Weighted Logistic Regression Model: Use your fit_weighted_model function to fit the model and determine the optimal threshold for binding_pvalue.
#
# Calculate Hypergeometric P-Value: For each combination of thresholds and alpha, calculate the hypergeometric p-value of the overlap between the expression and binding data.
#
# Find Optimal Parameters: Choose the set of parameters (expr_effect_size_threshold, expr_pvalue_thres, and alpha) that minimize the hypergeometric p-value.


# Generate example data
set.seed(123)
gene_data <- data.frame(
  response = rnorm(1000), # Example response variable
  expr_effect_size = rnorm(1000, 0, 2),
  expr_p_value = runif(1000, 0.0001, 1),
  binding_effect_size = rnorm(1000, 2, 1),
  binding_pvalue = runif(1000, 0.0001, 1)
)

# label genes as responsive (TRUE) or non-responsive (FALSE) based on
# expr_effect_size and expr_p_value
label_genes <- function(data, effect_size_thres, expr_pvalue_thres) {
  # check that required cols are in the data
  stopifnot(all(c("expr_effect_size", "expr_p_value") %in% colnames(data)))
  data %>%
    mutate(
      response = ifelse(
        expr_effect_size > effect_size_thres & expr_pvalue < expr_pvalue_thres,
        TRUE,
        FALSE
      )
    )
}

add_weight = function(data, alpha){
  # check that required cols are in the data
  stopifnot(all(c("binding_effect_size") %in% colnames(data)))
  data %>%
    mutate(rank = rank(binding_effect_size)) %>%
    mutate(weight = exp(-alpha * rank))
}

# Function to fit weighted logistic regression model and return optimal threshold
# on the binding_pvalue given the ROC curve
fit_weighted_model <- function(data) {
  # check that required cols are in the data
  stopifnot(all(c("response", "binding_pvalue", "weight") %in% colnames(data)))
  # Fit logistic regression model
  model <- glm(response ~ binding_pvalue, family = binomial(), data = data, weights = weight)

  # Predict probabilities
  predictions <- predict(model, type = "response")

  # Determine the optimal threshold
  # Here we'll use a simple method based on ROC curve analysis
  library(pROC)
  roc_curve <- roc(data$response, predictions)
  coords <- coords(roc_curve, "best", best.method="closest.topleft")

  # Return the optimal threshold
  return(coords$threshold)
}

# Function to calculate hypergeometric P-value
calculate_hypergeometric <- function(expr_data, chip_data, total_genes) {
  overlap <- length(intersect(expr_data$gene, chip_data$gene))
  phyper(overlap - 1, length(expr_data$gene), total_genes - length(expr_data$gene), length(chip_data$gene), lower.tail = FALSE)
}

# Iterate over thresholds
optimal_p_value <- 1
optimal_threshold <- NA

for (threshold in seq(0.1, 2, 0.1)) {
  significant_genes <- fit_weighted_model(gene_data, threshold)
  p_value <- calculate_hypergeometric(significant_genes, significant_genes, nrow(gene_data))
  if (p_value < optimal_p_value) {
    optimal_p_value <- p_value
    optimal_threshold <- threshold
  }
}

# Output optimal threshold and P-value
list(
  optimal_threshold = optimal_threshold,
  optimal_p_value = optimal_p_value
)
