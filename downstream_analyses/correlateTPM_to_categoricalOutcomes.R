# Load necessary libraries
library(dplyr)
library(tidyr)
library(ggplot2)

# Generate example data
set.seed(123)
gene_names <- paste0("gene", 1:100)
expression_data <- matrix(rnorm(300, mean = 5, sd = 1), ncol = 100)
colnames(expression_data) <- gene_names
outcome_data <- data.frame(
  outcome1 = sample(c("low", "medium", "high"), 100, replace = TRUE),
  outcome2 = sample(c("low", "medium", "high"), 100, replace = TRUE),
  outcome3 = sample(c("low", "medium", "high"), 100, replace = TRUE)
)
rownames(outcome_data) <- gene_names

# Convert expression data to tidy format
expression_data_tidy <- expression_data %>%
  as.data.frame() %>%
  rownames_to_column(var = "sample_id") %>%
  pivot_longer(-sample_id, names_to = "gene_name", values_to = "tpm")

# Join expression and outcome data
merged_data <- expression_data_tidy %>%
  left_join(outcome_data %>% rownames_to_column("gene_name"), by = c("gene_name" = "gene_name"))

merged_data


# Convert categorical outcomes to numeric values for correlation calculation
merged_data$outcome1_numeric <- ifelse(merged_data$outcome1 == "low", 1, ifelse(merged_data$outcome1 == "medium", 2, 3))
merged_data$outcome2_numeric <- ifelse(merged_data$outcome2 == "low", 1, ifelse(merged_data$outcome2 == "medium", 2, 3))
merged_data$outcome3_numeric <- ifelse(merged_data$outcome3 == "low", 1, ifelse(merged_data$outcome3 == "medium", 2, 3))

# Calculate correlations
correlations <- merged_data %>%
  select(-gene_name) %>%
  group_by(sample_id) %>%
  summarise(
    corr_outcome1 = cor(tpm, outcome1_numeric),
    corr_outcome2 = cor(tpm, outcome2_numeric),
    corr_outcome3 = cor(tpm, outcome3_numeric)
  )


# Visualize correlations
ggplot(correlations, aes(x = sample_id, y = corr_outcome1)) +
  geom_point() +
  ggtitle("Correlation between gene expression TPM and outcome1") +
  xlab("Sample ID") +
  ylab("Correlation coefficient") +
  theme_bw()

ggplot(correlations, aes(x = sample_id, y = corr_outcome2)) +
  geom_point() +
  ggtitle("Correlation between gene expression TPM and outcome2") +
  xlab("Sample ID") +
  ylab("Correlation coefficient") +
  theme_bw()

ggplot(correlations, aes(x = sample_id, y = corr_outcome3)) +
  geom_point() +
  ggtitle("Correlation between gene expression TPM and outcome3") +
  xlab("Sample ID") +
  ylab("Correlation coefficient") +
  theme_bw()
