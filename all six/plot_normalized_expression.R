library(Seurat)
library(ggplot2)

plot_normalized_expression <- function(s.object, gene) {
  # Fetch expression
  expr_data <- FetchData(s.object, vars = gene)
  
  # Rename column to a consistent name
  colnames(expr_data) <- "expression"
  
  # Check for negative values
  has_negatives <- any(expr_data$expression < 0)
  cat("Any negative expression values for", gene, "?", has_negatives, "\n")
  
  # Add cell name and index
  expr_data$cell <- rownames(expr_data)
  expr_data$index <- seq_len(nrow(expr_data))
  
  # Plot
  p <- ggplot(expr_data, aes(x = index, y = expression)) +
    geom_point(alpha = 0.6, color = "darkblue") +
    labs(title = paste(gene, "Expression Across Cells (log-normalized)"),
         x = "Cell Index",
         y = paste(gene, "Expression")) +
    theme_minimal()
  
  print(p)
}
