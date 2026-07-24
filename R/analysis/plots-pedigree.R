
# Load required libraries
library(tidygraph)
library(ggraph)
library(ggplot2)

pedigree_key <-
  data.frame(
    name = c(
      "YICKADEE", "BAROSSA", "HAYAC",
      "SHIRALEE", "Sapphire", "RAINBOW",
      "CHIKUZEN", "CHISAYA", "EUREKA",
      "Zephyr", "Bronowski", "MUTU",
      "Oscar", "AG-Castle", "Eyre",
      "CHARLTON", "Georgie", "KAROO",
      "MALUKA", "Monty", "SCOOP", "Mystic",
      "SIREN", "Surpass 600", "TI 10",
      "Tornado TT", "CBWA-Trilogy & Trilogy",
      "Tribune", "Wesbrook", "Wesroona",
      "Bronowski", "CRESOR", "DUNKELD",
      "ERGLU", "Tower"
    ),
    translation = c(
      "Yickadee", "Barossa", "Haya",
      "Shiralee", "Sapphire", "Rainbow",
      "Chikuzen", "Chisaya", "Eureka"
    ))

library(ggrepel)

nodes <-
  read.csv(
    "nodes.csv",
    header = TRUE,
    stringsAsFactors = FALSE
  )
edges <-
  read.csv(
    "edges.csv",
    header = TRUE,
    stringsAsFactors = FALSE
  )

# Function to abbreviate cross names while keeping cultivar names full
abbreviate_cross <- function(name) {
  if (grepl("\\[", name)) {  # It's a cross (contains brackets)
    # Replace each word with its first letter, keeping brackets and slashes
    gsub("([A-Za-z])[A-Za-z0-9-]+", "\\1", name)
  } else {
    name  # Keep cultivar/founder names as is
  }
}

# Use a layout that computes proper positions for the DAG structure
# Add abbreviated names to the nodes data
nodes$display_name <- sapply(nodes$name, abbreviate_cross)
nodes <- nodes |>
  group_by(year) |>
  mutate(
    x = (row_number() - (n() + 1) / 2) * 3,
    y = year
  ) |>
  ungroup()
# Create the graph with the updated nodes
pedigree_graph <- tbl_graph(nodes = nodes, edges = edges, directed = TRUE)
# Then create layout - it will automatically include display_name
base_layout <- create_layout(
  pedigree_graph, layout = "manual", x = nodes$x, y = nodes$y)

# Abbreviate cross names for display
base_layout$display_name <- sapply(base_layout$name, abbreviate_cross)

# Replace Y with year, but keep the computed X positions which respect graph structure
base_layout$y <- base_layout$year

# Scale X to spread nodes more horizontally, with extra spacing for years > 1969
base_layout$x <- ifelse(
  base_layout$year > 1969,
  base_layout$x * 10,  # More spacing for modern varieties
  base_layout$x * 4  # Standard spacing for founders
)

ggraph(base_layout) +
  geom_edge_diagonal(
    arrow = arrow(length = unit(2, 'mm'), type = "closed"),
    end_cap = circle(3, 'mm'),
    color = "gray50",
    alpha = 0.6
  ) +       
  geom_node_point(aes(shape = type), size = 15, fill = "white", color = "black") +  
  geom_node_text(aes(label = display_name), size = 2, hjust = 0.5, vjust = 0.5, color = "black") +
  theme_minimal() +
  scale_shape_manual(values = c("cross" = 23, "cultivar" = 21, "founder" = 22)) +  # Adjust based on your actual type values
  labs(y = "Year", x = "", title = "Pedigree DAG Over Time") +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "bottom"
  ) +
  scale_y_reverse()