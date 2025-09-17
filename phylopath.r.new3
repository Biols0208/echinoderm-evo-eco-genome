library(phylopath)
library(ape)
library(ggplot2)
 
# Tree
tree <- read.tree('ref.tree')
my_tree = drop.tip(tree, c("Branchiostoma_floridae", "Henricia", "Pecten_maximus"))

# Data
my_data <- read.csv("genome_Phenotype.csv")
rownames(my_data) <- my_data$Species

# Convert to ordered numeric variable
relative_size_ordered <- as.numeric(factor(my_data$Relative_size,
                                          levels = c("small", "medium_to_small", "medium", "large", "big"),ordered = TRUE))
habit_ordered <- as.numeric(factor(my_data$Habitat,
                                          levels = c("Soft", "Hard", "Mixed", "Specialized"),ordered = TRUE))
diet_ordered <- as.numeric(factor(my_data$Diet_Class,
										  levels = c("Suspension_Feeder", "Deposit_Feeder", "Herbivore", "Carnivore", "Omnivore"),ordered = TRUE))	

my_data$is_Asteroidea <- ifelse(my_data$group == 'Asteroidea', 1, 0)
my_data$is_Ophiuroidea <- ifelse(my_data$group == 'Ophiuroidea', 1, 0)
my_data$is_Crinoidea <- ifelse(my_data$group == 'Crinoidea', 1, 0)
my_data$is_Holothuroidea <- ifelse(my_data$group == 'Holothuroidea', 1, 0)
my_data$is_Echinoidea <- ifelse(my_data$group == 'Echinoidea', 1, 0)

# Create environmental adaptability indices
temp_adaptability <- my_data$suitable_temperature_range / pmax(my_data$suitable_temperature_avarage, 1)
sal_adaptability <- my_data$suitable_salinity_range / pmax(my_data$suitable_salinity_avarage, 1)
depth_adaptability <- my_data$suitable_depth_range / pmax(my_data$suitable_depth_avarage, 1)

# Environmental stress index
env_stress <- scale(temp_adaptability + sal_adaptability + depth_adaptability)[,1]

# Metabolic demand index
metabolic_demand <- scale(relative_size_ordered * habit_ordered * diet_ordered)[,1]

# Developmental complexity score
dev_complexity <- (ifelse(my_data$Reproductive_Mode == "planktotrophic", 3, 1) + ifelse(my_data$Adult_Strategy == "K", 0, 1))

# Niche breadth (environmental tolerance)
niche_breadth <- scale(temp_adaptability + sal_adaptability + depth_adaptability)[,1]

# Specialization index (inverse of niche breadth)
specialization <- -niche_breadth

## select and rename
my_data_ne = data.frame(my_data$genomesize, my_data$GC, my_data$DNA, my_data$LTR, my_data$LINE, my_data$Total,
                        relative_size_ordered, 
						as.numeric(ifelse(my_data$Reproductive_Mode == "planktotrophic", 2, 1)) + rnorm(length(my_data$Reproductive_Mode), 0, 0.001),
                        as.numeric(ifelse(my_data$Larval_Strategy == "R", 2, 1)) + rnorm(length(my_data$Reproductive_Mode), 0, 0.001),  
						as.numeric(ifelse(my_data$Adult_Strategy == "K", 0, 1)) + rnorm(length(my_data$Reproductive_Mode), 0, 0.001),
                        diet_ordered, habit_ordered, 

                        my_data$suitable_temperature_avarage, my_data$suitable_temperature_range,
                        my_data$limit_temperature_min, my_data$limit_temperature_max,

                        my_data$suitable_salinity_avarage, my_data$suitable_salinity_range,
                        my_data$limit_salinity_min, my_data$limit_salinity_max,

                        my_data$suitable_depth_avarage, my_data$suitable_depth_range,
                        my_data$limit_depth_max,
						
						temp_adaptability, sal_adaptability, depth_adaptability, env_stress, metabolic_demand, dev_complexity, niche_breadth, specialization,
						as.factor(my_data$is_Asteroidea), as.factor(my_data$is_Ophiuroidea), 
						as.factor(my_data$is_Crinoidea), as.factor(my_data$is_Holothuroidea), as.factor(my_data$is_Echinoidea))

rownames(my_data_ne) <- my_data$Species
# genomesize:GE Total:TE Relative_size:Rsize Reproductive_Mode:RM Larval_Strategy:LS Adult_Strategy:AS Diet_Class:DC Habitat:HA
# suitable_temperature_avarage:STA suitable_temperature_range:STR limit_temperature_min:LTN limit_temperature_max:LTX
# suitable_salinity_avarage:SSA suitable_salinity_range:SSR limit_salinity_min:LSN limit_salinity_max:LSX
# suitable_depth_avarage:SDA suitable_depth_range:SDR limit_depth_max:LDX
# temp_adaptability:TA sal_adaptability:SA depth_adaptability:DA env_stress:ES metabolic_demand:MD dev_complexity:DeC niche_breadth:NB specialization:SP
colnames(my_data_ne) <- c("GE", "GC", "DNA", "LTR", "LINE", "TE",
                          "RS", "RM", "LS", "AS", "DC", "HA",
                          "STA", "STR", "LTN", "LTX",
                          "SSA", "SSR", "LSN", "LSX",
                          "SDA", "SDR", "LDX",
						  "TA", "SA", "DA", "ES", "MD", "DeC", "NB", "SP",
						  "ISA", "ISO", "ISC", "ISH", "ISE")
#head(my_data_ne)

## check data quality 
## If the independent variables are highly correlated (e.g., |r| > 0.9), putting them together in a model may lead to unstable results. This check can help you determine whether you need to remove or combine certain variables when building the model.
check_data_quality <- function(data) {
  cat("Data Quality Check:\n")
  cat("------------------\n")

  # Check for missing values
  missing_per_col <- colSums(is.na(data))
  if (any(missing_per_col > 0)) {
    cat("Columns with missing values:\n")
    print(missing_per_col[missing_per_col > 0])
  } else {
    cat("No missing values detected\n")
  }

  # Check for zero variance columns
  numeric_cols <- sapply(data, is.numeric)
  zero_var <- sapply(data[, numeric_cols], var, na.rm = TRUE) == 0
  if (any(zero_var)) {
    cat("\nColumns with zero variance:\n")
    print(names(zero_var)[zero_var])
  }

  # Check for high correlations
  cor_matrix <- cor(data[, numeric_cols], use = "complete.obs")
  high_cor <- which(abs(cor_matrix) > 0.9 & cor_matrix != 1, arr.ind = TRUE)
  if (nrow(high_cor) > 0) {
    cat("\nHighly correlated variable pairs (|r| > 0.9):\n")
    for (i in 1:nrow(high_cor)) {
      if (high_cor[i,1] < high_cor[i,2]) {  # Avoid duplicates
        cat(sprintf("  %s - %s: %.3f\n",
                   rownames(cor_matrix)[high_cor[i,1]],
                   colnames(cor_matrix)[high_cor[i,2]],
                   cor_matrix[high_cor[i,1], high_cor[i,2]]))
      }
    }
  }

  cat("\nData dimensions:", nrow(data), "rows x", ncol(data), "columns\n")
}
check_data_quality(my_data_ne)

## Test Model
# genomesize:GE Total:TE Relative_size:Rsize Reproductive_Mode:RM Larval_Strategy:LS Adult_Strategy:AS Diet_Class:DC Habitat:HA
# suitable_temperature_avarage:STA suitable_temperature_range:STR limit_temperature_min:LTN limit_temperature_max:LTX
# suitable_salinity_avarage:SSA suitable_salinity_range:SSR limit_salinity_min:LSN limit_salinity_max:LSX
# suitable_depth_avarage:SDA suitable_depth_range:SDR limit_depth_max:LDX
# temp_adaptability:TA sal_adaptability:SA depth_adaptability:DA env_stress:ES metabolic_demand:MD dev_complexity:DeC niche_breadth:NB specialization:SP
# Highly correlated variable pairs:
# RM-LS:1.000   STA-LTX: 0.926   SDR-LDX:0.907   SSR-SA:0.996   RM-DeC:0.947   LS-DeC:0.947   ES-SP:-1.000   NB-SP:-1.000 

models_ne <- define_model_set(
  null = c(),
  # 假设：发育复杂度 -> 体型 -> 基因组大小
  LifeHistory_Simple = c(GE ~ RS, RS ~ DeC),
  # 假设：繁殖模式和食性共同影响体型，体型和繁殖模式再影响基因组大小
#  LifeHistory_Tradeoff = c(GE ~ RS + RM, RS ~ RM + DC),
  # 假设：体型和食性决定代谢需求，代谢需求通过影响TEs进而影响基因组大小
  Metabolic_Mediation_A = c(GE ~ TE, TE ~ RS + DC + HA + ISA + ISO + ISC + ISH + ISE), # 假说A: 纲影响TE
  Metabolic_Mediation_B = c(GE ~ TE + ISA + ISO + ISC + ISH + ISE, TE ~ RS + DC + HA), # 假说B: 纲影响GE
  # 假设：繁殖与成体策略决定发育复杂度，后者通过TEs影响基因组 (erro: RM-DeC:0.947   LS-DeC:0.947)
  Development_MediationA = c(GE ~ DeC + TE, TE ~ DeC + ISA + ISO + ISC + ISH + ISE),
  Development_MediationB = c(GE ~ DeC + TE + ISA + ISO + ISC + ISH + ISE, TE ~ DeC),
  # 假设：环境压力直接影响基因组大小,环境相关模型基准
#  Environment_Simple = c(GE ~ ES),
  # 假设：环境压力激活TEs，两者共同影响基因组大小
  Environment_Mediation_A = c(GE ~ ES + TE, TE ~ ES + ISA + ISO + ISC + ISH + ISE), # 假说A: 纲影响TE
  Environment_Mediation_B = c(GE ~ ES + TE + ISA + ISO + ISC + ISH + ISE, TE ~ ES), # 假说B: 纲影响GE
  # 假设：具体的环境梯度（温盐深）直接影响基因组大小，且这些梯度间存在关联
#  Environment_Gradient = c(GE ~ STA + SSA + SDA + ISA + ISO + ISC + ISH + ISE, STA ~ SDA, SSA ~ SDA),
  # 假设：更宽的生态位（更强的环境适应力）通过影响TEs，进而影响基因组
#  Niche_Breadth = c(GE ~ NB + TE, TE ~ NB),
  # 假设：基因组大小同时受体型和环境压力的影响
#  Integrated_Size_Env = c(GE ~ RS + ES + ISA + ISO + ISC + ISH + ISE, RS ~ RM + HA),
  # 整合多个生活史和环境因素
#  Comprehensive_LH_Web = c(GE ~ RS + RM + AS + ISA + ISO + ISC + ISH + ISE, RS ~ RM + DC + HA, RM ~ STA, AS ~ DC + HA),
  # 假设：物种能耐受的生理极限（最高/最低温盐等）直接塑造了基因组大小
#  Env_Limits_Direct = c(GE ~ LTX + LTN + LSX + LSN + LDX + ISA + ISO + ISC + ISH + ISE),
  # 假设：不是平均值或极限值，而是物种的生态适应范围大小（温度范围，盐度范围等）
#  Env_Ranges_Direct = c(GE ~ STR + SSR + SDR),
  # 假设：极端温盐深可能会激活TEs，进而影响基因组大小 
  Env_LimitsA = c(GE ~ TE, TE ~ LTX + LTN + LSX + LSN + LDX + ISA + ISO + ISC + ISH + ISE),
  Env_LimitsB = c(GE ~ TE + ISA + ISO + ISC + ISH + ISE, TE ~ LTX + LTN + LSX + LSN + LDX),
  Integrated_Eco_Evo_Devo = c(GE ~ TE + DeC + ISA + ISO + ISC + ISH + ISE, TE ~ DeC + LTX + LTN + LSX + LSN + LDX, DeC ~ LTX + LTN + LSX + LSN + LDX)
)

result <- phylo_path(models_ne, data = my_data_ne, tree = my_tree, model = 'lambda', method = "logistic_MPLE")
result

#result$d_sep$Development_MediationB
print(result$d_sep$Development_MediationB, n = Inf) 

#result$d_sep$Env_LimitsB
print(result$d_sep$Env_LimitsB, n = Inf) 

warnings()
s <- summary(result)
print(s, n = Inf)

pdf("out3.pdf", width = 16, height = 20)

# Safer plotting helpers to avoid ggraph arc errors
safe_print_plot <- function(p) {
  tryCatch({
    if (inherits(p, "gg") || inherits(p, "ggplot")) {
      print(p)
    } else {
      # Many phylopath plot() calls draw directly; try printing anyway
      try(print(p), silent = TRUE)
    }
  }, error = function(e) {
    message("Plot skipped due to error: ", e$message)
    plot.new(); title(main = "Plot skipped due to plotting error")
  })
}

# 1) Model set structure (robust to failures)
g_ms <- try(plot_model_set(models_ne), silent = TRUE)
if (!inherits(g_ms, "try-error")) {
  safe_print_plot(g_ms)
} else {
  plot.new(); title(main = "Model set plot failed; skipped")
}

# 2) Model comparison
safe_print_plot(try(plot(s), silent = TRUE))

# 3) Best model (force straight edges and robust layout)
best_model <- best(result, boot = 500)
#safe_print_plot(try(plot(best_model, curvature = 0.5), silent = TRUE))
safe_print_plot(try(plot(best_model, algorithm = 'nicely', curvature = 0.8), silent = TRUE))

# 4) Average model (force straight edges and robust layout)
average_model <- average(result)
#average_model_full <- average(result, avg_method = "full")
coef_plot(average_model, order_by = "strength", to = "GE") + ggplot2::coord_flip() + ggplot2::theme_bw()
coef_plot(average_model, order_by = "strength", to = "TE") + ggplot2::coord_flip() + ggplot2::theme_bw()
#coef_plot(average_model_full, order_by = "strength", to = "GE") + ggplot2::coord_flip() + ggplot2::theme_bw()
#safe_print_plot(try(plot(average_model, curvature = 0.2), silent = TRUE))
safe_print_plot(try(plot(average_model, algorithm = 'nicely', curvature = 0.8), silent = TRUE))
#safe_print_plot(try(plot(average_model, algorithm = 'star', curvature = 0.5), silent = TRUE))

#positions <- data.frame(
#  name = c("GE", "MD", "TE", "ISA", "ISO", "ISC", "ISH", "ISE", "RS", "DC", "LTX", "LTN", "LSX", "LSN", "LDX"),
#  x= c(5,6,4,3,4,5,6,8,8,8,3,4,5,6,7),
#  y= c(3,2,2,4,4,4,4,4,2.5,3.5,1,1,1,1,1))
#safe_print_plot(try(plot(average_model, algorithm = 'mds', curvature = 0, manual_layout = positions), silent = TRUE))
#safe_print_plot(try(plot(best_model, algorithm = 'mds', curvature = 0, manual_layout = positions), silent = TRUE))
dev.off()
