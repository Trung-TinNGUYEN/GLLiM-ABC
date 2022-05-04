#'@export
get_gllim_colors <- function(){
  return(c("MW2"="black", "MW2b"="black","L2"= "blue", "L2b"="blue", "SA"="green",
  "SAb"="green","GE"="red", "GEV"="red", "Gmixt"="green",
           "Hilbert" = "orange", "marginal" = "red", "Curve matching" = "darkblue",
           "MMD"="009E73", "Posterior"="black",
           "Summary"="orange",
           "Semi-auto"="orange",
           "SA + constraint"="red",
           "Wasserstein"="darkblue",
           "Euclidean"="#009E73",
           "Swap" = "purple",
           "Rej. Summary" = "purple",
           "W + constraint" = "darkgreen"
  ))
}