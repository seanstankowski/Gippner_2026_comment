################# plot a set of clines for variying values of s with fixed sigma ################# 

# transect
x <- seq(-100, 100, length.out = 1000)

# value of dispersal
sigma <- 2

# Different strengths of selection
s_values <- c(0.005, 0.005, 0.025, 0.1, 0.4) # first value is a dummy position used to set up the plot. Its plotted white

# Cline centre
centre <- 0

# for each s calculate the expected cline width assuming het disadvantage (coefficient proportionality of 8)
widths <- sqrt((8 * sigma^2) / s_values)

# Function for sigmoid cline
cline <- function(x, centre, width) {
  1 / (1 + exp(-4 * (x - centre) / width))
}

# assign colours for other clines
cols <- colorRampPalette(c("blue", "red"))(length(s_values))

# Plot first dummy cline white to set up plot
plot(
  x,
  cline(x, centre, widths[1]),
  type = "l",
  lwd = 3, col="white",
  ylim = c(0, 1),
  xlab = "Distance across transect",
  ylab = "Allele frequency",
)

# Add remaining clines
for (i in 2:length(s_values)) {
  lines(
    x,
    cline(x, centre, widths[i]),
    lwd = 3,
    lty = 1,
    col = cols[i]
  )
}

# Add cline centre
# abline(v = centre, lty = 2, col = "grey50")

# Legend
legend(
  "topleft",
  legend = paste0(
    "s = ", s_values,
    ", w = ", round(widths, 1)
  ),
  lwd = 3,
  lty = seq_along(s_values),
  bty = "n"
)


################ How staggered clines broaden the cline for the hybrid index ################# 

# transect
x <- seq(-100, 100, length.out = 1000)

# value of dispersal
sigma <- 2

# fixed strength of selection
s <- 0.4

# number of clines
n_clines <- 100

# sample cline centres around the midpoint of the transect with sd of 25
centres <- rnorm(
  n = n_clines,
  mean = 0,
  sd = 25
)

# order for nicer plotting with colour gradient
centres <- sort(centres)

# expected width for each individual cline
width <- sqrt((8 * sigma^2) / s)

# Function for sigmoid cline
cline <- function(x, centre, width) {
  1 / (1 + exp(-4 * (x - centre) / width))
}

# Make matrix of individual clines
cline_matrix <- sapply(centres, function(c) {
  cline(x, centre = c, width = width)
})

# Hybrid index: average across clines
hybrid_index <- rowMeans(cline_matrix)

# Colours for individual clines
cols <- colorRampPalette(c("blue", "red"))(length(centres))

# Plot blank canvas
plot(
  x,
  hybrid_index,
  type = "n",
  ylim = c(0, 1),
  xlab = "Distance across transect",
  ylab = "Allele frequency / hybrid index"
)

# Add individual clines
for (i in seq_along(centres)) {
  lines(
    x,
    cline_matrix[, i],
    lwd = 2,
    col = cols[i]
  )
}

# Add hybrid index
lines(
  x,
  hybrid_index,
  lwd = 5,
  col = "black"
)


######## geographic overlap of two taxa with no hybridization
## this is just achived via binomial samling along the sigmoid, not that the hybrids are added to the plot manually in adobe illustrator

# -----------------------------
# Parameters
# -----------------------------

n_demes <- 100
n_ind <- 20

centre <- 50
width <- 10 ##### this controls the amount of overlap of the lineages

# deme positions
demes <- 1:n_demes

# -----------------------------
# Expected cline
# -----------------------------

# probability of sampling a type-1 individual in each deme
p_type1 <- 1 / (1 + exp(-4 * (demes - centre) / width))

# -----------------------------
# Binomial sampling
# -----------------------------

# number of type-1 individuals per deme
n_type1 <- rbinom(
  n = n_demes,
  size = n_ind,
  prob = p_type1
)

# number of type-0 individuals per deme
n_type0 <- n_ind - n_type1

# observed allele frequency / mean hybrid index per deme
mean_type1 <- n_type1 / n_ind

# -----------------------------
# Individual-level data for plotting
# -----------------------------

deme_id <- rep(demes, each = n_ind)

type <- unlist(
  lapply(seq_along(demes), function(i) {
    c(rep(0, n_type0[i]), rep(1, n_type1[i]))
  })
)

# -----------------------------
# Plot individuals and cline
# -----------------------------

plot(
  deme_id + runif(length(deme_id), -0.25, 0.25),
  type,
  pch = 16,
  col = ifelse(type == 0, "dodgerblue3", "firebrick3"),
  ylim = c(-0.1, 1.1),
  xlab = "Deme",
  ylab = "p",
  yaxt = "n",
)

axis(
  side = 2,
  at = c(0, 1),
  labels = c("0", "1"),
  las = 1
)

# expected sigmoid cline
lines(
  demes,
  p_type1,
  lwd = 4,
  col = "black",
  lty = 1
)


