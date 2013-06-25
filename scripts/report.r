df <- read.csv("results_lithium.csv", header=TRUE)

indices <- df$Fusion.Events <= 0

ambient   <- df$Ambient.Radius
max_r     <- df$Maximum.Radius
wall_temp <- df$Wall.Temperature
min_r     <- df$Minimum.Radius
wall_v    <- df$Max.Wall.Velocity
ratio     <- df$Expansion.Ratio
fusions   <- df$Fusion.Events
expansion <- max_r / ambient

log_fusions <- log(fusions)
log_fusions[indices] = 0.0

reg <- lm(log_fusions ~ ambient + expansion + wall_temp)

summary(reg)

pdf("example1.pdf")
layout(matrix(c(1,2,3,4),2,2)) # optional 4 graphs/page
plot(ambient, log(fusions))
plot(expansion, log(fusions))
plot(wall_temp, log(fusions))
#plot(reg)
dev.off()
