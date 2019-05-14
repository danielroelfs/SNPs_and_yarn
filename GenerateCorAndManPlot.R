### Libraries ####

library(tidyverse)
library(danielR)
library(normentR)

### Create sample correlation data ####

# Create base parameters
nnodes <- 21
x_vector <- seq(1,nnodes,1)

radius<-rep(9,nnodes)

name<-as.character(1:nnodes)
labelnames <- c("DMN 1","DMN 2","RL","LL","SomMot","SomSens","PreMot","PrefCort","Audit 1","Audit 2","Cing","RCE",
                "VentStr","Front","Cereb","Thala","TempPar","Visu 1","Visu 2","Visu 3","Visu 4")

bdat <- data.frame(x_vector,radius,name)
bdat$label_radius <- bdat$radius + 1

# Create stats parameters
nedges <- 210
from <- as.character(sample(1:nnodes,nedges,replace = TRUE))
to <- as.character(sample(1:nnodes,nedges,replace = TRUE))
target <- sample(seq(1,nnodes,2),nedges,replace = TRUE)
strength <- rnorm(nedges,mean=1,sd=0.25)
to_x <- sample(seq(3,nnodes,2),nedges,replace = TRUE)
to_rad <- rep(9,nedges)
from_x <- sample(1:nnodes,nedges,replace = TRUE)
from_rad <- rep(9,nedges)

cdat <- data.frame(from,to,target,strength,to_x,to_rad,from_x,from_rad)

# Add radial points
bdat2 <- bdat %>%
  mutate(x_vector = as.integer(factor(x_vector))) %>%
  mutate(theta = x_vector / n_distinct(x_vector) * 2 * pi + pi / 2) %>%
  mutate(x = radius * cos(theta),
         y = radius * sin(theta),
         y.label = label_radius * sin(theta),
         name = as.character(name))

# Combine two data frames
cdat2 <- cdat %>%
  select(from, to, target, strength) %>%
  mutate_at(vars(from, to), as.character) %>%
  left_join(bdat2 %>% select(name, x, y), 
            by = c("from" = "name")) %>%
  rename(x.start = x, y.start = y) %>%
  left_join(bdat2 %>% select(name, x, y),
            by = c("to" = "name")) %>%
  rename(x.end = x, y.end = y)

# Remove points where nodes connect to themselves
plotdat <- cdat2[which(cdat2[ ,"from"] != cdat2[, "to"]), ]

# Set plot range
plot.range <- max(abs(c(bdat2$x, bdat2$y, bdat2$y.label))) * 1.1

### Plot network correlation ####

cor.plot <- ggplot(bdat2, aes(x = x, y = y)) +
  geom_curve(data = plotdat %>% filter(x.start > 0),
             aes(x = x.start, y = y.start, 
                 xend = x.end, yend = y.end, 
                 color = strength),
             size = 2, curvature = -0.3) +
  geom_curve(data = plotdat %>% filter(x.start <= 0),
             aes(x = x.start, y = y.start,
                 xend = x.end, yend = y.end,
                 color = strength),
             size = 2, curvature = 0.3) +
  geom_point(size = 8, color = "grey90") +
  #geom_text(aes(y = y.label, label = labelnames)) +
  expand_limits(x = c(-plot.range, plot.range),
                y = c(-plot.range, plot.range)) +
  scale_color_daniel(discrete = FALSE, palette = "vik") +
  coord_equal() +
  theme_blankcanvas() +
  theme(
    rect = element_rect(fill = "transparent"),
    legend.position = 'none',
    plot.margin = unit(rep(0, 4), "cm")
  )
cor.plot
ggsave('corplot4.png',cor.plot,width = 50, height = 50, 
       dpi = 300, units = 'cm',
       bg = "transparent")

### Create Manhattan data ####

gwas.dat <- simulateGWAS(nSNPs = 1e5, 
                         N = 2e4, 
                         AddSigSNPs = TRUE, 
                         SigCHR = NULL, 
                         nSigCols = 6)
gwas.dat <- gwas.dat[gwas.dat$P < 1, ]

# Prepare data for plotting
nCHR <- length(unique(gwas.dat$CHR))
gwas.dat$BPcum <- NA
s <- 0
nbp <- c()

# Add cumulative base pair position
for (i in unique(gwas.dat$CHR)){
  nbp[i] <- max(gwas.dat[gwas.dat$CHR == i,]$BP)
  gwas.dat[gwas.dat$CHR == i,"BPcum"] <- gwas.dat[gwas.dat$CHR == i,"BP"] + s
  s <- s + nbp[i]
}

# Set plot parameters
axis.set <- gwas.dat %>% group_by(CHR) %>% summarize(center=(max(BPcum) + min(BPcum)) / 2)
ylim <- abs(floor(log10(min(gwas.dat$P)))) + 2

# Add significance thresholds
supersig <- 0.05 / 1e6 / nedges
sig <- 0.05 / 1e6
ssig <- 0.1 / 1e4

### Plot Manhattan ####

# Set colors
getcols <- daniel_pal(palette = "vik")
getcols_rev <- daniel_pal(palette = "vik",reverse = TRUE)
#allcols <- getcols(256)
#c1 <- allcols[50]
#c2 <- allcols[256-50]
#c2 <- "#BD4648"
#c <- getcols(4)
c <- c(getcols(nCHR/2),getcols_rev((nCHR/2)+1))
c <- c[-12]

# Create plot object
manh.plot <- ggplot() +
  geom_point(gwas.dat, mapping = aes(x=BPcum, y=-log10(P),color=as.factor(CHR)), 
             alpha = 1, size = 1.25) +
  #scale_color_manual(values = rep(c(c1,c2), nCHR)) +
  scale_color_manual(values = c) +
  #scale_color_manual(values = getcols(nCHR)) +
  scale_x_continuous(label = axis.set$CHR, breaks = axis.set$center) +
  scale_y_continuous(expand = c(0,0), limits = c(-12,ylim)) +
  geom_hline(yintercept = -log10(supersig), color = "grey80", linetype = "dashed") + 
  geom_hline(yintercept = -log10(sig), color = "grey60", linetype = "dashed") + 
  geom_hline(yintercept = -log10(ssig), color = "grey40", linetype = "dashed") + 
  labs(x = NULL, y = "-log10(p)") + 
  coord_polar() +
  theme_blankcanvas() +
  theme( 
    rect = element_rect(fill = "transparent"),
    axis.text.x = element_blank(),
    legend.position = "none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    plot.margin = unit(rep(0, 4), "cm")
  )
manh.plot
ggsave('manhplot4.png',manh.plot,width = 50, height = 50, 
       dpi = 300, units = 'cm',
       bg = "transparent")
