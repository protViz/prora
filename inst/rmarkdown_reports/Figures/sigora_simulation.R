
# replicate(n = 100, {
#   seed <- round(runif(1,0,10000))
#   a1 <- genesFromRandomPathways(
#     GPSrepo = kegH,
#     np = 3, # From np random pathways
#     ng = 50,# Draw nq genes at random
#     seed = seed
#   )
#   cbind(
#     ora = nrow(ora(a1[["genes"]], kegH)),
#     sigora = sigora(
#       GPSrepo = kegH,
#       queryList = a1[["genes"]],
#       level = 4
#     )[[1]] %>% filter(Bonferroni <= 0.05) %>% nrow)
# }) %>% apply(., 2, cbind) -> sim_out

# save(sim_out, file = "~/Dropbox/FGCZ/Figures/simulation_sigora.Rda")

load(file = "~/Dropbox/FGCZ/Figures/simulation_sigora.Rda")
mp <- barplot(colMeans(sim_out), ylim = c(0,80), ylab = "Significantly enriched pathways after Bonferroni correction")
arrows(x0 = mp, y0 = colMeans(sim_out), x1 = mp, y1 = colMeans(sim_out)+apply(sim_out, 2, sd), angle = 90)
