

#' ### Riqueza capturada según estimadores
specpool(mc_apcyn_melic_saptc)
specpool(mc_apcyn_melic_saptc)[[1]]/specpool(mc_apcyn_melic_saptc)*100

#SCBD (species contribution to beta diversity) y LCBD (local contribution...)
mc_apcyn_melic_saptc_beta <- beta.div(mc_apcyn_melic_saptc, method = "hellinger", nperm = 9999)
mc_apcyn_melic_saptc_beta$SCBD[mc_apcyn_melic_saptc_beta$SCBD >= mean(mc_apcyn_melic_saptc_beta$SCBD)]
row.names(mc_apcyn_melic_saptc[which(mc_apcyn_melic_saptc_beta$p.LCBD <= 0.05),])
#
