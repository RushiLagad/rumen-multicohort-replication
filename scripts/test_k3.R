source("lib_replication.R")
tf <- function(r, df) r * sqrt(df) / sqrt(1 - r^2)
mk <- function(id, rn, rc, ra) rbind(
  if (!is.na(rn)) data.frame(id=id, cohort="nelore", W=tf(rn, 48), df=48, n=52),
  if (!is.na(rc)) data.frame(id=id, cohort="charolais", W=tf(rc, 26), df=26, n=30),
  if (!is.na(ra)) data.frame(id=id, cohort="alberta", W=tf(ra, 520), df=520, n=529))
est <- rbind(mk("K1", .50, .55, .45), mk("K2", .00, .00, .25), mk("K3", .02, -.03, .01), mk("K4", NA, .50, .40), mk("K5", -.50, -.55, -.45))
tg <- data.frame(id=paste0("K", 1:5), dir=c(1, 1, 1, 1, 1), testable=TRUE)
res <- classify_set(tg, est); print(res[, c("id", "k", "r", "q_bh", "loo_all_sig", "label")], digits=3, row.names=FALSE)
stopifnot(res$k == c(3, 3, 3, 2, 3), res$label[1] == "replicated", res$label[2] != "replicated", res$label[3] == "not replicated",
          res$label[4] == "replicated", res$label[5] == "opposite sign")
cat("K3 TESTS PASSED\n")
