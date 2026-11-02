source("lib_replication.R")
tf <- function(r, df) r * sqrt(df) / sqrt(1 - r^2)   # r -> t
stopifnot(abs(partial_r(tf(0.4, 48), 48) - 0.4) < 1e-12)
stopifnot(identical(dir_sign(c("+","positive","-","negative","higher","lower","x")), c(1,1,-1,-1,1,-1,NA)))
# cohorts: nelore n=52 df=48 ; charolais n=30 df=26
mk <- function(id, rn, rc) rbind(
  if (!is.na(rn)) data.frame(id=id, cohort="nelore", W=tf(rn,48), df=48, n=52),
  if (!is.na(rc)) data.frame(id=id, cohort="charolais", W=tf(rc,26), df=26, n=30))
# family of 6 (k=2 type): 
tg <- data.frame(id=paste0("L",1:6), dir=c(1,1,1,-1,1,1), testable=c(TRUE,TRUE,TRUE,TRUE,TRUE,FALSE))
est <- rbind(mk("L1", .70, .75),   # both strong, concordant -> replicated
             mk("L2", .60, .35),   # pooled sig, charolais alone not -> single-cohort support
             mk("L3", -.70,-.75),  # both strong, opposite to published(+) -> opposite sign
             mk("L4", .05, .02),   # null -> not replicated
             mk("L5", .70, NA))    # only nelore -> k=1 single-cohort support
res <- classify_set(tg, est)
print(res[, c("id","k","r","p","q_bh","family_size","loo_all_sig","loo_nominal","label")], digits=3)
exp <- c("replicated","single-cohort support","opposite sign","not replicated","single-cohort support","untestable")
stopifnot(identical(res$label, exp), all(res$family_size == 6), res$q_bh[6] == 1)
# k=1 exact t p
x <- pool_one(2.5, 26, 30); stopifnot(abs(x$p - 2*pt(-2.5,26)) < 1e-12)
# BH sizes 64 and 6
tg64 <- data.frame(id=paste0("N",1:64), dir=rep(1,64), testable=c(rep(TRUE,35), rep(FALSE,29)))
set.seed(1); e64 <- do.call(rbind, lapply(1:35, function(i) data.frame(id=paste0("N",i), cohort="charolais",
        W=ifelse(i<=3, 4.5, rnorm(1)), df=26, n=30)))
r64 <- classify_set(tg64, e64); stopifnot(all(r64$family_size==64), sum(r64$label=="untestable")==29)
print(table(r64$label))
# BH adjustment uses family size 64 not 35: q for N1 = p*64/rank
p <- r64$p; stopifnot(all.equal(r64$q_bh, p.adjust(p,"BH")))
cat("ALL TESTS PASSED\n")
