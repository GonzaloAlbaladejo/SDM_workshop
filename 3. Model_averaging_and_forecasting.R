# 3. Model Averaging and forecasting





Zli <- li[, 1:(n_speciation_axes + 1)]
f1 <- function(x) rep(x, prb)

Sli <- apply(Zli, 2, f1)
m_cord <- apply(Sli, 2, mean) # Specificity and marginality coordinates
cov <- t(as.matrix(Sli)) %*% as.matrix(Sli)/nrow(Sli)
maha <- data.frame(Mahalanobis.Dist = mahalanobis(Zli, center = m_cord, cov = cov))