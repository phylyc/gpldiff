library(io)
library(ggplot2)

# Explore model behaviour under different group assignments

lapply(list.files("../R/", "\\.R$", full.names=TRUE), source);

set.seed(5);
# simulate two non-linear data series with missingness
data0 <- rldiff(200);
f <- data0$f;

options(plot=list(height=10))
graphics.off()

fit0 <- gpldiff(data0)
qdraw({plot(fit0, data=data0)})

# g in {-0.25, 0.25}
data1 <- data0;
data1$g <- sign(data1$g) / 4;
fit1 <- gpldiff(data1)
qdraw({plot(fit1, data=data1)})

# g in {-1, 1}
data2 <- data0;
data2$g <- sign(data2$g);
fit2 <- gpldiff(data2)
qdraw({plot(fit2, data=data2)})
plot(fit2, data=data2)

# g is scaled s.t. mean(g) = 0 and var(g) = 1
data3 <- data0;
# g needs to be a vector or else diag(g) will misbehaviour in gpldiff
data3$g <- as.numeric(scale(sign(data3$g)));
fit3 <- gpldiff(data3)
qdraw({plot(fit3, data=data3)})

# g in {0, 1}
data4 <- data0;
data4$g <- ifelse(data4$g > 0, 1, 0);
fit4 <- gpldiff(data4)
qdraw({plot(fit4, data=data4)})

