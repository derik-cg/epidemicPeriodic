parameters <- c(s = 10, r = 28, b = 8/3)
state <- c(X = 0, Y = 1, Z = 1)

Lorenz <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    dX <- s * (Y - X)
    dY <- X * (r - Z) - Y
    dZ <- X * Y - b * Z
    list(c(dX, dY, dZ))
  })
}

times <- seq(0, 500, by = 0.01)
library(deSolve)
out <- ode(y = state, times = times, func = Lorenz, parms = parameters)

par(oma = c(0, 0, 3, 0))
plot(out, xlab = "time", ylab = "-")
plot(out[, "Y"], out[, "Z"], pch = ".", type = "l")
mtext(outer = TRUE, side = 3, "Lorenz model", cex = 1.5)