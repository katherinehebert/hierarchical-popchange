# Plot distributions for conceptual figure

set.seed(22)

df = data.frame(
"stable" = rnorm(100, 0, .1),
"unstable" = rnorm(100, 0, .3)
)

ggplot(data = df) +
  geom_density(aes(x = stable)) +
  coord_cartesian(xlim = c(-1,1))

ggplot(data = df) +
  geom_density(aes(x = unstable))+
  coord_cartesian(xlim = c(-1,1))


curve(dnorm(x), xlim = c(-10,10), ylim = c(0, 0.5))
abline(v = 0)
abline(v = c(-1,1), lty = 2)

curve(dnorm(x,sd = 2.5), xlim = c(-10,10), ylim = c(0, 0.5))
abline(v = 0)
abline(v = c(-2.5,2.5), lty = 2)

curve(dnorm(x,sd = 4), xlim = c(-10,10), ylim = c(0, 0.5))
abline(v = 0)
abline(v = c(-4,4), lty = 2)
