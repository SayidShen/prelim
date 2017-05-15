ggplot(diamonds, aes(carat, price)) + geom_point()
ggplot() + geom_point(aes(carat, price), diamonds)
ggplot(diamonds) + geom_point(aes(carat, price))
ggplot(diamonds, aes(carat)) + geom_point(aes(y = price))


plot(price~carat, data=diamonds)
qplot(carat, price, data=diamonds)


ggplot(diamonds, aes(carat, price)) + geom_point(aes(colour = 'navy'))

ggplot(diamonds, aes(carat, price)) + geom_point(colour = 'navy')
ggplot(diamonds, aes(carat, price)) + geom_point(alpha = 0.2)
ggplot(diamonds, aes(carat, price)) + geom_point(size = 0.2)
ggplot(diamonds, aes(carat, price)) + geom_point(shape = 1)


ggplot(diamonds, aes(carat, price)) + geom_point(size = 0.2) + geom_smooth()

ggplot(diamonds, aes(log10(carat), log10(price))) + geom_point(size = 0.25) + geom_smooth()


ggplot(diamonds, aes(log10(carat), log10(price))) + geom_point(size = 0.5)
ggplot(diamonds, aes(carat, price, log="xy")) + geom_point(size = 0.5)
ggplot(diamonds, aes(carat, price)) + geom_point(size = 0.5) + coord_trans(x = "log10", y = "log10")
ggplot(diamonds, aes(carat, price)) + geom_point(size = 0.5) + coord_trans(x = "log10") + coord_trans(y = "log10")


# color by a continuous variable
ggplot(diamonds, aes(carat, price, colour=depth)) + geom_point() + coord_trans(x = "log10", y = "log10")
# color by a factor variable
ggplot(diamonds, aes(carat, price, colour=color)) + geom_point() + coord_trans(x = "log10", y = "log10")
# color and shape by two different variables
ggplot(diamonds, aes(carat, price, shape=cut, colour=color)) + geom_point() + coord_trans(x = "log10", y = "log10")

# faceting
ggplot(diamonds, aes(carat, price)) + geom_point() + facet_wrap(~color, ncol=4)

ggplot(diamonds, aes(carat, price)) + geom_point() + facet_grid(~color, labeller=label_both)


# boxplot
ggplot(diamonds, aes(1, depth)) + geom_boxplot()

ggplot(diamonds, aes(cut, depth)) + geom_boxplot()


# histogram
ggplot(diamonds, aes(depth)) + geom_histogram()

ggplot(diamonds, aes(depth)) + geom_histogram(binwidth=0.2) + xlim(56, 67)

ggplot(diamonds, aes(depth)) + geom_histogram(binwidth = 0.2) + facet_wrap(~cut) + xlim(56, 67)
ggplot(diamonds, aes(depth, fill=cut)) + geom_histogram(binwidth=0.2) + xlim(56,67)


