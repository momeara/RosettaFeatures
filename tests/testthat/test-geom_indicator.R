context("geom_indicator")

library(ggplot2)

test_that("geom_indicator can be used in a plot without error", {

    data <- rbind(
        data.frame(x= c(rnorm(200)+5,rnorm(170)+2),strata=factor(1),slice =factor(1)),
        data.frame(x= c(rnorm(500)+5,rnorm(32 )+2),strata=factor(2),slice =factor(1)),
        data.frame(x= c(rnorm(356)+5,rnorm(120)+2),strata=factor(2),slice =factor(2)))
    p <- ggplot(data=data, aes(x=x, colour=strata, indicator=strata))
    p <- p + geom_density()
    p <- p + geom_indicator()
    p <- p + facet_wrap( ~ slice )
    p
    
})
                
