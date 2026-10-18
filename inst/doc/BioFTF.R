## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>",
                      fig.width = 6.5, fig.height = 4.5)

## -----------------------------------------------------------------------------
library(BioFTF)

x <- matrix(c(54, 20, 17,  6,  3,
              25, 25, 25, 25,  0,
              35, 35, 30,  0,  0,
              70, 10, 10,  5,  5),
            nrow = 4, byrow = TRUE)
rownames(x) <- c("A", "B", "C", "D")

hill_profile(x, domain = c(0, 1, 2, 10, 100, Inf))

## -----------------------------------------------------------------------------
hill_profile_plot(x, from = 0, to = 10, n = 501)

## -----------------------------------------------------------------------------
hill_profile(x, domain = c(-1, 0), allow_negative = TRUE)

## -----------------------------------------------------------------------------
hill_first_derivative(x, domain = c(0, 1, 2))
hill_second_derivative(x, domain = c(0, 1, 2))
hill_derivative_plot(x, order = 1)
hill_curvature_plot(x)

## -----------------------------------------------------------------------------
hill_area(x)
hill_cum_plot(x)
hill_arc(x)
hill_ranking(x)

## -----------------------------------------------------------------------------
hill_ranking(x, from = 0, to = 0.5)[, c("hill_area", "area_rank")]

## -----------------------------------------------------------------------------
years <- 2018:2027
z <- matrix(c(40, 35, 15,  5,  5,  0,
              38, 34, 16,  6,  4,  2,
              35, 30, 18,  8,  5,  4,
              30, 28, 20, 10,  7,  5,
              25, 26, 21, 12,  9,  7,
              20, 24, 22, 14, 11,  9,
              16, 21, 23, 16, 13, 11,
              12, 18, 24, 18, 15, 13,
               9, 15, 24, 20, 17, 15,
               6, 12, 23, 22, 19, 18),
            nrow = length(years), byrow = TRUE)
rownames(z) <- years

hill_total_change(z)

## ----fig.height = 5.5---------------------------------------------------------
hill_surface_plot(z, from = 0, to = 10, n = 121, plot = "static")

## ----eval = FALSE-------------------------------------------------------------
#  hill_surface_plot(z, from = 0, to = 10, n = 121)

## ----fig.width = 7.5, fig.height = 4------------------------------------------
w <- z
w[, 1] <- z[, 1] + 15
w[, 6] <- z[, 6] + 4

hill_surface_compare(z, w, from = 0, to = 10, n = 121,
                     labels = c("Area 1", "Area 2"), plot = "static")

## ----fig.height = 5.5---------------------------------------------------------
hill_surface_compare(z, w, from = 0, to = 10, n = 121,
                     labels = c("Area 1", "Area 2"), plot = "static",
                     mode = "difference")

## ----eval = FALSE-------------------------------------------------------------
#  hill_surface_compare(z, w, from = 0, to = 10, n = 121,
#                       labels = c("Area 1", "Area 2"))
#  hill_surface_compare(z, w, from = 0, to = 10, n = 121,
#                       labels = c("Area 1", "Area 2"), mode = "difference")

