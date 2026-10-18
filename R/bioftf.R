## ---------------------------------------------------------------------------
## 1. Input checks and grids
## ---------------------------------------------------------------------------

.check_abundance_matrix <- function(x) {
  if (is.data.frame(x)) {
    if (!all(vapply(x, is.numeric, logical(1L)))) {
      stop("'x' must contain numeric columns only.", call. = FALSE)
    }
    x <- as.matrix(x)
  }
  if (is.numeric(x) && is.null(dim(x))) {
    x <- matrix(x, nrow = 1L, dimnames = list(NULL, names(x)))
  }
  if (!is.matrix(x) || !is.numeric(x)) {
    stop("'x' must be a numeric matrix or data frame.", call. = FALSE)
  }
  storage.mode(x) <- "double"
  if (anyNA(x)) {
    stop("'x' must not contain missing values.", call. = FALSE)
  }
  if (any(!is.finite(x))) {
    stop("'x' must contain only finite values.", call. = FALSE)
  }
  if (any(x < 0)) {
    stop("'x' must contain non-negative abundances.", call. = FALSE)
  }
  if (any(rowSums(x) <= 0)) {
    stop("Each community must have a positive total abundance.", call. = FALSE)
  }
  if (is.null(rownames(x))) {
    rownames(x) <- paste("community", seq_len(nrow(x)), sep = ".")
  }
  if (is.null(colnames(x))) {
    colnames(x) <- paste("species", seq_len(ncol(x)), sep = ".")
  }
  x
}

.relative_abundance <- function(x) {
  x <- .check_abundance_matrix(x)
  x / rowSums(x)
}

.hill_domain <- function(domain = NULL, from = 0, to = 2, n = 101,
                         allow_negative = FALSE) {
  if (is.null(domain)) {
    if (length(n) != 1L || !is.numeric(n) || !is.finite(n) || n < 3) {
      stop("'n' must be a single number greater than or equal to 3.",
           call. = FALSE)
    }
    if (length(from) != 1L || length(to) != 1L ||
        !is.finite(from) || !is.finite(to)) {
      stop("'from' and 'to' must be finite numbers. Use 'domain' to ",
           "evaluate the profile at Inf.", call. = FALSE)
    }
    if (to <= from) {
      stop("'to' must be greater than 'from'.", call. = FALSE)
    }
    domain <- seq(from, to, length.out = as.integer(round(n)))
  }
  if (!is.numeric(domain) || length(domain) < 1L || anyNA(domain)) {
    stop("'domain' must be a numeric vector without missing values.",
         call. = FALSE)
  }
  domain <- as.numeric(domain)
  if (any(domain < 0) && !allow_negative) {
    stop("Negative orders are rejected by default: for q < 0 the Hill ",
         "number exceeds species richness and is driven by the rarest ",
         "species. Set 'allow_negative = TRUE' to compute it anyway.",
         call. = FALSE)
  }
  domain
}

## Grid required by tools that integrate or differentiate along q.
.finite_grid <- function(q, what) {
  if (any(!is.finite(q))) {
    stop(what, " require finite orders.", call. = FALSE)
  }
  if (length(q) < 2L || any(diff(q) <= 0)) {
    stop(what, " require at least two strictly increasing orders.",
         call. = FALSE)
  }
  invisible(q)
}

.q_labels <- function(q) format(q, digits = 10, trim = TRUE)

.as_profile <- function(values, q, communities) {
  out <- matrix(values, nrow = length(q), ncol = length(communities),
                dimnames = list(.q_labels(q), communities))
  attr(out, "q") <- q
  out
}

.trapz <- function(x, y) {
  sum(diff(x) * (y[-length(y)] + y[-1L]) / 2)
}

.cumtrapz <- function(x, y) {
  c(0, cumsum(diff(x) * (y[-length(y)] + y[-1L]) / 2))
}

## ---------------------------------------------------------------------------
## 2. Hill number and its analytical derivatives
## ---------------------------------------------------------------------------

## Returns c(D, D', D'') at a finite order q for a vector of relative
## abundances. With L(q) = log(sum(p^q)) and g(q) = L(q) / (1 - q), the Hill
## number is D = exp(g). L' and L'' are the mean and the variance of log(p)
## under the weights p^q / sum(p^q). For |q - 1| < 1e-3 the closed form is
## replaced by its Taylor expansion, whose coefficients are the cumulants of
## log(p) under p.
.hill_core <- function(p, q) {
  p <- p[p > 0]
  lp <- log(p)
  e <- q - 1
  if (abs(e) < 1e-3) {
    m1 <- sum(p * lp)
    d <- lp - m1
    mu2 <- sum(p * d^2)
    mu3 <- sum(p * d^3)
    mu4 <- sum(p * d^4)
    mu5 <- sum(p * d^5)
    mu6 <- sum(p * d^6)
    k4 <- mu4 - 3 * mu2^2
    k5 <- mu5 - 10 * mu3 * mu2
    k6 <- mu6 - 15 * mu4 * mu2 - 10 * mu3^2 + 30 * mu2^3
    g <- -(m1 + mu2 * e / 2 + mu3 * e^2 / 6 + k4 * e^3 / 24 +
             k5 * e^4 / 120 + k6 * e^5 / 720)
    g1 <- -(mu2 / 2 + mu3 * e / 3 + k4 * e^2 / 8 + k5 * e^3 / 30 +
              k6 * e^4 / 144)
    g2 <- -(mu3 / 3 + k4 * e / 4 + k5 * e^2 / 10 + k6 * e^3 / 36)
  } else {
    a <- q * lp
    amax <- max(a)
    w <- exp(a - amax)
    sw <- sum(w)
    w <- w / sw
    ## log1p/expm1 keep the relative accuracy of L when q is close to 1.
    L <- if (abs(e) < 0.5) log1p(sum(p * expm1(e * lp))) else amax + log(sw)
    m <- sum(w * lp)
    v <- sum(w * (lp - m)^2)
    u <- 1 - q
    g <- L / u
    g1 <- m / u + L / u^2
    g2 <- v / u + 2 * m / u^2 + 2 * L / u^3
  }
  D <- exp(g)
  c(D, D * g1, D * (g1^2 + g2))
}

.hill_number <- function(p, q) {
  p <- p[p > 0]
  if (is.infinite(q)) {
    return(if (q > 0) 1 / max(p) else 1 / min(p))
  }
  if (q == 0) {
    return(length(p))
  }
  .hill_core(p, q)[1L]
}

## which: 1 = profile, 2 = first derivative, 3 = second derivative.
.hill_matrix <- function(x, q, which = 1L) {
  p <- .relative_abundance(x)
  values <- vapply(seq_len(nrow(p)), function(i) {
    vapply(q, function(value) {
      if (which == 1L) .hill_number(p[i, ], value)
      else .hill_core(p[i, ], value)[which]
    }, numeric(1L))
  }, numeric(length(q)))
  .as_profile(values, q, rownames(p))
}

## ---------------------------------------------------------------------------
## 3. Data summaries
## ---------------------------------------------------------------------------

datirel <- function(x) {
  .relative_abundance(x)
}

summary_species <- function(x) {
  x <- .check_abundance_matrix(x)
  data.frame(
    richness = rowSums(x > 0),
    total_abundance = rowSums(x),
    min_abundance = apply(x, 1L, min),
    max_abundance = apply(x, 1L, max),
    row.names = rownames(x)
  )
}

summary_species_relative <- function(x) {
  p <- .relative_abundance(x)
  data.frame(
    richness = rowSums(p > 0),
    min_relative_abundance = apply(p, 1L, min),
    max_relative_abundance = apply(p, 1L, max),
    row.names = rownames(p)
  )
}

## ---------------------------------------------------------------------------
## 4. Profiles and functional tools
## ---------------------------------------------------------------------------

hill_profile <- function(x, from = 0, to = 2, n = 101, domain = NULL,
                         allow_negative = FALSE) {
  q <- .hill_domain(domain = domain, from = from, to = to, n = n,
                    allow_negative = allow_negative)
  .hill_matrix(x, q, which = 1L)
}

hill <- hill_profile

hill_derivative <- function(x, order = 1, from = 0, to = 2, n = 101,
                            domain = NULL, allow_negative = FALSE) {
  if (length(order) != 1L || !order %in% c(1, 2)) {
    stop("'order' must be 1 or 2.", call. = FALSE)
  }
  q <- .hill_domain(domain = domain, from = from, to = to, n = n,
                    allow_negative = allow_negative)
  if (any(!is.finite(q))) {
    stop("Derivatives require finite orders.", call. = FALSE)
  }
  .hill_matrix(x, q, which = as.integer(order) + 1L)
}

hill_first_derivative <- function(x, from = 0, to = 2, n = 101,
                                  domain = NULL, allow_negative = FALSE) {
  hill_derivative(x, order = 1, from = from, to = to, n = n,
                  domain = domain, allow_negative = allow_negative)
}

hill_second_derivative <- function(x, from = 0, to = 2, n = 101,
                                   domain = NULL, allow_negative = FALSE) {
  hill_derivative(x, order = 2, from = from, to = to, n = n,
                  domain = domain, allow_negative = allow_negative)
}

hill_curvature <- function(x, from = 0, to = 2, n = 101, domain = NULL,
                           allow_negative = FALSE) {
  d1 <- hill_derivative(x, order = 1, from = from, to = to, n = n,
                        domain = domain, allow_negative = allow_negative)
  d2 <- hill_derivative(x, order = 2, from = from, to = to, n = n,
                        domain = domain, allow_negative = allow_negative)
  abs(d2) / (1 + d1^2)^(3 / 2)
}

hill_radius <- function(x, from = 0, to = 2, n = 101, domain = NULL,
                        allow_negative = FALSE) {
  curv <- hill_curvature(x, from = from, to = to, n = n, domain = domain,
                         allow_negative = allow_negative)
  out <- 1 / curv
  out[curv < 1e-8] <- NA_real_
  out
}

.area_vector <- function(y, baseline) {
  q <- attr(y, "q")
  .finite_grid(q, "Profile areas")
  out <- apply(y - baseline, 2L, function(value) .trapz(q, value))
  names(out) <- colnames(y)
  out
}

.arc_vector <- function(y) {
  q <- attr(y, "q")
  .finite_grid(q, "Arc lengths")
  out <- apply(y, 2L, function(value) sum(sqrt(diff(q)^2 + diff(value)^2)))
  names(out) <- colnames(y)
  out
}

.check_baseline <- function(baseline) {
  if (length(baseline) != 1L || !is.numeric(baseline) ||
      !is.finite(baseline)) {
    stop("'baseline' must be a single finite number.", call. = FALSE)
  }
  baseline
}

hill_area <- function(x, from = 0, to = 2, n = 101, domain = NULL,
                      baseline = 0, allow_negative = FALSE) {
  baseline <- .check_baseline(baseline)
  y <- hill_profile(x, from = from, to = to, n = n, domain = domain,
                    allow_negative = allow_negative)
  out <- .area_vector(y, baseline)
  data.frame(hill_area = unname(out), row.names = names(out))
}

hill_cum <- function(x, from = 0, to = 2, n = 101, domain = NULL,
                     baseline = 0, allow_negative = FALSE) {
  baseline <- .check_baseline(baseline)
  y <- hill_profile(x, from = from, to = to, n = n, domain = domain,
                    allow_negative = allow_negative)
  q <- attr(y, "q")
  .finite_grid(q, "Cumulative integrals")
  out <- apply(y - baseline, 2L, function(value) .cumtrapz(q, value))
  .as_profile(out, q, colnames(y))
}

hill_arc <- function(x, from = 0, to = 2, n = 101, domain = NULL,
                     allow_negative = FALSE) {
  y <- hill_profile(x, from = from, to = to, n = n, domain = domain,
                    allow_negative = allow_negative)
  out <- .arc_vector(y)
  data.frame(arc_length = unname(out), row.names = names(out))
}

hill_ranking <- function(x, from = 0, to = 2, n = 101, domain = NULL,
                         baseline = 0, allow_negative = FALSE) {
  baseline <- .check_baseline(baseline)
  y <- hill_profile(x, from = from, to = to, n = n, domain = domain,
                    allow_negative = allow_negative)
  aa <- .area_vector(y, baseline)
  ar <- .arc_vector(y)
  key <- .hill_matrix(x, c(0, 1, 2), which = 1L)
  data.frame(
    richness = key[1L, ],
    shannon = key[2L, ],
    simpson = key[3L, ],
    hill_area = unname(aa),
    area_rank = rank(-aa, ties.method = "average"),
    arc_length = unname(ar),
    arc_rank = rank(ar, ties.method = "average"),
    row.names = names(aa)
  )
}

hill_total_change <- function(x, from = 0, to = 2, n = 101, domain = NULL,
                              baseline = 0, allow_negative = FALSE) {
  baseline <- .check_baseline(baseline)
  y <- hill_profile(x, from = from, to = to, n = n, domain = domain,
                    allow_negative = allow_negative)
  aa <- .area_vector(y, baseline)
  k <- length(aa)
  if (k < 2L) {
    stop("At least two communities or time points are required.",
         call. = FALSE)
  }
  data.frame(
    from = names(aa)[-k],
    to = names(aa)[-1L],
    area_from = unname(aa[-k]),
    area_to = unname(aa[-1L]),
    absolute_change = unname(diff(aa)),
    relative_change = unname(diff(aa) / aa[-k]),
    stringsAsFactors = FALSE
  )
}

## ---------------------------------------------------------------------------
## 5. Profile plots
## ---------------------------------------------------------------------------

.curve_colours <- function(k) {
  if (k == 1L) "black" else grDevices::hcl.colors(k, "Dark 3")
}

.plot_curves <- function(y, xlab, ylab, main, legend = "topright", ...) {
  q <- attr(y, "q")
  if (any(!is.finite(q))) {
    stop("Plots require finite orders.", call. = FALSE)
  }
  k <- ncol(y)
  cols <- .curve_colours(k)
  ltys <- rep_len(1:6, k)
  graphics::matplot(q, y, type = "l", lty = ltys, col = cols, lwd = 2,
                    xlab = xlab, ylab = ylab, main = main, ...)
  if (!is.null(legend)) {
    graphics::legend(legend, legend = colnames(y), lty = ltys, col = cols,
                     lwd = 2, cex = 0.8, bty = "n")
  }
  invisible(y)
}

hill_profile_plot <- function(x, from = 0, to = 2, n = 101, domain = NULL,
                              allow_negative = FALSE, legend = "topright",
                              ...) {
  .plot_curves(
    hill_profile(x, from = from, to = to, n = n, domain = domain,
                 allow_negative = allow_negative),
    xlab = "Order q", ylab = "Hill number",
    main = "Hill-number diversity profiles", legend = legend, ...
  )
}

hill_plot <- hill_profile_plot

hill_derivative_plot <- function(x, order = 1, from = 0, to = 2, n = 101,
                                 domain = NULL, allow_negative = FALSE,
                                 legend = "bottomright", ...) {
  y <- hill_derivative(x, order = order, from = from, to = to, n = n,
                       domain = domain, allow_negative = allow_negative)
  label <- if (order == 1) "First derivative" else "Second derivative"
  .plot_curves(y, xlab = "Order q", ylab = label,
               main = paste(label, "of Hill-number profiles"),
               legend = legend, ...)
}

hill_curvature_plot <- function(x, from = 0, to = 2, n = 101,
                                domain = NULL, allow_negative = FALSE,
                                legend = "topright", ...) {
  .plot_curves(
    hill_curvature(x, from = from, to = to, n = n, domain = domain,
                   allow_negative = allow_negative),
    xlab = "Order q", ylab = "Curvature",
    main = "Curvature of Hill-number profiles", legend = legend, ...
  )
}

hill_radius_plot <- function(x, from = 0, to = 2, n = 101,
                             domain = NULL, allow_negative = FALSE,
                             legend = "topleft", ...) {
  .plot_curves(
    hill_radius(x, from = from, to = to, n = n, domain = domain,
                allow_negative = allow_negative),
    xlab = "Order q", ylab = "Radius of curvature (log scale)",
    main = "Radius of curvature of Hill-number profiles",
    legend = legend, log = "y", ...
  )
}

hill_cum_plot <- function(x, from = 0, to = 2, n = 101, domain = NULL,
                          baseline = 0, allow_negative = FALSE,
                          legend = "topleft", ...) {
  .plot_curves(
    hill_cum(x, from = from, to = to, n = n, domain = domain,
             baseline = baseline, allow_negative = allow_negative),
    xlab = "Order q", ylab = "Cumulative integral",
    main = "Cumulative integral functions of Hill-number profiles",
    legend = legend, ...
  )
}

## ---------------------------------------------------------------------------
## 6. Biodiversity surfaces
## ---------------------------------------------------------------------------

.surface_time <- function(time, labels) {
  if (is.null(time)) {
    time <- suppressWarnings(as.numeric(labels))
    if (anyNA(time) || any(diff(time) <= 0)) {
      time <- seq_along(labels)
    }
  }
  if (!is.numeric(time) || length(time) != length(labels) ||
      any(!is.finite(time))) {
    stop("'time' must be a finite numeric vector with one value for each ",
         "row of the data.", call. = FALSE)
  }
  if (length(time) < 2L || any(diff(time) <= 0)) {
    stop("A surface requires at least two rows and strictly increasing ",
         "'time' values.", call. = FALSE)
  }
  as.numeric(time)
}

hill_surface <- function(x, from = 0, to = 2, n = 101, domain = NULL,
                         time = NULL, allow_negative = FALSE) {
  y <- hill_profile(x, from = from, to = to, n = n, domain = domain,
                    allow_negative = allow_negative)
  q <- attr(y, "q")
  .finite_grid(q, "Surfaces")
  time <- .surface_time(time, colnames(y))
  z <- matrix(as.vector(y), nrow = length(q),
              dimnames = list(.q_labels(q), colnames(y)))
  list(q = q, time = time, z = z)
}

## Static perspective plot with facets coloured by height.
.persp_surface <- function(s, xlab, ylab, zlab, main, theta, phi,
                           palette = "viridis", zlim = range(s$z),
                           diverging = FALSE, ...) {
  nq <- length(s$q)
  nt <- length(s$time)
  facet <- (s$z[-1L, -1L] + s$z[-1L, -nt] + s$z[-nq, -1L] +
              s$z[-nq, -nt]) / 4
  if (diff(zlim) <= 0) {
    zlim <- zlim + c(-0.5, 0.5)
  }
  crange <- if (diverging) c(-1, 1) * max(abs(zlim)) else zlim
  cols <- grDevices::hcl.colors(100L, palette)
  index <- 1L + floor(99 * (facet - crange[1L]) / diff(crange))
  index <- pmin(pmax(index, 1L), 100L)
  pmat <- graphics::persp(
    x = s$q, y = s$time, z = s$z, zlim = zlim, col = cols[index],
    border = NA, theta = theta, phi = phi, ticktype = "detailed",
    cex.axis = 0.8, xlab = paste0("\n", xlab), ylab = paste0("\n", ylab),
    zlab = paste0("\n", zlab), main = main, ...
  )
  for (j in seq_len(nt)) {
    graphics::lines(grDevices::trans3d(s$q, s$time[j], s$z[, j], pmat),
                    col = "grey20", lwd = 0.6)
  }
  for (i in unique(round(seq(1, nq, length.out = min(nq, 11L))))) {
    graphics::lines(grDevices::trans3d(s$q[i], s$time, s$z[i, ], pmat),
                    col = "grey20", lwd = 0.6)
  }
  invisible(pmat)
}

.need_plotly <- function() {
  if (!requireNamespace("plotly", quietly = TRUE)) {
    stop("Interactive surfaces require the 'plotly' package. Install it ",
         "or use plot = \"static\".", call. = FALSE)
  }
}

.scene <- function(xlab, ylab, zlab) {
  list(xaxis = list(title = xlab), yaxis = list(title = ylab),
       zaxis = list(title = zlab))
}

hill_surface_plot <- function(x, from = 0, to = 2, n = 101, domain = NULL,
                              time = NULL,
                              plot = c("interactive", "static"),
                              allow_negative = FALSE,
                              xlab = "Order q", ylab = "Time",
                              zlab = "Hill number",
                              main = "Hill-number surface",
                              theta = 40, phi = 25, ...) {
  plot <- match.arg(plot)
  s <- hill_surface(x, from = from, to = to, n = n, domain = domain,
                    time = time, allow_negative = allow_negative)
  if (plot == "interactive") {
    .need_plotly()
    p <- plotly::plot_ly(x = s$q, y = s$time, z = t(s$z), type = "surface",
                         colorscale = "Viridis", ...)
    return(plotly::layout(p, title = main,
                          scene = .scene(xlab, ylab, zlab)))
  }
  .persp_surface(s, xlab = xlab, ylab = ylab, zlab = zlab, main = main,
                 theta = theta, phi = phi, ...)
  invisible(s)
}

hill_surface_compare <- function(x, y, from = 0, to = 2, n = 101,
                                 domain = NULL, time = NULL,
                                 labels = c("Surface 1", "Surface 2"),
                                 plot = c("interactive", "static"),
                                 mode = c("overlay", "difference"),
                                 allow_negative = FALSE,
                                 main = "Comparison of Hill-number surfaces",
                                 theta = 40, phi = 25, ...) {
  plot <- match.arg(plot)
  mode <- match.arg(mode)
  if (length(labels) != 2L) {
    stop("'labels' must have length 2.", call. = FALSE)
  }
  if (NROW(x) != NROW(y)) {
    stop("'x' and 'y' must have the same number of rows.", call. = FALSE)
  }
  sx <- hill_surface(x, from = from, to = to, n = n, domain = domain,
                     time = time, allow_negative = allow_negative)
  sy <- hill_surface(y, from = from, to = to, n = n, domain = domain,
                     time = sx$time, allow_negative = allow_negative)
  sd <- list(q = sx$q, time = sx$time, z = sy$z - sx$z)
  out <- list(surface1 = sx, surface2 = sy, difference = sd)
  dmain <- paste0(main, ": ", labels[2L], " minus ", labels[1L])
  if (plot == "interactive") {
    .need_plotly()
    if (mode == "difference") {
      lim <- max(abs(sd$z))
      if (lim <= 0) lim <- 1
      p <- plotly::plot_ly(x = sd$q, y = sd$time, z = t(sd$z),
                           type = "surface", colorscale = "RdBu",
                           cmin = -lim, cmax = lim, ...)
      return(plotly::layout(p, title = dmain,
                            scene = .scene("Order q", "Time", "Difference")))
    }
    flat <- function(colour) list(list(0, colour), list(1, colour))
    p <- plotly::plot_ly(...)
    p <- plotly::add_surface(p, x = sx$q, y = sx$time, z = t(sx$z),
                             name = labels[1L], colorscale = flat("#1F77B4"),
                             opacity = 0.75, showscale = FALSE,
                             showlegend = TRUE)
    p <- plotly::add_surface(p, x = sy$q, y = sy$time, z = t(sy$z),
                             name = labels[2L], colorscale = flat("#D95F02"),
                             opacity = 0.75, showscale = FALSE,
                             showlegend = TRUE)
    return(plotly::layout(p, title = main,
                          scene = .scene("Order q", "Time", "Hill number")))
  }
  if (mode == "difference") {
    .persp_surface(sd, xlab = "Order q", ylab = "Time", zlab = "Difference",
                   main = dmain, theta = theta, phi = phi,
                   palette = "Blue-Red 3", diverging = TRUE, ...)
  } else {
    oldpar <- graphics::par(mfrow = c(1, 2), mar = c(2, 1, 3, 1))
    on.exit(graphics::par(oldpar), add = TRUE)
    zlim <- range(sx$z, sy$z)
    .persp_surface(sx, xlab = "Order q", ylab = "Time", zlab = "Hill number",
                   main = labels[1L], theta = theta, phi = phi,
                   zlim = zlim, ...)
    .persp_surface(sy, xlab = "Order q", ylab = "Time", zlab = "Hill number",
                   main = labels[2L], theta = theta, phi = phi,
                   zlim = zlim, ...)
  }
  invisible(out)
}

## ---------------------------------------------------------------------------
## 7. Wrappers
## ---------------------------------------------------------------------------

export_plots <- function(x, file, from = 0, to = 2, n = 101) {
  if (missing(file) || !is.character(file) || length(file) != 1L) {
    stop("'file' must be the path of the PDF file to create.",
         call. = FALSE)
  }
  grDevices::pdf(file)
  on.exit(grDevices::dev.off(), add = TRUE)
  hill_profile_plot(x, from = from, to = to, n = n)
  hill_derivative_plot(x, order = 1, from = from, to = to, n = n)
  hill_derivative_plot(x, order = 2, from = from, to = to, n = n)
  hill_curvature_plot(x, from = from, to = to, n = n)
  hill_radius_plot(x, from = from, to = to, n = n)
  hill_cum_plot(x, from = from, to = to, n = n)
  invisible(file)
}

alltools <- function(x, from = 0, to = 2, n = 101) {
  list(
    relative_abundances = datirel(x),
    species_summary = summary_species(x),
    hill_profile = hill_profile(x, from = from, to = to, n = n),
    first_derivative = hill_first_derivative(x, from = from, to = to, n = n),
    second_derivative = hill_second_derivative(x, from = from, to = to,
                                               n = n),
    curvature = hill_curvature(x, from = from, to = to, n = n),
    radius = hill_radius(x, from = from, to = to, n = n),
    area = hill_area(x, from = from, to = to, n = n),
    cumulative_integral = hill_cum(x, from = from, to = to, n = n),
    arc = hill_arc(x, from = from, to = to, n = n),
    ranking = hill_ranking(x, from = from, to = to, n = n)
  )
}

## ---------------------------------------------------------------------------
## 8. Deprecated names of BioFTF 1.x
## ---------------------------------------------------------------------------

.legacy <- function(old, new) {
  .Deprecated(new, package = "BioFTF", old = old, msg = paste0(
    "'", old, "()' is deprecated; use '", new, "()' instead. ",
    "Since BioFTF 3.0.0 this tool is computed on the Hill-number profile ",
    "(order q in [0, 2] by default) and no longer on the Patil-Taillie ",
    "beta profile (beta in [-1, 1]), so the values differ from those of ",
    "BioFTF 1.x."
  ))
}

beta <- function(x, n = 101, ...) {
  .legacy("beta", "hill_profile")
  hill_profile(x, n = n, ...)
}

beta_plot <- function(x, n = 101, ...) {
  .legacy("beta_plot", "hill_profile_plot")
  hill_profile_plot(x, n = n, ...)
}

betaprime <- function(x, n = 101, ...) {
  .legacy("betaprime", "hill_first_derivative")
  hill_first_derivative(x, n = n, ...)
}

betaprime_plot <- function(x, n = 101, ...) {
  .legacy("betaprime_plot", "hill_derivative_plot")
  hill_derivative_plot(x, order = 1, n = n, ...)
}

betasecond <- function(x, n = 101, ...) {
  .legacy("betasecond", "hill_second_derivative")
  hill_second_derivative(x, n = n, ...)
}

betasecond_plot <- function(x, n = 101, ...) {
  .legacy("betasecond_plot", "hill_derivative_plot")
  hill_derivative_plot(x, order = 2, n = n, ...)
}

curvature <- function(x, n = 101, ...) {
  .legacy("curvature", "hill_curvature")
  hill_curvature(x, n = n, ...)
}

curvature_plot <- function(x, n = 101, ...) {
  .legacy("curvature_plot", "hill_curvature_plot")
  hill_curvature_plot(x, n = n, ...)
}

radius <- function(x, n = 101, ...) {
  .legacy("radius", "hill_radius")
  hill_radius(x, n = n, ...)
}

radius_plot <- function(x, n = 101, ...) {
  .legacy("radius_plot", "hill_radius_plot")
  hill_radius_plot(x, n = n, ...)
}

arc <- function(x, n = 101, ...) {
  .legacy("arc", "hill_arc")
  hill_arc(x, n = n, ...)
}

area <- function(x, n = 101, ...) {
  .legacy("area", "hill_area")
  hill_area(x, n = n, ...)
}

ranking <- function(x, n = 101, ...) {
  .legacy("ranking", "hill_ranking")
  hill_ranking(x, n = n, ...)
}
