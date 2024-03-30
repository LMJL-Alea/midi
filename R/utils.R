bvalue <- function(small_delta, big_delta, G) {
  gyromagnetic_ratio()^2 * small_delta^2 * G^2 * (big_delta - small_delta / 3)
}

besselJ_derivative <- function(x, n) {
  ifelse(
    n == 0,
    -besselJ(x, 1),
    (besselJ(x, n - 1) - besselJ(x, n + 1)) / 2
  )
}

# this is 2 exp(-x) M(x, 0.5, 1.5)
kummer <- Vectorize(function(x) {
  integrate(\(t) exp(-x * (1 - t)) * t^-0.5, 0, 1)$value
})
# this is 2 exp(-x) M'(x, 0.5, 1.5)
kummer_derivative <- Vectorize(function(x) {
  integrate(\(t) exp(-x * (1 - t)) * t^0.5, 0, 1)$value
})

concentration_index <- function(k) {
  kummer_derivative(k) / kummer(k)
}

slice_triangles <- function(z, n, i, j, k, l) {
  list(c(z, j, i), c(i, j, l), c(l, j, k), c(k, n, l))
}

cylinder_mesh <- function(r, xs, ys, zs, h, n_slices = 40) {
  theta <- seq(0, 2 * pi, length.out = n_slices + 1)
  x <- xs + r * cos(theta)
  y <- ys + r * sin(theta)
  z1 <- zs + rep(0, length(x))
  z2 <- (zs + h) * rep(1, length(x))

  n <- n_slices * 2 + 1

  triangles <- list()
  for (s in 1:n_slices) {
    j <- ifelse(s <= n_slices - 1, s + 1, 1)
    k <- ifelse(s <= n_slices - 1, j + n_slices, n_slices + 1)
    l <- s + n_slices
    triangles <- c(triangles, slice_triangles(0, n, s, j, k, l))
  }
  triangles <- do.call(rbind, triangles)

  x_coords <- c(xs, x[-length(x)], x[-length(x)], xs)
  y_coords <- c(ys, y[-length(y)], y[-length(y)], ys)
  z_coords <- c(zs, z1[-length(z1)], z2[-length(z2)], zs + h)
  vertices <- cbind(x_coords, y_coords, z_coords)

  list(
    vertices = vertices, triangles = triangles,
    x = x, y = y, z1 = z1, z2 = z2
  )
}

cylinder_traces <- function(r, xs, ys, zs, h,
                            n_slices = 20L,
                            show_circular_mesh = TRUE,
                            n_circles = 2L,
                            show_linear_mesh = FALSE,
                            surface_kw = list(),
                            line_kw = list()) {
  res <- cylinder_mesh(r, xs, ys, zs, h, n_slices)
  vertices <- res$vertices
  triangles <- res$triangles
  x <- res$x
  y <- res$y

  surface <- rlang::list2(
    type = "mesh3d",
    x = vertices[, 1],
    y = vertices[, 2],
    z = vertices[, 3],
    i = triangles[, 1],
    j = triangles[, 2],
    k = triangles[, 3],
    !!!surface_kw
  )

  traces <- list(surface)

  line_kw$showlegend <- FALSE

  if (show_circular_mesh) {
    zsubs <- seq(zs, zs + h, length.out = n_circles)
    for (zc in zsubs) {
      traces <- c(traces, list(rlang::list2(
        type = "scatter3d",
        x = x, y = y, z = rep(zc, length(x)),
        mode = "lines", !!!line_kw
      )))
    }
  }

  if (show_linear_mesh) {
    for (i in 1:length(x)) {
      traces <- c(traces, list(rlang::list2(
        type = "scatter3d",
        x = c(x[i], x[i]),
        y = c(y[i], y[i]),
        z = c(zs, zs + h),
        mode = "lines",
        !!!line_kw
      )))
    }
  }

  traces
}
