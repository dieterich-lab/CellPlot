#library(miscset)

cp_within <- function (x, r) {
  x >= min(r) & x <= max(r)
}

cp_which_next_to <- function(x, n) {
  x <- x[1]
  which.min(abs(n-x))
}

cp_colorpal <- function (x, b = range(x), n = 100, cols = c("deepskyblue2", "white", "coral"), centered = TRUE) {
  if (n < length(x))
    stop("invalid size n")
  if (centered) {
    b <- max(abs(b))
    b <- c(-b, b)
  }
  if (n%%2 > 0)
    n <- n - 1
  if (length(cols)%%2 > 1)
    n <- n + 1
  x <- scaler(x, c(1, n), b = b)
  ix <- sapply(x, cp_which_next_to, seq(n))
  pal <- colorRampPalette(cols)(n)
  pal[ix]
}

#x <- c(.3,1,.1)
#.which_next(1,x)
#barplot(x,,col=.color3(x))
#barplot(x,col = .color3(x, c(0,.5)))

cp_cells <- function (hi, pi, ci, horiz, w, h, p, u, cols, cellrange, centered) {
  nci <- length(ci)
  sh <- seq(0, hi, length.out = nci)
  sw <- c(pi - w, pi + w)
  cols <- cp_color3(ci, cols=cols, centered=centered)
  if (horiz) {
    sh <- scaler(sh, u[1:2], b = c(0, max(h)))
    sw <- scaler(sw, u[3:4], b = c(0 - w, max(p) + w))
    rect(sh[-nci], sw[1], sh[-1], sw[2], col = cols)
  } else {
    sh <- scaler(sh, u[3:4])
    sw <- scaler(sw, u[1:2], b = c(0, max(h)+1))
    rect(sw[1], sh[-nci], sw[2], sh[-1], col = cols)
  }
}

