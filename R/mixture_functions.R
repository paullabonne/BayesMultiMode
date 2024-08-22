#' @importFrom sn dst
#' @importFrom sn dsn
#' @importFrom stats dnorm
#'

# Mixture of pdf_func
#' @keywords internal
pdf_func_mix <- function(x, pars, pdf_func) {
  pdf_func <- match.fun(pdf_func) # solves NOTE "pdf_func is undefined"

  mixture <- 0
  for (i in 1:nrow(pars)) {
    mixture <- mixture + pars[i, "eta"] * pdf_func(x, pars[i, ])
  }
  return(mixture)
}

#' @keywords internal
test_and_export <- function(p, pdf_func, dist, pars_names, dist_type, loc) {
  par_type <- deparse(substitute(p))
  par_type <- str_extract(par_type, "[a-z]+")
  list_func <- list()

  K <- length(p) / length(pars_names)
  assert_that(K %% 1 == 0,
    msg = paste("All variables in", par_type, "must have the same number of components.")
  )

  assert_that("eta" %in% pars_names,
    msg = paste(
      par_type,
      "should include a variable named eta representing mixture proportions."
    )
  )

  pars_mat <- vec_to_mat(p, pars_names)
  pars_mat <- na.omit(pars_mat)

  assert_that(round(sum(pars_mat[, "eta"], na.rm = T), 2) == 1,
    msg = "The mixture proportions, eta, should sum to one."
  )


  if (!is.null(pdf_func)) {
    assert_that(is.function(pdf_func),
      msg = "pdf_func must be a function"
    )
    func_test <- try(pdf_func(1, pars_mat[1, ]), silent = T)
    assert_that(!("try-error" %in% class(func_test)),
      msg = paste0(
        "running pdf_func returns an error; ",
        "\n pdf_func should have two arguments :",
        "\n the first argument represents the observation where the pdf is evaluated;",
        "\n the second argument is a named vector representing the function parameters;",
        "\n for instance: \n pdf_func <- function(x, pars) dnorm(x, pars['mu'], pars['sigma'])"
      )
    )
    assert_that(!is.na(pdf_func(1, pars_mat[1, ])),
      msg = "running pdf_func for the first component returns NA"
    )
    assert_that(!is.na(dist_type),
      msg = "dist_type must be provided when argument pdf_func is used; \n i.e. 'continuous' or 'discrete'"
    )
    assert_that(dist_type %in% c("continuous", "discrete"),
      msg = "dist_type should be either 'discrete' or 'continuous'"
    )
  }

  msg_0 <- paste0("variable names in ", par_type, " should be ")

  if (!is.na(dist)) {
    if (dist == "poisson") {
      assert_that(sum(pars_names %in% c("eta", "lambda")) == 2,
        msg = paste0(msg_0, "eta and lambda when dist = poisson")
      )
      pdf_func <- function(x, pars) dpois(x, pars["lambda"])
    }

    if (dist == "shifted_poisson") {
      assert_that(sum(pars_names %in% c("eta", "kappa", "lambda")) == 3,
        msg = paste0(msg_0, "eta and lambda when dist = shifted_poisson")
      )
      pdf_func <- function(x, pars) dpois(x - pars["kappa"], pars["lambda"])
    }

    if (dist == "normal") {
      assert_that(sum(pars_names %in% c("eta", "mu", "sigma")) == 3,
        msg = paste0(msg_0, "eta, mu and sigma when dist = normal")
      )
      pdf_func <- function(x, pars) dnorm(x, pars["mu"], pars["sigma"])
    }

    if (dist == "skew_normal") {
      assert_that(sum(pars_names %in% c("eta", "xi", "omega", "alpha")) == 4,
        msg = paste0(msg_0, "eta, xi, omega and alpha when dist = skew_normal")
      )
      pdf_func <- function(x, pars) dsn(x, pars["xi"], pars["omega"], pars["alpha"])
      loc <- "xi"
    }

    if (dist %in% c("normal", "skew_normal")) {
      dist_type <- "continuous"
    } else if (dist %in% c("poisson", "shifted_poisson")) {
      dist_type <- "discrete"
    } else {
      stop("Unsupported distribution; dist should be either normal, skew_normal, poisson or shifted_poisson")
    }
  }

  if ((is.na(dist) | !(dist %in% c("normal", "skew_normal"))) & dist_type != "discrete") {
    assert_that(is.string(loc) && !is.na(loc),
      msg = paste0(
        "loc argument must be given when using a continuous distribution other than the normal distribution;",
        "\n loc should be the location parameter of pdf_func"
      )
    )
    assert_that(loc %in% pars_names[pars_names != "eta"],
      msg = paste0("loc must a parameter included in ", par_type, " other than eta")
    )
  }

  list_func$dist_type <- dist_type
  list_func$pdf_func <- pdf_func
  list_func$loc <- loc
  list_func$K <- K

  return(list_func)
}

#' @keywords internal
vec_to_mat <- function(pars, pars_names) {
  pars_mat <- c()

  for (i in 1:length(pars_names)) {
    pars_mat <- cbind(pars_mat, pars[grep(pars_names[i], names(pars))])
  }

  pars_mat <- matrix(pars_mat,
    ncol = length(pars_names),
    dimnames = list(NULL, pars_names)
  )
  # colnames(pars_mat) <- pars_names

  return(pars_mat)
}
