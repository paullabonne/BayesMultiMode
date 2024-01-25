#' Print method for \code{BayesMode} objects
#' 
#' @param x An object of class \code{BayesMode}.
#' @param max_length maximum number of elements (for vector) or rows (for matrices) to show.
#' @param max_width maximum number of columns to show (for matrices).
#' @param ... Not used.
#' 
#' @importFrom utils head
#' 
#' @export
print.BayesMode <- function(x, max_length = 5, max_width = 5, ...) {
  print_list(x, max_length, max_width)
}


#' Print method for \code{Mode} objects
#' 
#' @param x An object of class \code{Mode}.
#' @param max_length maximum number of elements (for vector) or rows (for matrices) to show.
#' @param max_width maximum number of columns to show (for matrices).
#' @param ... Not used.
#' 
#' @export
print.Mode <- function(x, max_length = 5, max_width = 5, ...) {
  print_list(x, max_length, max_width)
}


#' Print method for \code{BayesMixture} objects
#' 
#' @param x An object of class \code{BayesMixture}.
#' @param max_length maximum number of elements (for vector) or rows (for matrices) to show.
#' @param max_width maximum number of columns to show (for matrices).
#' @param ... Not used.
#' 
#' @export
print.BayesMixture <- function(x, max_length = 5, max_width = 5, ...) {
  print_list(x, max_length, max_width)
}

#' Print method for \code{Mixture} objects
#' 
#' @param x An object of class \code{Mixture}.
#' @param max_length maximum number of elements (for vector) or rows (for matrices) to show.
#' @param max_width maximum number of columns to show (for matrices).
#' @param ... Not used.
#' 
#' @export
print.Mixture <- function(x, max_length = 5, max_width = 5, ...) {
  print_list(x, max_length, max_width)
}


print_list <- function(x, max_length = 5, max_width = 5) {
  assert_that(max_length >= 1,
              max_width >= 1)
  # Check for data type and print accordingly
  # Print list details
  for (i in 1:length(x)) {
    cat("\n", names(x)[i])
    head_print(x[[i]], max_length, max_width)
  }
}

head_print <- function(x, max_length = 5, max_width = 5) {
  n = as.integer(max_length)
  m = as.integer(max_width)
  
  if (is.vector(x)) {
    
    cat(paste0(" (", class(x), sprintf(" vector, dim %d", length(x)),"):"),"\n")
    if (length(x) <= n) {
      print(x)
    } else {
      print(head(x, n))
      cat(sprintf("... (%d more elements)\n", length(x) - n))
    }
    
  } else if (is.matrix(x)) {
    
    # Print matrix
    rows_to_print <- min(nrow(x), n)
    cols_to_print <- min(ncol(x), m)
    cat(paste0(sprintf(" (matrix, dim %dx%d", nrow(x),ncol(x)),"):"),"\n")
    if (nrow(x) <= n & ncol(x) <= m) {
      print(x)
    } else {
      print(head(x, c(rows_to_print, cols_to_print)))
      if (nrow(x) > n & ncol(x) > m) {
        cat(sprintf("... (%d more rows and %d more columns)\n", nrow(x) - n, ncol(x) - m))  
      }
      if (nrow(x) > n & ncol(x) <= m) {
        cat(sprintf("... (%d more rows)\n", nrow(x) - n))  
      }
      if (nrow(x) <= n & ncol(x) > m) {
        cat(sprintf("... (%d more columns)\n", ncol(x) - m))  
      }
    }
    
  } else if (is.null(x)) {
    
    cat(" (NULL)\n")
    
  } else if (is.function(x)) {
    cat(":\n", head(args(x))[1])
    print(body(x))
  } else {
    # Fallback for other types
    print(x)
    cat(sprintf("Type: %s\n", class(x)))
  }
}
