match_from_dots <- function(dots, fn){
  arg <- match(names(formals(fn)), names(dots))
  dots[arg[!is.na(arg)]]
}

wrap <- function(...){
  dots <- list(...)

  checkmate::assert_named(dots)

  list(
    fn1 = do.call("fn1", match_from_dots(dots, fn1)),
    fn2 = do.call("fn2", match_from_dots(dots, fn2))
  )
}
