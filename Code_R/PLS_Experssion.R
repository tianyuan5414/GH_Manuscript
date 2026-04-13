#' Declare a function to generate a experrsion function based given component number
#' and given pls model
  #' @param compI single numeric value, indicating the number of component
  #' @param plsModel the pls model output by pls::plsr
lm_eqn2 <- function(compI, plsModel){
  m <- lm(plsModel$model$Y_train ~ plsModel$fitted.values[,,compI]);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2,
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3)))
  as.character(as.expression(eq));
}