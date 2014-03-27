convertToProp <- function(y, T0 = NA, Ctrl = NA){
  if(is.na(Ctrl)) Ctrl <- max(y, na.rm = TRUE)
  if(is.na(T0))
    return(y/Ctrl)
  else return((y-T0)/(Ctrl-T0))
}
