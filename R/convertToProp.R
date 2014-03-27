convertToProp <- function(y, T0=min(y, na.rm=TRUE), Ctrl=max(y, na.rm=TRUE)){
  # If nothing specified: T0 is min(y), Ctrl is max(y)
  if(is.na(T0) | is.na(Ctrl))
    warning(call.=FALSE, "Please, provide non-NA T0 and Ctrl values.")
  return((y-T0)/(Ctrl-T0))
}
