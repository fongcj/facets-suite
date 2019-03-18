###########################
#### nguyenb@mskcc.org ####
###########################

# first calculate the size of each intergral copy
get_intergral_size = function(fit) {
  cncf = fit$cncf
  cncf$size = cncf$end - cncf$start
  tt = table(cncf$tcn.em)
  if(length(tt) == 1) {
    integral_size = sum(cncf$size, na.rm = T)
    names(integral_size) = names(tt)
  } else {
    tt2 = names(tt)
    integral_size=vector()
    for(y in 1:length(tt2)){
      integral_size[y] = sum(cncf$size[which(cncf$tcn.em == tt2[y])])
    }
    names(integral_size) = tt2
  }
  names(integral_size) = as.integer(names(integral_size))
  return(integral_size)
}

# two function to calculate FGA and WGD
get_FGA = function(integral_size) {
  1 - integral_size[which.max(integral_size)] / sum(integral_size, na.rm = T)
}

is_WGD = function(integral_size) {
  cut_off = as.integer(names(integral_size)) > 2
  sum(integral_size[cut_off], na.rm = T)/(sum(integral_size, na.rm = T)) > .5
}


