###########################
#### nguyenb@mskcc.org ####
###########################
# Tue Apr 16 12:41:16 2019 ------------------------------

# function to calculate the size of each intergral copy
get_integral_size = function(fit, method = 'em') {
  cncf = fit$cncf
  cncf$size = cncf$end - cncf$start
  if(method == 'cncf') {
    tt = table(cncf$tcn)
    }
  else {
    tt = table(cncf$tcn.em)
    }
  if(length(tt) == 1) {
    integral_size = sum(cncf$size, na.rm = T)
    names(integral_size) = names(tt)
  } else {
    tt2 = names(tt)
    integral_size=vector()
    for(y in 1:length(tt2)){
      if(method == 'cncf') {
        integral_size[y] = sum(cncf$size[which(cncf$tcn == tt2[y])])
      } else {
        integral_size[y] = sum(cncf$size[which(cncf$tcn.em == tt2[y])])
      }
    }
    names(integral_size) = tt2
  }
  names(integral_size) = as.integer(names(integral_size))
  return(integral_size)
}

# main function to calculate FGA / GAIN / LOSS / LOH / WGD
get_FGA_BN = function(fit, out, method = 'em', include_loh = T) {
  cncf = fit$cncf
  cncf$size = cncf$end - cncf$start
  integral_size = get_integral_size(fit, method)
  cut_off = as.numeric(names(integral_size)) > 2
  is_WGD = sum(integral_size[cut_off], na.rm = T)/(sum(integral_size, na.rm = T)) > .5
  major_cn = as.numeric(names(which.max(integral_size)))
  gain = sum(integral_size[as.numeric(names(integral_size)) > major_cn], na.rm = T)/sum(integral_size, na.rm = T)
  loss = sum(integral_size[as.numeric(names(integral_size)) < major_cn], na.rm = T)/sum(integral_size, na.rm = T)
  if(method == 'cncf') {
    LOH = sum(cncf$size[which(cncf$tcn == major_cn & cncf$lcn == 0)], na.rm = T)/sum(integral_size, na.rm = T)
  } else {
    LOH = sum(cncf$size[which(cncf$tcn.em == major_cn & cncf$lcn.em == 0)], na.rm = T)/sum(integral_size, na.rm = T)
  }
  if(include_loh) {
    FGA = gain + loss + LOH
  } else {
    FGA = gain + loss
  }
  sample_name = as.character(out$IGV[1,1])
  if(is.null(fit$emflags)) { fit$emflags = '' }
  ouput = data.frame('PURITY' = fit$purity,'PLOIDY' = fit$ploidy,'FGA' = FGA, 'GAIN' = gain, 'LOSS'= loss, 'LOH' = LOH, 'IS_WGD' = is_WGD, 'EM_FLAGS' = fit$emflags, row.names = sample_name)
  return(ouput)
}

get_FChrA = function(fit, out, method = 'em') {
  cncf = fit$cncf
  cncf$size = cncf$end - cncf$start
  integral_size = get_integral_size(fit, method)
  major_cn = as.numeric(names(which.max(integral_size)))
  FChrA = vector(length = length(unique(cncf$chrom)))
  if(method == 'cncf') {
    for(i in 1:length(unique(cncf$chrom))) {
      cncf_chr = cncf[cncf$chrom == unique(cncf$chrom)[i],]
      FChrA[i] = sum(cncf_chr$size[which(cncf_chr$tcn != major_cn)], na.rm = T)/sum(cncf_chr$size, na.rm = T)
    }} else {
      for(i in 1:length(unique(cncf$chrom))) {
        cncf_chr = cncf[cncf$chrom == unique(cncf$chrom)[i],]
        FChrA[i] = sum(cncf_chr$size[which(cncf_chr$tcn.em != major_cn)], na.rm = T)/sum(cncf_chr$size, na.rm = T)
    }
    }
  sample_name = as.character(out$IGV[1,1])
  
  ouput = data.frame(t(FChrA), row.names = sample_name)
  colnames(ouput) = paste0('Chr', unique(cncf$chrom))
  return(ouput)
}


# function to transfrom in calls
get_Calls = function(fit, out, method = 'em') {
  cncf = fit$cncf
  cncf$size = cncf$end - cncf$start
  integral_size = get_integral_size(fit, method)
  major_cn = as.numeric(names(which.max(integral_size)))
  if(method == 'cncf') { 
    cncf$seg.intergral = NA
    cncf$seg.intergral[which(cncf$tcn.em == major_cn)] = 0
    cncf$seg.intergral[which(cncf$tcn.em < major_cn)] = -1
    cncf$seg.intergral[which(cncf$tcn.em > major_cn)] = 1
    cncf$seg.intergral[which(cncf$tcn.em < major_cn/2)] = -2
    cncf$seg.intergral[which(cncf$tcn.em > major_cn*2)] = 2
  } else {
    cncf$seg.intergral = NA
    cncf$seg.intergral[which(cncf$tcn == major_cn)] = 0
    cncf$seg.intergral[which(cncf$tcn < major_cn)] = -1
    cncf$seg.intergral[which(cncf$tcn > major_cn)] = 1
    cncf$seg.intergral[which(cncf$tcn < major_cn/2)] = -2
    cncf$seg.intergral[which(cncf$tcn > major_cn*2)] = 2
  }
  seg_file = out$IGV
  seg_file$seg.mean = cncf$seg.intergral
  return(seg_file)
}

