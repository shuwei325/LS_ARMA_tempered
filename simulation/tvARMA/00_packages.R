requiredPackages = c('car','doSNOW','doParallel',
                     'foreach','LSTS','MixedTS','moments',
                     'parallel','stabledist','tcltk','tseries','xtable')
for(p in requiredPackages){
  if(!require(p,character.only = TRUE)) install.packages(p)
  library(p,character.only = TRUE)
}





