requiredPackages = c('car','doSNOW','doParallel',
                     'foreach','MixedTS','moments',
                     'parallel','stabledist','tcltk','xtable')
for(p in requiredPackages){
  if(!require(p,character.only = TRUE)) install.packages(p)
  library(p,character.only = TRUE)
}







