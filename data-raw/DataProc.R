## code to prepare `DATASET` dataset goes here
#set.seed(307)
MZI = generate_MZIGRM_data(N = 500, J = 10, K = 4,rho = .3,
                           seed_person_params = 307,
                           seed_item_params = 308,
                           seed_response = 309)
MZI_dat = data.frame(MZI$item_resps,
                     MZI$person_params)
colnames(MZI_dat) = c(paste0("item",1:10),
                      "Presence","Severity")



usethis::use_data(MZI_dat, overwrite = TRUE)

load("~/Research/Data/CAT-SRP Data/STB.RData")
usethis::use_data(STB,overwrite=TRUE)

#Generate MZI proportion ZI table
bg = expand.grid(b0=seq(-6,6,by=.05), #grid of item intercepts
                 b1=seq(-6,6,by=.05))
bg = as.matrix(bg)
pZI = vector(length=nrow(bg))
for(i in 1:nrow(bg)){
  pZI[i] = propZI(b0=bg[i,1],b1=bg[i,2])
}

#Large table of expected proportion of zero inflation due to theta0 (assumes BVN(0,1))
propZI_table = data.frame(propZI=pZI,b0=bg[,1],b1=bg[,2])
row.names(propZI_table) <- NULL
usethis::use_data(propZI_table,internal=TRUE,compress="xz",overwrite=T)
