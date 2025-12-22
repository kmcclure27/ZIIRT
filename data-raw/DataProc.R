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
