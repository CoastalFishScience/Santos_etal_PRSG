#' """ Santos et al., Puerto Rico SeaGrant Final Report 
#'     @authors: Valentina Bautista and W.R. James 
#'     date: """6/30/2025


# Load Libraries  ---------------------------------------------------------

library(tidyverse)
library(hypervolume)
library(dplyr)

conflicted::conflict_prefer('select','dpylr')


# Puerto Rico Coral Reef Monitoring Program Dataset  ------------------------------------------------------------
# load trait data 
df_tr = read_csv('data/TraitTable_Imput.csv') 

# load community data 
df_ben <- read.csv('data/benthicdata_imput.csv') |> 
  filter(site %in% c("Carlos Rosario",  "Dakiti", "Cayo Diablo", ####Select for NEMC sites
                     "Palominitos", "Palominos", "Canal Luis Pena")) |> 
  group_by(species, YEAR, site) |> 
  summarize(
    pc = mean(percentcover, na.rm = TRUE),
    pc_sd = sd(percentcover, na.rm = TRUE),
    .groups = 'drop') |> 
  filter(pc > 0)


## Build HV ---------------------------------------------
df_hv = df_ben |> 
  left_join(df_tr, by = 'species') |> 
  ungroup() |>
  mutate(percentcover = pc, 
         `corallite diameter` = `mean_corallite diameter`,          
         `growth rate` = `mean_growth rate`,
         `skeletal density` = `mean_skeletal density`,
         `symbiodinium density` = `mean_symbiodinium density`,
         `colony maximum diameter` = `mean_colony maximum diameter`) |> 
  dplyr::select(site,YEAR, percentcover:`colony maximum diameter`) |> 
  group_by(site,YEAR) |>  ## trait values previously scaled in trait imputation 
  nest(weight = percentcover, data = `corallite diameter`:`colony maximum diameter`) |> 
  mutate(hv = map2(data,weight, \(data,weight) hypervolume_gaussian(data, 
                                                                    name = paste(`site`, YEAR,sep = '_'),
                                                                    weight = weight$percentcover, ###weighted by percent cover 
                                                                    samples.per.point = 1000,
                                                                    kde.bandwidth = estimate_bandwidth(data), 
                                                                    sd.count = 3, 
                                                                    quantile.requested = 0.95, 
                                                                    quantile.requested.type = "probability", 
                                                                    chunk.size = 1000, 
                                                                    verbose = F)),
         hv_size = map_dbl(hv, \(hv) get_volume(hv)),
         centroid = map(hv, \(hv) get_centroid(hv)))
#########STOPPP###
##saveRDS(df_hv, 'data/NEMC_hvs_noboot.rds')
##BIG FILE DO NOT OPEN COMPUTER WILL CRASH, save on desktop 

# SG PR Project Dataset --------------------------------------------------------------
  
## 2023 Hypervolumes: Site-level ------------------------------------------------------------
# Load community data (Taglab Area from LAI)
taglab_df =read.csv("data/SGPR_TaglabAnnot_2023.csv")

SGPR_2023 = taglab_df |> 
    left_join(df_tr, by = 'species') |>  ##same trait data previously used 
    ungroup() |>
    mutate(CoverArea = `TagLab.Area`,
           `corallite diameter` = `mean_corallite diameter`,          
           `growth rate` = `mean_growth rate`,
           `skeletal density` = `mean_skeletal density`,
           `symbiodinium density` = `mean_symbiodinium density`,
           `colony maximum diameter` = `mean_colony maximum diameter`) |> 
    dplyr::select(Site,CoverArea:`colony maximum diameter`) |> 
    group_by(Site) |>  ## trait values previously scaled in trait imputation 
    nest(weight = CoverArea, data = `corallite diameter`:`colony maximum diameter`) |> 
    mutate(hv = map2(data,weight, \(data,weight) hypervolume_gaussian(data, 
                                                                      name = paste(`Site`,sep = '_'),
                                                                      weight = weight$CoverArea, ###weighted by species total area 
                                                                      samples.per.point = 1000,
                                                                      kde.bandwidth = estimate_bandwidth(data), 
                                                                      sd.count = 3, 
                                                                      quantile.requested = 0.95, 
                                                                      quantile.requested.type = "probability", 
                                                                      chunk.size = 1000, 
                                                                      verbose = F)),
           hv_size = map_dbl(hv, \(hv) get_volume(hv)),
           centroid = map(hv, \(hv) get_centroid(hv)))
  
  
#########STOPPP###
#saveRDS(SGPR_2023, "data/SGPR2023_hvs.rds")
##BIG FILE DO NOT OPEN COMPUTER WILL CRASH, save on desktop 
  