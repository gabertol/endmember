library(endmember)
library(tidyverse)

abc<-read.csv('D:/Documentos Drive/Academicos/Publicações/20XX_EM_GZ_botu/Code/data.csv',sep=';')

abcd<-abc %>% select(Coarse.Sand,Medium.Sand,Fine.Sand,Very.Fine.Sand)

abcde<-RECA(abcd)

