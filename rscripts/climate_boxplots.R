setwd("C:/Users/Michel/Documents/eukaryotes/data/")
load("200_all_data_long_export_filtered-_climate.RData")
f <- psob_molten
rm(psob_molten)
colnames(f)
attach(f)
precip <- as.numeric(RACMO_precip_mm_35to1km )
tmp <- as.numeric(RACMO_tmp_2m_35to1km)
windsp <- as.numeric(RACMO_windsp_10m_35to1km)
head(precip)

par(mfrow=c(2,2))
boxplot(precip ~ Location,main="precipitation")
boxplot(tmp ~ Location,main="temperature")
boxplot(windsp ~ Location,main="windspeed")
boxplot(as.numeric(SLOPE) ~ Location, ,main="slope")

