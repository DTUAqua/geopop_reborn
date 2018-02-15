## apt-get install php-cli
library(DATRAS)

## Overview:
downloadExchange()

## All years:
## downloadExchange("NS-IBTS")
downloadExchange("NS-IBTS", 2000:2003)

d <- readExchangeDir()

save(d, file="all.RData")
