devtools::install_github("alko989/gridConstruct/gridConstruct", ref = "lonlatgrid")
## Default inputs
SPECIES   <- "Gadus morhua"
## QUARTER   <- 4
# KM        <- 10
LONSTEP   <- 0.9
LATSTEP   <- 1.5
MINSIZE   <- 4
MAXSIZE   <- 120
MINYEAR   <- 2000
MAXYEAR   <- 2015
BY        <- 1
YEARCLASS <- 2001
MINAGE    <- 1
MAXAGE    <- 8
DATFILE   <- "IBTS.RData"
OUTFILE   <- paste0("results", YEARCLASS, ".RData")

## For scripting
input <- parse(text=Sys.getenv("SCRIPT_INPUT"))
print(input)
eval(input)

## Load data
library(DATRAS)
d <- local({
    load(DATFILE)
    stopifnot(length(ls()) == 1)
    get(ls())
})
stopifnot( class(d) == "DATRASraw" )
d <- subset(d, StatRec != "38G2")

## Make grid
library(gridConstruct)
grid <- gridConstruct(d,icesSquare=TRUE, type = "lonlatGrid", lonstep = LONSTEP, latstep = LATSTEP)
plot(grid)
map("worldHires",add=TRUE)
d <- subset(d, Species == SPECIES)
## Data subset
## d <- addSpectrum(d,cm.breaks=seq(MINSIZE,MAXSIZE,by=BY))
d <- addSpectrum(d, by=BY)

## TODO: Casper sugested to use yearly ALKs 
d <- addNage(d, ages = MINAGE:MAXAGE)
d$haulid <- d$haul.id
##d <- subset(d, Quarter == QUARTER, Gear != "GRT")
d <- subset(d, Year %in% MINYEAR:MAXYEAR )
d <- subset(d, 25<HaulDur & HaulDur<35 )
d <- as.data.frame(d, response="Nage")
d$N <- round(d$Nage)

## Calculate year class
a <- as.numeric(gsub("\\+","",as.character(d$ageGroup)))
y <- as.numeric(as.character(d$Year)) 
d$yearClass <- y - a

## Subset
d <- subset(d, yearClass == YEARCLASS)

## Discretize time - equidistant
breaks <- seq(floor(min(d$abstime)),
              ceiling(max(d$abstime)),
              by= 1/52 )
d$time <- cut(d$abstime, breaks)

library(mapdata)

## Set up time factor (careful with empty factor levels ! )
##years <- as.numeric(levels(d$Year))
##d$time <- factor(d$Year, levels=min(years):max(years))

## Set up spatial factor and grid
## grid <- gridConstruct(d,nearestObs=100)
d$position <- gridFactor(d,grid)
Q0 <- -attr(grid,"pattern")
diag(Q0) <- 0; diag(Q0) <- -rowSums(Q0)
I <- .symDiagonal(nrow(Q0))

## Set up design matrix
## TODO: Gear by sizeGroup
##A <- sparse.model.matrix( ~ sizeGroup:time + Gear - 1, data=d)

d$fac <- factor(paste(d$Year, d$Quarter))
A <- sparse.model.matrix( ~ fac - 1 , data=d)
## A0 <- sparse.model.matrix( ~ sizeGroup:time - 1, data=d)
## A <- sparse.model.matrix( ~ Gear - 1, data=d)
## A <- A[,-which.max(table(d$Gear)),drop=FALSE]

##B <- cbind2(A,A0); B <- t(B)%*%B
B <- t(A)%*%A
if(min(eigen(B)$val)<1e-8)stop("Singular B")

data <- as.list(d)
data$A <- A
data$I <- I
data$Q0 <- Q0

data <- data[!sapply(data,is.character)]
data <- data[!sapply(data,is.logical)]

library(TMB)
compile("model.cpp")
dyn.load(dynlib("model"))
obj <- MakeADFun(
    data=data,
    parameters=list(
        logdelta= 0 ,
        logkappa= 0 ,
        tphi_time= 0 ,
##        tphi_size = 0 ,
        logsigma= 0 ,
        beta= rep(0, ncol(A)) ,
        eta= array(0, c(nrow(Q0), nlevels(d$time) ) ),
        etanug= array(0, c(nlevels(d$haulid) ) )
        ),
    DLL="model",
    random=c("eta","etanug","beta")
    )

print(obj$par)
runSymbolicAnalysis(obj)
system.time(obj$fn())

system.time(opt <- nlminb(obj$par, obj$fn, obj$gr))

system.time(sdr <- sdreport(obj))
pl <- as.list(sdr, "Estimate")
names(pl$beta) <- colnames(data$A)

## for(i in 100)image(grid, concTransform(pl$eta[,i]))
library(animation)
animation::saveVideo({for(i in 1:416) {image(grid, concTransform(pl$eta[,i])); title(levels(d$time)[i])}})

save(grid, sdr, pl, file=OUTFILE)

