## Note, there must be at least one age sample for each age group for each year
## There is a function in the "surveyIndex" package called "fixAgeGroup" to take care of this
applyALK<-function(d,spatial, ages){
    d.ysplit = split(d,d$Year) ## estimate separate ALK for each year

    ## Declare settings for ALK model
    mf = ""
    ack=TRUE;
    useBICs=TRUE;
    varCofs=FALSE;
    maxKs=50;
    
    if(!spatial) d.ALK= lapply(d.ysplit,fitALK,minAge=min(ages),maxAge=max(ages),method=1)
    if(spatial) d.ALK= lapply(d.ysplit,fitALK,minAge=min(ages),maxAge=max(ages),autoChooseK=ack,useBIC=useBICs,varCof=varCofs,maxK=maxKs)


    ## Predict numbers-at-age from ALK and lengths
    d.Nage=lapply(d.ALK,predict)
    for(i in 1:length(d.ALK)) d.ysplit[[i]]$Nage=d.Nage[[i]];
    dd <- do.call("c",d.ysplit)
    d.ysplit=NULL;
    
    dd
}
