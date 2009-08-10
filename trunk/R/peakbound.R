peakbound=function(x){ 
    ##assuming 500 bp window and flanking 500 bp windows
    ##TO DO: 
    ##1)make general so can handle any window size
    ##2)make to possibly handle multiple peaks
    ##3)make so can look for max in flanks if which.max=border
    ##############################
    #left bound
    #find max of center window
    max=which.max(x[501:1000])+500
    hold=matrix(0, 1,500)
    for(bound in 500:50){
        X=(max-bound):max
        Y=x[(max(max-bound,1)):max]
        B=sum((X-mean(X))*Y)/sum((X-mean(X))^2)
        hold[bound]=(B^2)*sum((X-mean(X))^2)/sum((Y-mean(Y))^2)
    }

    fit=which.max(hold)+1
    peakstart=max-fit
    ##############################
    #right bound
    hold=matrix(0, 1,500)
    for(bound in 50:500){
        X=max:(max+bound)
        Y=x[(max:(max+bound))]
        B=sum((X-mean(X))*Y)/sum((X-mean(X))^2)
        hold[bound]=(B^2)*sum((X-mean(X))^2)/sum((Y-mean(Y))^2)
    }
    fit=which.max(hold)
    peakend=max+fit
    ##############################
    return(c(peakstart, peakend))
}

