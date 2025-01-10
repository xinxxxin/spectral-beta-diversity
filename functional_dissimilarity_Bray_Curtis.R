# A new parametric measure of functional dissimilarity: Bridging the gap between the Bray-Curtis dissimilarity and the Euclidean distance
# Functional dissimilarity function from Ricotta and Pavoine (2022) (https://www.sciencedirect.com/science/article/pii/S0304380022000084#sec0008)

spedisparam <-
function(comm, method = c("D", "Delta", "L", "M"), abundance = c("relative", "absolute"), alpha = 2, tol = 1e-8, ...)
{
    q <- alpha[1]
    if(!is.numeric(q)) stop("q must be a numeric")
    if(!inherits(comm, "matrix") && !inherits(comm, "data.frame")) stop("Incorrect definition of comm")
    method <- method[1]
    if(!method%in%c("D", "Delta", "L", "M")) stop("Incorrect definition for method")
    abundance <- abundance[1]
    if(!abundance%in%c("relative", "absolute"))
stop("Incorrect definition for abundance")
    dataset <- t(comm)
    total <- colSums(dataset)
    if(abundance == "relative")
        abu <- sweep(dataset, 2, total, "/")
    else abu <- dataset
    num.plot <- dim(dataset)[2]
    num.sp <- dim(dataset)[1]
    names<-list(colnames(dataset), colnames(dataset))
    dis.matrix<-matrix(0, nrow=num.plot, ncol=num.plot, dimnames=names)
    for (i in 2:num.plot) {
    for (j in 1:(i-1)) {
        Zik <- abu[, j]
        Zih <- abu[, i]
        tabZ <- rbind.data.frame(Zik, Zih)
        garde <- apply(tabZ, 2, sum)>tol
	   Zik <- Zik[garde]
        Zih <- Zih[garde]
        tabZ <- tabZ[, garde]

        NUM <- (sum(((abs(Zik-Zih))^q)))^(1/q)
         
        if(method == "D")
            DEN <- (sum(((Zik+Zih)^q)))^(1/q)
        else if(method == "Delta")
            DEN <- ( sum((Zik^q))+ sum((Zih^q)) )^(1/q)
        else if(method == "M")
            DEN <- 1
        else
            DEN <- (sum((Zik^q)))^(1/q) + (sum((Zih^q)))^(1/q)
        index <- NUM/DEN
        dis.matrix[i, j] <- index
    }
    }
    dis.matrix <- dis.matrix + t(dis.matrix)
return(as.dist(dis.matrix, ...))
}

fundisparam <-
function(comm, dis, method = c("D", "Delta", "L", "M"), abundance = c("relative", "absolute"), alpha = 2, tol = 1e-8, ...)
{
    q <- alpha[1]
    if(!is.numeric(q)) stop("q must be a numeric")
    if(inherits(dis, "dist")) dis <- as.matrix(dis)
    if(!inherits(dis, "matrix")) stop("Incorrect definition of dis")
    if(any(dis< (-tol))) stop("Incorrect definition of dis")
    dis[dis<0] <- 0
    if(!inherits(comm, "matrix") && !inherits(comm, "data.frame")) stop("Incorrect definition of comm")
    if(any(dis>1)){
        warning("dissimilarities in dis are not in the range 0-1. They have been normalized by the maximum")
        dis <- dis/max(dis)
    }
    if(any(!colnames(comm) %in%rownames(dis))) stop("At least one species in the matrix of abundances is missing in the matrix of dissimilarities")
    if(any(!colnames(comm) %in%colnames(dis))) stop("At least one species in the matrix of abundances is missing in the matrix of dissimilarities")
    method <- method[1]
    if(!method%in%c("D", "Delta", "L", "M")) stop("Incorrect definition for method")
    abundance <- abundance[1]
    if(!abundance%in%c("relative", "absolute"))
stop("Incorrect definition for abundance")
dis <- dis[colnames(comm), colnames(comm)]
    dataset <- t(comm)
    similarities <- 1-as.matrix(dis)
    total <- colSums(dataset)
    if(abundance == "relative")
        abu <- sweep(dataset, 2, total, "/")
    else abu <- dataset
    num.plot <- dim(dataset)[2]
    num.sp <- dim(dataset)[1]
    names<-list(colnames(dataset), colnames(dataset))
    dis.matrix<-matrix(0, nrow=num.plot, ncol=num.plot, dimnames=names)
    for (i in 2:num.plot) {
    for (j in 1:(i-1)) {
        mat_folk <- similarities*abu[, j]
        mat_folk2 <- similarities*abu[, i]
        Zik <- colSums(mat_folk)
        Zih <- colSums(mat_folk2)
        tabZ <- rbind.data.frame(Zik, Zih)
        garde <- apply(tabZ, 2, sum)>tol & apply(abu[,c(i,j)], 1, sum)>tol
	   Zik <- Zik[garde]
        Zih <- Zih[garde]
        tabZ <- tabZ[, garde]
        wk <- (abu[, j]+abu[, i])/sum(abu[, c(i,j)])
        wk <- wk[garde]

        NUM <- (sum((wk*(abs(Zik-Zih))^q)))^(1/q)
         
        if(method == "D")
            DEN <- (sum((wk*(Zik+Zih)^q)))^(1/q)
        else if(method == "Delta")
            DEN <- ( sum(wk*(Zik^q))+ sum(wk*(Zih^q)) )^(1/q)
        else if(method == "M")
            DEN <- 1
        else
            DEN <- (sum(wk*(Zik^q)))^(1/q) + (sum(wk*(Zih^q)))^(1/q)
        index <- NUM/DEN
        dis.matrix[i, j] <- index
    }
    }
    dis.matrix <- dis.matrix + t(dis.matrix)
return(as.dist(dis.matrix, ...))
}
