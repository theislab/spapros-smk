# Performs feature/gene selection
#
# More detailed description
#
# @param data: Numeric matrix where samples are in row and features are in column.
# @param k: Number of features to be selected.
#
# @return Vector containing ID of the features.
#
# @examples
# GeneID = SelfE(Data, 50)

## Set libPaths for jobs on icb cluster
#if (.libPaths() == "/opt/R/lib/R/library") {
#    .libPaths(c("/home/louis.kuemmerle/bin",.libPaths()))
#}

SelfE <- function(data, k)
{
    idx <- vector()
    Y <- as.matrix(data)
    Phi <- as.matrix(data)
    R <- as.matrix(Y)

    for (iter in 1:k)
    {
        print(iter)
        K <- abs((t(R)) %*% R)
        
        c <- vector()
        for (i in 1:ncol(Phi))
        {
            c[i] <- norm(as.matrix(K[i, ]), '2')
        }
        pos <- order(-c)
        idx <- append(idx, pos[1])

        PhiS <- as.matrix(Phi[, idx])
        Yiter <- PhiS %*% pinv(PhiS) %*% Y
        R <- Y - Yiter
    }
    return(idx)
}
