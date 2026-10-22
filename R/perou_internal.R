# Internal helper functions ported VERBATIM from the Perou Lab classifier
# (vendor/perou-lab/PAM50_classifier/PAM50_R/subtypePrediction_functions.R, commit 45121d8). Kept as-is to
# preserve call-identity; only the plotting helpers (myplot/myHeatmap/pcaEA/cols,
# which pull in ctc/heatmap.plus) are split into perou_plots.R. Do not "clean up"
# these — any change risks altering subtype calls.
#
# @keywords internal
# @noRd

collapseIDs <- function(x, allids = row.names(x), method = "mean") {
  allids <- as.vector(allids)
  ids <- levels(as.factor(allids))
  x.col <- NULL

  if (length(ids) == dim(x)[1]) {
    dimnames(x)[[1]] <- allids
    return(x)
  }

  for (i in 1:length(ids)) {
    if (sum(allids == ids[i]) > 1) {
      indices <- allids == ids[i]
      if (method == "mean") {
        vals <- apply(x[indices, ], 2, mean, na.rm = TRUE)
      }
      if (method == "median") {
        vals <- apply(x[indices, ], 2, median, na.rm = TRUE)
      }
      if (method == "stdev") {
        temp <- x[indices, ]
        stdevs <- apply(temp, 1, sd, na.rm = TRUE)
        vals <- temp[match(max(stdevs), stdevs), ]
      }
      if (method == "iqr") {
        temp <- x[indices, ]
        iqrs <- apply(temp, 1, function(x) {
          quantile(x, .75, na.rm = TRUE) - quantile(x, .25, na.rm = TRUE)
        })
        vals <- temp[match(max(iqrs), iqrs), ]
      }
      x.col <- rbind(x.col, vals)
    } else {
      x.col <- rbind(x.col, x[allids == ids[i], ])
    }
  }

  dimnames(x.col) <- list(ids, dimnames(x)[[2]])
  return(x.col)
}

overlapSets <- function(x, y) {
  x <- x[dimnames(x)[[1]] %in% dimnames(y)[[1]], ]
  y <- y[dimnames(y)[[1]] %in% dimnames(x)[[1]], ]
  x <- x[sort.list(row.names(x)), ]
  y <- y[sort.list(row.names(y)), ]
  return(list(x = x, y = y))
}

standardize <- function(x) {
  annAll <- dimnames(x)
  x <- scale(x)
  dimnames(x) <- annAll
  return(x)
}

medianCtr <- function(x) {
  annAll <- dimnames(x)
  medians <- apply(x, 1, median, na.rm = TRUE)
  x <- t(scale(t(x), center = medians, scale = FALSE))
  dimnames(x) <- annAll
  return(x)
}

# Between/within sum-of-squares for a feature across classes. Only used by the
# training/feature-selection path of sspPredict (centroids = FALSE), which the
# PAM50 published-centroid workflow never exercises; defined for completeness.
bwss <- function(x, cl) {
  gm <- mean(x, na.rm = TRUE)
  wss <- 0
  bss <- 0
  for (k in unique(cl)) {
    xi <- x[cl == k]
    m <- mean(xi, na.rm = TRUE)
    wss <- wss + sum((xi - m)^2, na.rm = TRUE)
    bss <- bss + length(xi) * (m - gm)^2
  }
  list(bss = bss, wss = wss)
}

readarray <- function(dataFile, designFile = NA, hr = 1, impute = TRUE, method = "mean") {
  headerRows <- hr
  x <- read.table(dataFile, sep = "\t", header = FALSE, fill = TRUE, stringsAsFactors = FALSE)

  if (headerRows == 1) {
    sampleNames <- as.vector(t(x[1, -1]))
    x <- x[-1, ]
    classes <- NULL
    ids <- x[, 1]
    xd <- x[, -1]
    xd <- apply(xd, 2, as.numeric)
    xd <- collapseIDs(xd, ids, method)
  } else {
    sampleNames <- as.vector(t(x[1, -1]))
    x <- x[-1, ]
    classes <- x[1:(headerRows - 1), ]
    dimnames(classes)[[1]] <- classes[, 1]
    classes <- classes[, -1]
    classes[classes == ""] <- NA
    classes <- t(classes)
    rownames(classes) <- sampleNames
    classes <- as.data.frame(classes)
    xd <- x[(-1:-(headerRows - 1)), ]
    ids <- as.vector(t(xd[, 1]))
    xd <- xd[, -1]
    xd <- apply(xd, 2, as.numeric)
    xd <- collapseIDs(xd, ids, method)
  }

  features <- dim(xd)[1]
  samples <- dim(xd)[2]
  geneNames <- rownames(xd)
  xd <- apply(xd, 2, as.numeric)
  rownames(xd) <- geneNames
  colnames(xd) <- sampleNames

  if (!is.na(designFile)) {
    x <- read.table(designFile, sep = "\t", header = TRUE, row.names = 1, fill = TRUE, stringsAsFactors = FALSE)
    xd <- xd[, sort.list(colnames(xd))]
    xd <- xd[, colnames(xd) %in% rownames(x)]
    x <- x[rownames(x) %in% colnames(xd), ]
    x <- x[sort.list(rownames(x)), ]
    classes <- as.data.frame(x)
  }

  if (sum(apply(xd, 2, is.na)) > 0 & impute) {
    allAnn <- dimnames(xd)
    data.imputed <- impute::impute.knn(as.matrix(xd))$data
    xd <- data.imputed[1:features, ]
    dimnames(xd) <- allAnn
  }

  return(list(xd = xd, classes = classes, nfeatures = features, nsamples = samples, fnames = geneNames, snames = sampleNames))
}

sspPredict <- function(x, classes = "", y, nGenes = "", priors = "equal", std = FALSE, distm = "euclidean", centroids = FALSE) {
  dataMatrix <- x
  features <- dim(x)[1]
  samples <- dim(x)[2]
  sampleNames <- dimnames(x)[[2]]
  featureNames <- dimnames(x)[[1]]

  tdataMatrix <- y
  tfeatures <- dim(y)[1]
  tsamples <- dim(y)[2]
  tsampleNames <- dimnames(y)[[2]]
  tfeatureNames <- dimnames(y)[[1]]

  temp <- overlapSets(dataMatrix, tdataMatrix)
  dataMatrix <- temp$x
  tdataMatrix <- temp$y
  sfeatureNames <- row.names(dataMatrix)

  if (std) {
    dataMatrix <- standardize(dataMatrix)
    tdataMatrix <- standardize(tdataMatrix)
  }

  if (!centroids) {
    thisClass <- as.vector(classes[, 1])
    nClasses <- nlevels(as.factor(thisClass))
    classLevels <- levels(as.factor(thisClass))
    for (j in 1:nClasses) {
      thisClass[thisClass == classLevels[j]] <- j
    }
    thisClass <- as.numeric(thisClass)
    dataMatrix <- dataMatrix[, !(is.na(thisClass))]
    thisClass <- thisClass[!(is.na(thisClass))]

    scores <- apply(dataMatrix, 1, bwss, thisClass)
    trainscores <- vector()
    for (j in 1:dim(dataMatrix)[1]) {
      trainscores[j] <- scores[[row.names(dataMatrix)[j]]]$bss / scores[[row.names(dataMatrix)[j]]]$wss
    }

    dataMatrix <- dataMatrix[sort.list(trainscores, decreasing = TRUE), ]
    tdataMatrix <- tdataMatrix[sort.list(trainscores, decreasing = TRUE), ]

    if (nGenes == "") {
      nGenes <- dim(dataMatrix)[1]
    }
    print(paste("Number of genes used:", nGenes))

    dataMatrix <- dataMatrix[1:nGenes, ]
    tdataMatrix <- tdataMatrix[1:nGenes, ]

    centroids <- matrix(, nrow = nGenes, ncol = nClasses)
    for (j in 1:nClasses) {
      centroids[, j] <- apply(dataMatrix[, thisClass == j], 1, mean)
    }
    dimnames(centroids) <- list(row.names(dataMatrix), NULL)
  } else {
    nGenes <- dim(dataMatrix)[1]
    print(paste("Number of genes used:", nGenes))
    centroids <- dataMatrix
    nClasses <- dim(centroids)[2]
    classLevels <- dimnames(centroids)[[2]]
  }

  distances <- matrix(ncol = nClasses, nrow = dim(tdataMatrix)[2])
  for (j in 1:nClasses) {
    if (distm == "euclidean") {
      distances[, j] <- dist(t(cbind(centroids[, j], tdataMatrix)))[1:(dim(tdataMatrix)[2])]
    }
    if (distm == "correlation" | distm == "pearson") {
      distances[, j] <- t(-1 * cor(cbind(centroids[, j], tdataMatrix), use = "pairwise.complete.obs"))[2:(dim(tdataMatrix)[2] + 1)]
    }
    if (distm == "spearman") {
      distances[, j] <- t(-1 * cor(cbind(centroids[, j], tdataMatrix), method = "spearman", use = "pairwise.complete.obs"))[2:(dim(tdataMatrix)[2] + 1)]
    }
  }

  scores <- apply(distances, 1, min)
  prediction <- vector(length = tsamples)
  for (i in 1:tsamples) {
    prediction[i] <- classLevels[match(scores[i], distances[i, ])]
  }
  names(prediction) <- tsampleNames

  return(list(predictions = prediction, testData = tdataMatrix, distances = distances, centroids = centroids))
}
