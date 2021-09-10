addPHATE <- function (ArchRProj = NULL, reducedDims = "IterativeLSI", method = "phateR",decay=15,
    name = "PHATE", dimsToUse = NULL, scaleDims = NULL, corCutOff = 0.75, saveModel = FALSE,
    verbose = TRUE, seed = 1, force = FALSE, threads = max(floor(getArchRThreads()/2),
        1), ...)
{
    .validInput(input = ArchRProj, name = "ArchRProj", valid = c("ArchRProj"))
    .validInput(input = reducedDims, name = "reducedDims", valid = c("character"))
    .validInput(input = method, name = "method", valid = c("character"))
    .validInput(input = name, name = "name", valid = c("character",
        "null"))
    .validInput(input = dimsToUse, name = "dimsToUse", valid = c("integer",
        "null"))
    .validInput(input = decay, name = "decay", valid = c("integer",
        "null"))
    .validInput(input = scaleDims, name = "scaleDims", valid = c("boolean",
        "null"))
    .validInput(input = corCutOff, name = "corCutOff", valid = c("numeric",
        "null"))
    .validInput(input = verbose, name = "verbose", valid = c("boolean"))
    .validInput(input = seed, name = "seed", valid = c("integer"))
    .validInput(input = force, name = "force", valid = c("boolean"))
    .validInput(input = threads, name = "threads", valid = c("integer"))
    if (name %in% names(ArchRProj@embeddings)) {
        if (!force) {
            stop("Embedding Already Exists! Either set force = TRUE or use a different name!")
        }
    }
    embeddingParams <- list(...)
    embeddingParams$X <- getReducedDims(ArchRProj = ArchRProj,
        reducedDims = reducedDims, dimsToUse = dimsToUse, corCutOff = corCutOff,
        scaleDims = scaleDims)
    if (tolower(method) != "phater") {
        message("Methods other than phateR not currently supported, defaulting to Rtsne!")
        method <- "phateR"
    }
    if (tolower(method) == "phater") {
        .requirePackage("ReductionWrappers", source = "cran")
        embeddingParams$X <- getReducedDims(ArchRProj = ArchRProj,
            reducedDims = reducedDims, dimsToUse = dimsToUse,
            corCutOff = corCutOff, scaleDims = scaleDims)
        phate_pos <- phateR::phate(embeddingParams$X,n.jobs=-1,decay=decay)
        dfEmbedding <- data.frame(phate_pos$embedding)
        colnames(dfEmbedding) <- paste0(reducedDims, "#PHATE_Dimension_",
            seq_len(ncol(dfEmbedding)))
        rownames(dfEmbedding) <- rownames(embeddingParams$X)
    } else {
        stop("PHATE Method Not Currently Supported!")
    }
    embeddingParams$X <- NULL
    ArchRProj@embeddings[[name]] <- SimpleList(df = dfEmbedding,
        params = embeddingParams)
    return(ArchRProj)
}
environment(addPHATE) <- asNamespace('ArchR')

addDiffusionMap <- function (ArchRProj = NULL, reducedDims = "IterativeLSI", method = "destiny",sigma.use=30,
    name = "DiffusionMap", dimsToUse = NULL, scaleDims = NULL, corCutOff = 0.75, saveModel = FALSE,
    verbose = TRUE, seed = 1, force = FALSE, threads = max(floor(getArchRThreads()/2),
        1), ...)
{
    .validInput(input = ArchRProj, name = "ArchRProj", valid = c("ArchRProj"))
    .validInput(input = reducedDims, name = "reducedDims", valid = c("character"))
    .validInput(input = method, name = "method", valid = c("character"))
    .validInput(input = name, name = "name", valid = c("character",
        "null"))
    .validInput(input = dimsToUse, name = "dimsToUse", valid = c("integer",
        "null"))
    .validInput(input = sigma.use, name = "sigma.use", valid = c("integer",
        "null"))
    .validInput(input = scaleDims, name = "scaleDims", valid = c("boolean",
        "null"))
    .validInput(input = corCutOff, name = "corCutOff", valid = c("numeric",
        "null"))
    .validInput(input = verbose, name = "verbose", valid = c("boolean"))
    .validInput(input = seed, name = "seed", valid = c("integer"))
    .validInput(input = force, name = "force", valid = c("boolean"))
    .validInput(input = threads, name = "threads", valid = c("integer"))
    if (name %in% names(ArchRProj@embeddings)) {
        if (!force) {
            stop("Embedding Already Exists! Either set force = TRUE or use a different name!")
        }
    }
    embeddingParams <- list(...)
    embeddingParams$X <- getReducedDims(ArchRProj = ArchRProj,
        reducedDims = reducedDims, dimsToUse = dimsToUse, corCutOff = corCutOff,
        scaleDims = scaleDims)
    if (tolower(method) != "destiny") {
        message("Methods other than destiny not currently supported, defaulting to Rtsne!")
        method <- "destiny"
    }
    if (tolower(method) == "destiny") {
        .requirePackage("ReductionWrappers", source = "cran")
        embeddingParams <- list()
        embeddingParams$X <- getReducedDims(ArchRProj = ArchRProj,
            reducedDims = reducedDims, dimsToUse = dimsToUse,
            corCutOff = corCutOff, scaleDims = scaleDims)
        knn <- destiny::find_dm_k(ncol(embeddingParams$X))
        print(paste("Using provided global sigma", sigma.use))
        print(paste("destiny will use", knn, "nearest neighbors."))
        dm <- destiny::DiffusionMap(embeddingParams$X, sigma = sigma.use, k = knn,
            n_eigs = 200, density_norm = TRUE, distance = "euclidean")
        dfEmbedding <- data.frame(DC1=dm$DC1,DC2=dm$DC2,DC3=dm$DC3,DC4=dm$DC4,DC5=dm$DC5,
            DC6=dm$DC6,DC7=dm$DC7,DC8=dm$DC8,DC9=dm$DC9,DC10=dm$DC10)
        colnames(dfEmbedding) <- paste0(reducedDims, "#DiffusionMap_Dimension_",
            seq_len(ncol(dfEmbedding)))
        rownames(dfEmbedding) <- rownames(embeddingParams$X)
    } else {
        stop("DiffusionMap Method Not Currently Supported!")
    }
    embeddingParams$X <- NULL
    ArchRProj@embeddings[[name]] <- SimpleList(df = dfEmbedding,
        params = embeddingParams)
    return(ArchRProj)
}
environment(addDiffusionMap) <- asNamespace('ArchR')
