#' @import GiottoUtils
#' @import tilework
#' @import methods
#' @importFrom terra crop ext geomtype head tail window `window<-` rasterize centroids spin rescale flip expanse relate
#' @importFrom GiottoClass affine spatShift shear XY calculateOverlap overlapToMatrix setGiotto createExprObj createSpatLocsObj spatUnit featType prov as.terra createGiottoPoints createGiottoPolygon processData analyzeData filterData reduceData spatIDs featIDs spatRelate defaultViewCoordinator resolveRecipe resolveKeep project_region getCellMetadata getFeatureMetadata getSpatialLocations getSpatialEnrichment
#' @importClassesFrom GiottoClass affine2d giotto filterParam reduceParam analyzeParam labelProportionsParam cellStatsParam featStatsParam viewCoordinator dataTableCoordinator cellMetaObj featMetaObj spatLocsObj spatEnrObj giottoPolygon giottoPoints exprObj dimObj
#' @importClassesFrom Giotto binarizeThreshParam libraryNormParam logNormParam covGroupsParam covLoessParam varParam scranMarkersParam pcaParam randomPcaParam irlbaPcaParam exactPcaParam autoPcaParam enrichParam pageEnrichParam rankEnrichParam hyperEnrichParam motifParam autoMotifParam smotifParam
#' @importFrom GiottoUtils %null%
#' @importClassesFrom terra SpatVector
#' @importFrom Matrix rowSums colSums rowMeans colMeans
#' @importClassesFrom Matrix Matrix
#' @import data.table
NULL
