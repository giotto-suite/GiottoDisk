#' @import GiottoUtils
#' @import tilework
#' @import methods
#' @importFrom terra crop ext geomtype head tail window `window<-` rasterize centroids spin rescale flip expanse relate
#' @importFrom GiottoClass affine spatShift shear XY calculateOverlap overlapToMatrix setGiotto createExprObj createSpatLocsObj spatUnit featType prov as.terra createGiottoPoints createGiottoPolygon processData analyzeData filterData reduceData spatIDs featIDs spatRelate defaultViewCoordinator resolveRecipe resolveKeep project_region getCellMetadata getFeatureMetadata getSpatialLocations getSpatialEnrichment spatial_smooth_operator feat_cor_from_moments
#' @importClassesFrom GiottoClass affine2d giotto filterParam reduceParam analyzeParam labelProportionsParam cellStatsParam featStatsParam edgeStatsParam spatialSmoothParam featCorParam viewCoordinator dataTableCoordinator cellMetaObj featMetaObj spatLocsObj spatEnrObj giottoPolygon giottoPoints giottoBinPoints exprObj dimObj
#' @importClassesFrom Giotto binarizeThreshParam libraryNormParam logNormParam covGroupsParam covLoessParam varParam scranMarkersParam pcaParam randomPcaParam irlbaPcaParam exactPcaParam autoPcaParam enrichParam pageEnrichParam rankEnrichParam hyperEnrichParam smotifParam dwlsDeconvParam
#' @importFrom Giotto histThreshold
#' @importFrom GiottoUtils %null%
#' @importClassesFrom terra SpatVector
#' @importFrom Matrix rowSums colSums rowMeans colMeans crossprod
#' @importClassesFrom Matrix Matrix
#' @import data.table
NULL
