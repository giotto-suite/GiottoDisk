#' @import GiottoUtils
#' @import tilework
#' @import methods
#' @importFrom terra crop ext geomtype head tail window `window<-` rasterize centroids spin rescale flip expanse relate
#' @importFrom GiottoClass affine spatShift shear XY calculateOverlap overlapToMatrix setGiotto createExprObj createSpatLocsObj spatUnit featType prov as.terra createGiottoPoints createGiottoPolygon processData analyzeData filterData reduceData spatIDs spatRelate defaultViewCoordinator resolveSubobject prepareIds getCellMetadata getFeatureMetadata getSpatialLocations getSpatialEnrichment
#' @importClassesFrom GiottoClass affine2d giotto filterParam reduceParam analyzeParam labelProportionsParam viewCoordinator dataTableCoordinator cellMetaObj featMetaObj spatLocsObj spatEnrObj giottoPolygon giottoPoints exprObj dimObj
#' @importClassesFrom Giotto binarizeThreshParam cellStatsParam featStatsParam libraryNormParam logNormParam covGroupsParam covLoessParam varParam scranMarkersParam pcaParam randomPcaParam irlbaPcaParam exactPcaParam autoPcaParam
#' @importFrom GiottoUtils %null%
#' @importClassesFrom terra SpatVector
#' @importClassesFrom Matrix Matrix
#' @import data.table
NULL
