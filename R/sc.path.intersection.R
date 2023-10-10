#' Perform UpSet Plot on Intersection of Significant Genes from SCMP Object
#'
#' This function takes a Single Cell MultiPathway (SCMP) object, extracts the list of 
#' significant genes for each pathway, and performs an UpSet plot to visualize their intersections.
#'
#' @param scmpObj An object of class SCMP, which should contain the @siggenes@summary slot
#' filled with the summary of significant genes for each pathway.
#'
#' @return An UpSet plot visualizing the intersections of significant genes across pathways.
#' 
#' @examples
#' \dontrun{
#' # Assuming 'scmp_object' is a pre-processed SCMP object with the relevant slots filled
#' sc.path.intersection(scmp_object)
#' }
#'
sc.path.intersection <- function(scmpObj) {
    
    # Extract the list of significant genes for each pathway
    gene.list <- lapply(colnames(scmpObj@siggenes@summary), function(path) {
        return(na.omit(scmpObj@siggenes@summary[[path]]))
    })
    names(gene.list) <- colnames(scmpObj@siggenes@summary)
    
    # Perform UpSet plot
    upset(
        fromList(gene.list),
        main.bar.color = "#F58A53",
        matrix.color = "#15918A",
        line.size = 1.5,
        point.size = 3,
        shade.color = "purple",
        text.scale = 1.5,
        sets.x.label = "Number of Features",
        sets.bar.color = "#EE446F"
    )
}
