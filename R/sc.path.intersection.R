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
#' \donttest{
#' # Assuming 'scmp_object' is a pre-processed SCMP object with the relevant slots filled
#' sc.path.intersection(scmp_object)
#' }
#' @importFrom S4Vectors isEmpty
#'
sc.path.intersection <- function(scmpObj) {
    
    library(ggupset)
    
    # Check the data
    assertthat::assert_that(
        is(scmpObj, "scMaSigProClass"),
        msg = "Please supply an object of the class 'scMaSigPro'")
    
    # Check if siggenes results exist for groups
    assertthat::assert_that(!S4Vectors::isEmpty(scmpObj@sig.genes@summary),
                            msg = "'sig.genes@Summary' slot is empty, please run 'sc.get.siggenes'")
    
    # Check if more 1 path exist 
    assertthat::assert_that(ncol(scmpObj@sig.genes@summary %>% as.data.frame()) >=2,
                            msg = "'sig.genes@Summary' slot is empty, please run 'sc.get.siggenes'")
    
    # Create nested vector list for the genes
    gene_list <- lapply(scmpObj@sig.genes@summary, function(path){
        # Drop Empty Sets
        path.genes <- path[!(path == " ")]
        return(path.genes)
    })
    
    # Create a unique list of all genes
    all_genes <- unique(unlist(gene_list))
    
    # Initialize the data frame
    gene_df <- data.frame(gene = all_genes)
    
    # Add columns for each pathway
    for (pathway in names(gene_list)) {
        gene_df[[pathway]] <- gene_df$gene %in% gene_list[[pathway]]
    }
    
    # Binarize varibales
    gene_df[,-1] <- apply(gene_df[, -1], 2, FUN = function(X){
        
        if(as.logical(X)){
            return(1)
        }else{
            return(0)
        }
    })
    
    
    
    View(gene_df)
    stop()
    # Convert the logical columns to factors for ggupset
    gene_df[-1] <- lapply(gene_df[-1], factor)
    
    gene_df <- t(gene_df)
    
    colnames(gene_df) <- gene_df[1,]
    gene_df <- as.matrix(gene_df[-1,])
    
    View(gene_df)
    
    df <- gene_df
    
    # Since your dataframe is transposed (genes as columns), you might need to transpose it back
    df <- as.data.frame(t(df))
    
    # Convert the row names to a column if they are not already
    df$pathway <- rownames(df)
    rownames(df) <- NULL
    
    # Convert the TRUE/FALSE columns to factors
    df[-ncol(df)] <- lapply(df[-ncol(df)], factor)
    
    View(df)
    
    # Now use the ComplexUpset package to create the UpSet plot
    p <- upset(df, intersect=colnames(df)[-ncol(df)]) 
    
    
    print()
    
    
    
    
    print(p)
    stop()
    
    # Now you can create your UpSet plot
    library(ggplot2)
    library(ggupset)
    
    View(gene_df)
    p <- ggplot(gene_df, aes(x = as.list(gene_df[-1]))) +
        geom_bar() +
        scale_x_upset(sets = names(gene_list))
    
    print(p)
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    stop()
    gene_df <- stack(gene.list)
    gene_df$ind <- as.factor(gene_df$ind)
    
    p <- ggplot(gene_df, aes(x = values)) +
        geom_bar() +
        scale_x_upset(n_intersections = 5)
    
    print(p)
    View(gene_df)
    
    stop()
    
    #Upset Frame
    upset.frame <- data.frame(geneName = NA,
                              path = NA,
                              expressed = NA)
    upset.frame
    
    # Fill fra,e
    for (i in names(annotated_list)){
        
        # get the df
        df <- annotated_list[[i]]
        
        # Initate temp
        upset.frame.tmp <- data.frame(geneName = NA,
                                  path = NA,
                                  expressed = NA)
        
        # Traverse the dataframe
        for (j in c(1:nrow(df))){
            
            # Extract Annotation
            vec <- df[j,, drop = FALSE]
            
            # Create entry
            entry <- c(geneName = vec[["gene"]],
                       path = i, expressed = vec[["anno"]])
            
            # Add to Dataframe
            upset.frame.tmp <- rbind(upset.frame.tmp, entry)
            
        }
        upset.frame <- rbind(upset.frame, 
                             upset.frame.tmp[-1,])
        
    }    
    upset.frame <- upset.frame[-1,]
    rownames(upset.frame) <- NULL
    upset.frame$expressed <- as.logical(upset.frame$expressed)
    
    
    
    df <- upset.frame
    
    
    # Convert to logical type if necessary
    df$expressed <- as.logical(df$expressed)
    
    # Make sure geneName and path are factors and include all possible levels
    df$geneName <- factor(df$geneName)
    df$path <- factor(df$path)
    
    # Ensure all combinations of geneName and path are present
    complete_df <- df %>%
        expand(geneName, path) %>%
        left_join(df, by = c("geneName", "path"))
    
    # Replace NA with FALSE in the expressed column for missing combinations
    complete_df$expressed[is.na(complete_df$expressed)] <- FALSE
    
    # Cast the complete data frame to a matrix format
    df_matrix <- dcast(complete_df, path ~ geneName, value.var = "expressed", drop=FALSE)
    
    # View the matrix
    View(df_matrix)
    
    stop()
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    complete_df <- upset.frame %>%
        complete(geneName, path, fill = list(expressed = FALSE))
    
    upset_matrix <- dcast(complete_df, path ~ geneName, value.var = "expressed")
    
    View(upset_matrix)
    
    
    stop()
    
    # Assuming your data frame is named df
    all_combinations <- expand.grid(geneName = unique(upset.frame$geneName),
                                    path = unique(upset.frame$path))
    full_df <- merge(all_combinations, upset.frame, all = TRUE)
    full_df$expressed[is.na(full_df$expressed)] <- FALSE
    
    upset_matrix <- dcast(full_df, path ~ geneName, value.var = "expressed")
    
    
    View(upset_matrix)
    
    stop()
    
    stop()
    
    annotated_list <- imap(annotated_list, ~mutate(.x, path = .y))
    
    annotated_list <- reduce(annotated_list, full_join, by = "gene")
    
    View(annotated_list)
    
    stop()
    
    # Get list of all genes
    gene.vector <- lapply(colnames(summary_df), function(group){
        gene.vector <- unique(summary_df[[group]])
        return(gene.vector)
    })
    
    
    View(as.list(scmpObj@sig.genes@summary))
    
    stop()
    
    # Check if atleast two paths exist
    assertthat::assert_that(nrow(summary_df) >=2,
                            msg = "Please run 'sc.get.siggenes' with 'each' or 'groups'")
    
    
    print(head(summary_df[c(1:5),]))
    
    summary_df  <- summary_df %>% 
        pivot_longer(cols = everything(), names_to = "Path", values_to = "Gene") %>%
        mutate(Value = TRUE)
    
    View(summary_df)
    
    stop()
    # Create Path Memebership
    unique_genes <- sort(unique(c(summary_df$Path1, summary_df$Path2vsPath1)))
    pathway_membership <- data.frame(Gene = unique_genes)
    
    # Determine membership in pathways
    pathway_membership$Path1 <- pathway_membership$Gene %in% summary_df$Path1
    pathway_membership$Path2vsPath1 <- pathway_membership$Gene %in% summary_df$Path2vsPath1
    
    # Convert to long format for plotting
    tidy_pathway_membership <- pathway_membership %>%
        gather(Pathway, Member, -Gene) %>%
        filter(Member)
    
    print(tidy_pathway_membership)
    
    p <- ggplot(tidy_pathway_membership, aes(x = as.factor(Gene), fill = Pathway)) +
        geom_bar() +
        scale_x_discrete(drop = FALSE) +
        scale_fill_manual(values = c("Path1" = "blue", "Path2vsPath1" = "red")) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        labs(x = "Gene", y = "Count", fill = "Pathway") +
        scale_x_upset(n_intersections = 10)
    
    print(p)
    
    stop()
    
    # get gene list per path
    path 
    
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
