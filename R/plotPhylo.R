#' Plot ggtree-based phylogenetic tree
#'
#' @author Domingos Cardoso
#'
#' @description Function to plot a [ggtree](https://bioconductor.org/packages/release/bioc/html/ggtree.html)-based
#' phylogenetic tree using with the the addition of node supports from different
#' phylogenetic estimation methods.
#'
#' @usage
#' plotPhylo(tree = NULL,
#'           add_raxml_tree = NULL,
#'           add_parsi_tree = NULL,
#'           highlight_clade = NULL,
#'           fill_gradient = NULL,
#'           highlight_taxa = NULL,
#'           highlight_color = NULL,
#'           understate_taxa = NULL,
#'           replace_taxa = NULL,
#'           prune_taxa = NULL,
#'           size_tiplab = NULL,
#'           xlimtree = NULL,
#'           gene_labels = NULL,
#'           ylimgene = NULL,
#'           phylogram_side = FALSE,
#'           phylogram_supports = FALSE,
#'           phylogram_height = NULL,
#'           save = FALSE,
#'           dpi = 600,
#'           dir = "RESULTS_edited_tree",
#'           filename = NULL,
#'           format = "pdf")
#'
#' @param tree A phylo object.
#'
#' @param add_raxml_tree A RaxML-based phylo object.
#'
#' @param add_parsi_tree A parsimony-based phylo object.
#'
#' @param highlight_clade A vector of two tip labels to highlight their clade or
#' a list of multiple vectors with two tip labels each to highlight their clades.
#'
#' @param fill_gradient Any color to highlight an specific clade in the tree.
#'
#' @param highlight_taxa Define a vector with specific tip labels or the name of
#' any taxa (species or genus name) present in the tree to highlight them with
#' an specific color; defaults to "black".
#'
#' @param highlight_color Color to highlight specific tip labels in the tree.
#'
#' @param understate_taxa Define a vector with specific tip labels or the name of
#' any taxa (species or genus name) present in the tree to understate them with
#' a gray color. All other tips will be in black, unless you have set the
#' the arguments \code{highlight_taxa} \code{highlight_color} that define taxa
#' to be highlight with an specific color.
#'
#' @param replace_taxa A vector of tip labels to be substituted with specific
#' alternative names. This is provided in the form of, for instance,
#' c("Harpalyce_formosa" = "Harpalyce_riparia",
#' "Harpalyce_cf_brasiliana" = "Harpalyce_magnibracteata"),
#' where the initial name represents the current tip label to be exchanged with
#' the corresponding alternative name.
#'
#' @param prune_taxa A vector with specific tip labels to be dropped from the tree.
#'
#' @param size_tiplab the size of tip labels; defaults to 2.
#'
#' @param xlimtree Set x axis limits specially for tree panel.
#'
#' @param gene_labels A title for the tree; usually a name of an specific gene
#' (or set of genes) from which the tree was reconstructed.
#'
#' @param ylimgene Set y axis limits for gene labels.
#'
#' @param phylogram_side if \code{TRUE}, a small phylogram with branch lengths
#' will be plotted on the upper left.
#'
#' @param phylogram_supports if \code{TRUE}, the side phylogram will include the
#' same symbols of node support as those plotted in the main edited tree.
#'
#' @param phylogram_height A number to adjust the size of the side phylogram.
#'
#' @param save Logical, if \code{TRUE}, the edited tree will be saved on disk.
#'
#' @param dpi One number in the range of 72-4000 referring to the image
#' resolution in the format of dots per inch in the output file. Default is to
#' create an output with 600 dpi.
#'
#' @param dir Pathway to the computer's directory, where the edited tree file will
#' be saved provided that the argument \code{save} is set up in \code{TRUE}. The
#' default is to create a directory named **RESULTS_edited_tree** and the tree
#' will be saved within a subfolder named after the current date.
#'
#' @param filename Name of the output file to be saved. The default is to
#' create a file entitled **tree**.
#'
#' @param format A character vector related to the file format of the global
#' map to be saved. The default is "jpg" to save the output in Joint
#' Photographic Experts Group (.jpg), but you can also choose "pdf" to save in
#' Portable Document Format (.pdf), "tiff" to save in Tag Image File Format
#' (.tiff) or "png" to save in Portable Network Graphics (.png).
#'
#' @return One or a list of objects of class "ggtree" "gg" "ggplot".
#'
#' @examples
#' \dontrun{
#' library(catGenes)
#'
#' data(Harpalyce_bayes_tree)
#' data(Harpalyce_parsimony_tree)
#' data(Harpalyce_raxml_tree)
#'
#' Harpalyce_clade <- c("Harpalyce_brasiliana_Cardoso2510",
#'                      "Harpalyce_formosa_Hughes2109")
#' outgroup_taxa <- c("Staminodianthus_duckei",
#'                    "Guianodendron_praeclarum",
#'                    "Diplotropis_martiusii",
#'                    "Leptolobium_dasycarpum",
#'                    "Leptolobium_panamense",
#'                    "Bowdichia_virgilioides",
#'                    "Clathrotropis_nitida",
#'                    "Dermatophyllum_secundiflorum")
#'
#' plotPhylo(tree = Harpalyce_bayes_tree,
#'           add_raxml_tree = Harpalyce_raxml_tree,
#'           add_parsi_tree = Harpalyce_parsimony_tree,
#'           highlight_clade = Harpalyce_clade,
#'           fill_gradient = "#D53E4F",
#'           understate_taxa = outgroup_taxa,
#'           size_tiplab = 4,
#'           gene_labels = c("ITS/5.8S", "ETS", "matK", "trnL intron"),
#'           phylogram_side = TRUE,
#'           phylogram_supports = TRUE,
#'           phylogram_height = 25,
#'           save = TRUE,
#'           dir = "results_edited_phylogeny",
#'           filename = "Harpalyce_edited_tree",
#'           format = "pdf")
#'}
#'
#' @importFrom ggtree ggtree MRCA geom_tiplab geom_hilight xlim_tree geom_text2 geom_point2 geom_treescale geom_tree
#' @importFrom treeio isTip
#' @importFrom ggplot2 annotate annotation_custom aes scale_colour_manual ggplotGrob element_rect theme
#' @importFrom cowplot save_plot
#' @importFrom phangorn Descendants
#' @importFrom phytools getDescendants
#' @importFrom ape Ntip
#' @importFrom tibble tibble
#' @importFrom stringr str_extract
#'
#' @export
#'

plotPhylo <- function(tree = NULL,
                      add_raxml_tree = NULL,
                      add_parsi_tree = NULL,
                      highlight_clade = NULL,
                      fill_gradient = NULL,
                      highlight_taxa = NULL,
                      highlight_color = NULL,
                      understate_taxa = NULL,
                      replace_taxa = NULL,
                      prune_taxa = NULL,
                      size_tiplab = NULL,
                      xlimtree = NULL,
                      gene_labels = NULL,
                      ylimgene = NULL,
                      phylogram_side = FALSE,
                      phylogram_supports = FALSE,
                      phylogram_height = NULL,
                      save = FALSE,
                      dpi = 600,
                      dir = "RESULTS_edited_tree",
                      filename = NULL,
                      format = "pdf") {

  #-----------------------------------------------------------------------------
  # Create a new directory to save the plot
  # If there is no directory... make one!

  if (save) {
    if (is.null(filename)) {
      filename <- paste0("tree_",
                         gsub("[/]|[.]", "", paste0(gene_labels, collapse = "_")))
    }

    foldername <- paste0(dir, "/", format(Sys.time(), "%d%b%Y"))
    fullname <- paste0(foldername, "/", filename, ".", format)
    if (!dir.exists(dir)) {
      dir.create(dir)
    }
    if (!dir.exists(foldername)) {
      dir.create(foldername)
    }
  }

  #-----------------------------------------------------------------------------
  if (!any(names(tree@data) %in% "prob")) {
    names(tree@data)[1] = "prob"
  }

  intree <- tree

  if (!is.null(prune_taxa)) {
    tree <- treeio::drop.tip(tree, prune_taxa)
  }

  #-----------------------------------------------------------------------------
  # Delete all taxa from an specific clade. Here example with Mimosoids
  # mims_node <- ggtree::MRCA(prunned_tree, c("Mimosa", "Adenanthera"))
  # del <- prunned_tree$tip.label[getDescendants(prunned_tree, mims_node)]
  # del <- del[!is.na(del)]
  # prunned_mims_tree <- ape::drop.tip(prunned_tree, del[!is.na(del)])
  # # save the prunned tree
  # write.tree(prunned_mims_tree, file = "Data/legs_mimosoids_prunned.tre")

  #-----------------------------------------------------------------------------

  if (!is.null(add_raxml_tree)) {
    raxml_tree <- add_raxml_tree
  } else {
    raxml_tree <- NULL
  }
  if (!is.null(add_parsi_tree)) {
    parsi_tree <- add_parsi_tree
  } else {
    parsi_tree <- NULL
  }

  if (!is.null(add_raxml_tree) |
      !is.null(add_parsi_tree)) {
    # Add support values from different phylogenetic estimation methods
    tree <- .phyloSupports(tree, raxml_tree, parsi_tree)
  }

  # Update specific tip names
  if (!is.null(replace_taxa)) {
    for (i in seq_along(replace_taxa)) {
      tree@phylo$tip.label <- gsub(names(replace_taxa)[i],
                                   replace_taxa[i], tree@phylo$tip.label)
    }
  }

  # Color highlight specific tips
  if (!is.null(highlight_taxa)) {
    nodes <- grep(paste0(highlight_taxa, collapse = "|"), tree@phylo$tip.label)
    tree@data$highlight <- NA
    tree@data$highlight[tree@data$node %in% nodes] <- "highlight_color"
    tree@data$highlight[is.na(tree@data$highlight)] <- "default_color"
  }

  # Understate the color of specific tips
  if (!is.null(understate_taxa)) {
    nodes <- grep(paste0(understate_taxa, collapse = "|"), tree@phylo$tip.label)
    tree@data$understate <- NA
    tree@data$understate[tree@data$node %in% nodes] <- "understated_color"
    tree@data$understate[is.na(tree@data$understate)] <- "default_color"
  }

  if (!is.null(highlight_clade)) {
    if (inherits(highlight_clade, "character")) {
      highlight_clade <- list(highlight_clade)
    }
    for (i in seq_along(highlight_clade)) {
      highlight_clade[[i]] <- gsub("_", " ", highlight_clade[[i]])
    }
  }

  if (!is.null(gene_labels)) {
    gene_labels <- paste0(gene_labels, "\n", collapse = "")
  } else {
    gene_labels <- ""
  }

  if (!any(names(tree@data) %in% "prob")) {
    names(tree@data)[1] = "prob"
  }

  #-----------------------------------------------------------------------------
  # Plotting the figure

  tree@phylo$tip.label <- gsub("_", " ", tree@phylo$tip.label)

  if (is.null(size_tiplab)) size_tiplab = 2

  tree_plot <- ggtree(tree, branch.length = "none", layout = "rectangular",
                      ladderize = TRUE, size = 0.5)

  if (is.null(xlimtree)) {
    xlimtree <- max(tree_plot$data$x)+5
  }

  tree_plot <- tree_plot +
    xlim_tree(xlimtree)

  if (!is.null(highlight_clade)) {
    for (i in seq_along(highlight_clade)) {
      tree_plot <- tree_plot +
        geom_hilight(node = ggtree::MRCA(tree, highlight_clade[[i]]),
                     fill = fill_gradient[i], gradient = TRUE,
                     gradient.direction = "tr", alpha = 0.8) +
        geom_tree(layout="rectangular")
    }
  }

  if (!is.null(highlight_taxa) & is.null(understate_taxa)) {
    tree_plot <- tree_plot +
      geom_tiplab(offset = 0.02, size = size_tiplab, fontface = "italic", aes(color=highlight),
                  show.legend = FALSE) +
      scale_colour_manual(values = c("black", highlight_color))
  }

  if (!is.null(understate_taxa) & is.null(highlight_taxa)) {
    tree_plot <- tree_plot +
      geom_tiplab(offset = 0.02, size = size_tiplab, fontface = "italic", aes(color=understate),
                  show.legend = FALSE) +
      scale_colour_manual(values = c("black", "gray60"))
  }

  if (!is.null(understate_taxa) & !is.null(highlight_taxa)) {
    tf <- tree_plot$data$understate %in% "understated_color"
    tree_plot$data$highlight[tf] <- "understated_color"
    tree_plot <- tree_plot +
      geom_tiplab(offset = 0.02, size = size_tiplab, fontface = "italic", aes(color=highlight),
                  show.legend = FALSE) +
      scale_colour_manual(values = c("black", highlight_color, "gray60"))
  }

  if (is.null(understate_taxa) & is.null(highlight_taxa)) {
    tree_plot <- tree_plot +
      geom_tiplab(offset = 0.02, size = size_tiplab, fontface = "italic",
                  show.legend = FALSE)
  }

  # Add support values below branches
  if (!is.null(add_raxml_tree) |
      !is.null(add_parsi_tree)) {
    if (!any(intree@data$prob > 1)) {
      tree_plot <- tree_plot +
        geom_text2(aes(subset = !isTip & gsub("[/].*", "", prob) >= 0.7,
                       label = prob),
                   size = 3, color = "gray40", hjust = 1.1, vjust = 1.5)
    } else {
      tree_plot <- tree_plot +
        geom_text2(aes(subset = !isTip & gsub("[/].*", "", prob) >= 70,
                       label = prob),
                   size = 3, color = "gray40", hjust = 1.1, vjust = 1.5)
    }

  } else {
    if (!any(intree@data$prob > 1)) {
      tree_plot <- tree_plot +
        geom_text2(aes(subset = !isTip & prob>=0.7 & prob<1,
                       label = stringr::str_extract(prob, "[[:digit:]][.][[:digit:]][[:digit:]]")),
                   size = 3, color = "gray40", hjust = 1.5, vjust = 1.5)
    } else {
      tree_plot <- tree_plot +
        geom_text2(aes(subset = !isTip & prob>=70 & prob<100,
                       label = prob),
                   size = 3, color = "gray40", hjust = 1.5, vjust = 1.5)
    }

  }

  if (is.null(ylimgene)) {
    ylimgene <- (max(tree_plot$data$y)-(max(tree_plot$data$y)*25/100))-10
  }

  if (!any(intree@data$prob > 1)) {
    tree_plot <- tree_plot +
      # Add point when support for node is greater than 0.95 posterior probability)
      geom_point2(aes(subset = prob>=1 & isTip == FALSE), size = 3, shape = 22,
                  fill = "black", alpha = 0.8, stroke = 0.05) +
      geom_point2(aes(subset = prob>=0.9 & prob<1), size = 3, shape = 22,
                  fill = "#0072b2", alpha = 0.8, stroke = 0.05) +
      geom_point2(aes(subset = prob>=0.8 & prob<0.9), size = 3, shape = 22,
                  fill = "#e69f00", alpha = 0.8, stroke = 0.05) +
      geom_point2(aes(subset = prob>=0.7 & prob<0.8), size = 3, shape = 22,
                  fill = "#009e73", alpha = 0.8, stroke = 0.05)
  } else {
    tree_plot <- tree_plot +
      # Add point when support for node is greater than 0.95 posterior probability)
      geom_point2(aes(subset = prob>=100 & isTip == FALSE), size = 3, shape = 22,
                  fill = "black", alpha = 0.8, stroke = 0.05) +
      geom_point2(aes(subset = prob>=90 & prob<100), size = 3, shape = 22,
                  fill = "#0072b2", alpha = 0.8, stroke = 0.05) +
      geom_point2(aes(subset = prob>=80 & prob<90), size = 3, shape = 22,
                  fill = "#e69f00", alpha = 0.8, stroke = 0.05) +
      geom_point2(aes(subset = prob>=70 & prob<80), size = 3, shape = 22,
                  fill = "#009e73", alpha = 0.8, stroke = 0.05)
  }

  tree_plot <- tree_plot +
    # Add gene labels
    annotate(geom = "text", x = 1, y = ylimgene, label = gene_labels,
             color = "gray60", size = 15, fontface = "italic")

  # Create phylogram
  if (phylogram_side) {
    phylogram <- ggtree(tree, layout = "rectangular", ladderize=TRUE, size = 0.1) +
      theme(
        panel.background = element_rect(fill = 'transparent', color = NA),
        plot.background = element_rect(fill = 'transparent', color = NA),
      )

    phylogram <- phylogram +
      geom_treescale(x = 0, y = max(phylogram$data$y)/1.1, color = "gray60",
                     fontsize = 3, linesize = 0.5, offset = 1)

    if (!is.null(highlight_clade)) {
      for (i in seq_along(highlight_clade)) {
        phylogram <- phylogram +
          geom_hilight(node = ggtree::MRCA(tree, highlight_clade[[i]]),
                       fill = fill_gradient[i], gradient = TRUE,
                       gradient.direction = "tr", alpha = 0.8,
                       size = 0.1) +
          geom_tree(layout="rectangular", size = 0.1)
      }
    }

    if (phylogram_supports) {
      if (!any(intree@data$prob > 1)) {
        phylogram <- phylogram +
          geom_point2(aes(subset = prob>=1 & isTip == FALSE), size = 1.5, shape = 22,
                      fill = "black", alpha=0.8, stroke = 0.05) +
          geom_point2(aes(subset = prob>=0.9 & prob<1), size = 1.5, shape = 22,
                      fill = "#0072b2", alpha=0.8, stroke = 0.05) +
          geom_point2(aes(subset = prob>=0.8 & prob<0.9), size = 1.5, shape = 22,
                      fill = "#e69f00", alpha=0.8, stroke = 0.05) +
          geom_point2(aes(subset = prob>=0.7 & prob<0.8), size = 1.5, shape = 22,
                      fill = "#009e73", alpha=0.8, stroke = 0.05)
      } else {
        phylogram <- phylogram +
          geom_point2(aes(subset = prob>=100 & isTip == FALSE), size = 1.5, shape = 22,
                      fill = "black", alpha = 0.8, stroke = 0.05) +
          geom_point2(aes(subset = prob>=90 & prob<100), size = 1.5, shape = 22,
                      fill = "#0072b2", alpha = 0.8, stroke = 0.05) +
          geom_point2(aes(subset = prob>=80 & prob<90), size = 1.5, shape = 22,
                      fill = "#e69f00", alpha = 0.8, stroke = 0.05) +
          geom_point2(aes(subset = prob>=70 & prob<80), size = 1.5, shape = 22,
                      fill = "#009e73", alpha = 0.8, stroke = 0.05)
      }
    }

    if (is.null(phylogram_height)) phylogram_height = 25

    fulltree <- tree_plot + ggplot2::annotation_custom(
      ggplot2::ggplotGrob(phylogram),
      xmin = -1,
      xmax = max(tree_plot$data$x) - (max(tree_plot$data$x)*70/100),
      ymin = max(tree_plot$data$y) - (max(tree_plot$data$y)*phylogram_height/100),
      ymax = max(tree_plot$data$y) - (max(tree_plot$data$y)*1/100))

  } else {
    fulltree <- tree_plot
  }

  if (save) {
    .save_tree(fulltree = fulltree,
               fullname = fullname,
               dpi = dpi,
               format = format,
               base_height = max(tree_plot$data$y)*15/100,
               base_width = max(tree_plot$data$x))
  }


  return(fulltree)
}


#-------------------------------------------------------------------------------
# Auxiliary function to add support values from different phylogenetic estimation
# methods into the same node of a reference tree

.phyloSupports <- function (tree = NULL,
                            raxml_tree = NULL,
                            parsi_tree = NULL) {

  #-----------------------------------------------------------------------------
  # Extract Bayesian posterior values for all internal nodes
  if (is.null(tree)) {
    stop("You must provide a Bayesian tree generated by Mr. Bayes.
          Find help also at DBOSLab-UFBA (Domingos Cardoso; cardosobot@gmail.com)")
  }

  data_tree <- tibble::tibble(tips = NA,
                              isTip = NA,
                              nodes = as.numeric(tree@data$node),
                              prob = tree@data$prob,
                              descendants = NA)

  tf <- data_tree$prob < 1
  data_tree$prob[tf] <- stringr::str_extract(data_tree$prob[tf],
                                             "[[:digit:]][.][[:digit:]][[:digit:]]")

  tf <- data_tree$nodes %in% 1:ape::Ntip(tree@phylo)
  data_tree$tips[tf] <- tree@phylo$tip.label
  data_tree$isTip[tf] <- TRUE
  data_tree$isTip[!tf] <- FALSE
  data_tree$prob[tf] <- NA

  tf <- data_tree$isTip == FALSE
  innodes <- data_tree$nodes[tf]
  for (i in seq_along(innodes)) {
    t <- tree@phylo$tip.label[phangorn::Descendants(tree@phylo,
                                                    innodes[i],
                                                    type = "tips")[[1]]]

    data_tree$descendants[tf][i] <- list(t)
  }


  #-----------------------------------------------------------------------------
  # Extract RAxML bootstrap values for all internal nodes
  if (!is.null(raxml_tree)) {
    data_raxml <- tibble::tibble(tips = NA,
                                 isTip = NA,
                                 nodes = raxml_tree@phylo[["edge"]][,2],
                                 bs = NA,
                                 descendants = NA)

    tf <- data_raxml$nodes %in% 1:ape::Ntip(raxml_tree@phylo)

    data_raxml$tips[tf] <- raxml_tree@phylo$tip.label
    data_raxml$isTip[tf] <- TRUE
    data_raxml$isTip[!tf] <- FALSE
    data_raxml$bs[tf] <- NA

    tf <- data_raxml$nodes %in% raxml_tree@data$node

    data_raxml$bs[tf] <- raxml_tree@data$bootstrap[-1]
    data_raxml$bs[data_raxml$bs < 50] <- NA

    tf <- data_raxml$isTip == FALSE & !is.na(data_raxml$bs)
    innodes <- data_raxml$nodes[tf]
    for (i in seq_along(innodes)) {
      t <- raxml_tree@phylo$tip.label[phangorn::Descendants(raxml_tree@phylo,
                                                            innodes[i],
                                                            type = "tips")[[1]]]

      data_raxml$descendants[tf][i] <- list(t)
    }
  }

  #-----------------------------------------------------------------------------
  # Extract parsimony bootstrap values for all internal nodes
  if (!is.null(parsi_tree)) {
    data_parsi <- tibble::tibble(tips=NA,
                                 isTip = NA,
                                 nodes = parsi_tree@phylo[["edge"]][,2],
                                 bs = round(parsi_tree@phylo[["edge.length"]]),
                                 descendants = NA)

    tf <- data_parsi$nodes %in% 1:ape::Ntip(parsi_tree@phylo)

    data_parsi$tips[tf] <- parsi_tree@phylo$tip.label
    data_parsi$isTip[tf] <- TRUE
    data_parsi$isTip[!tf] <- FALSE
    data_parsi$bs[tf] <- NA

    tf <- data_parsi$isTip == FALSE
    innodes <- data_parsi$nodes[tf]
    for (i in seq_along(innodes)) {
      t <- parsi_tree@phylo$tip.label[phangorn::Descendants(parsi_tree@phylo,
                                                            innodes[i],
                                                            type = "tips")[[1]]]

      data_parsi$descendants[tf][i] <- list(t)
    }
  }

  #-----------------------------------------------------------------------------
  # Combining the branch support values
  tf <- data_tree$isTip == FALSE
  innodes <- data_tree$nodes[tf]

  if (!is.null(raxml_tree)) {
    # Combine raxml values into the bayes tree
    desc <- data_raxml$descendants[!is.na(data_raxml$descendants)]
    node <- data_raxml$nodes[!is.na(data_raxml$descendants)]
    for (i in seq_along(innodes)) {

      temp <- list()
      for (l in seq_along(desc)) {
        tt <- desc[[l]] %in% data_tree$descendants[tf][i][[1]]

        if (length(which(tt == TRUE)) == length(tt) &
            length(tt) == length(data_tree$descendants[tf][i][[1]]))
          temp[[l]] <- node[l]
      }

      bs_raxml <- data_raxml$bs[data_raxml$nodes %in% unlist(temp)]
      if (length(bs_raxml) == 0) {
        bs_raxml <- "-"
      }

      data_tree$prob[tf][i] <- paste0(data_tree$prob[tf][i], "/", bs_raxml)
    }
  }

  if (!is.null(parsi_tree)) {
    # Combine parsimony values into the bayes tree
    desc <- data_parsi$descendants[!is.na(data_parsi$descendants)]
    node <- data_parsi$nodes[!is.na(data_parsi$descendants)]
    for (i in seq_along(innodes)) {

      temp <- list()
      for (l in seq_along(desc)) {
        tt <-  desc[[l]] %in% data_tree$descendants[tf][i][[1]]

        if (length(which(tt == TRUE)) == length(tt) &
            length(tt) == length(data_tree$descendants[tf][i][[1]]))
          temp[[l]] <- node[l]
      }

      bs_parsi <- data_parsi$bs[data_parsi$nodes %in% unlist(temp)]
      if (length(bs_parsi) == 0) {
        bs_parsi <- "-"
      }

      data_tree$prob[tf][i] <- paste0(data_tree$prob[tf][i], "/", bs_parsi)
    }
  }
  tree@data$prob[tf] <- data_tree$prob[tf]

  return(tree)
}


#-------------------------------------------------------------------------------
# Auxiliary function to save the edited tree
.save_tree <- function(fulltree,
                       fullname,
                       dpi,
                       format,
                       base_height,
                       base_width) {

  if (format == "pdf") {
    cowplot::save_plot(fullname,
                       fulltree,
                       ncol = 1, nrow = 1,
                       base_height = base_height,
                       base_width = base_width,
                       base_aspect_ratio = 1.5,
                       limitsize = FALSE)
  } else {
    cowplot::save_plot(fullname,
                       fulltree,
                       ncol = 1, nrow = 1,
                       base_height = base_height,
                       base_width = base_width,
                       base_aspect_ratio = 1.5,
                       limitsize = FALSE,
                       dpi = dpi,
                       bg = "white")
  }
}
