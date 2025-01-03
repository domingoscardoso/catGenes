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
#'           layout = "rectangular",
#'           branch.width = 0.5,
#'           branch.supports = TRUE,
#'           add.raxml.tree = NULL,
#'           add.parsi.tree = NULL,
#'           highlight.clade = NULL,
#'           fill.gradient = NULL,
#'           show.tip.label = TRUE,
#'           size.tip.label = 2,
#'           fontface.tip.label = "italic",
#'           fancy.tip.label = FALSE,
#'           abbrev.tip.label = FALSE,
#'           highlight.taxa = NULL,
#'           highlight.color = NULL,
#'           understate.taxa = NULL,
#'           replace.taxa = NULL,
#'           prune.taxa = NULL,
#'           xlim.tree = NULL,
#'           hexpand = NULL,
#'           gene.label = NULL,
#'           ylim.gene.label = NULL,
#'           phylogram.side = FALSE,
#'           phylogram.supports = FALSE,
#'           phylogram.height = NULL,
#'           save = FALSE,
#'           dpi = 600,
#'           dir = "RESULTS_edited_tree",
#'           filename = NULL,
#'           format = "pdf",
#'           ...)
#'
#' @param tree A phylo object.
#'
#' @param layout One of 'rectangular', 'dendrogram', 'slanted', 'ellipse',
#' 'roundrect', 'fan', 'circular', 'inward_circular', 'radial', 'equal_angle',
#' 'daylight' or 'ape'
#'
#' @param branch.width A numeric vector giving the width of the branches of
#' the plotted phylogeny. These are taken to be in the same order than the
#' component \code{edge} of \code{phylo}. If fewer widths are given than the
#' length of edge, then these are recycled. Defaults to 0.5.
#'
#' @param branch.supports A logical indicating whether to show the node supports
#' below the branches of the phylogeny. Defaults to TRUE, i.e. the statistical
#' supports are shown.
#'
#' @param add.raxml.tree A RaxML-based phylo object.
#'
#' @param add.parsi.tree A parsimony-based phylo object.
#'
#' @param highlight.clade A vector of two tip labels to highlight their clade or
#' a list of multiple vectors with two tip labels each to highlight their clades.
#'
#' @param fill.gradient One or more colors to highlight specific clades in the
#' tree.
#'
#' @param show.tip.label A logical indicating whether to show the tip labels on
#' the phylogeny. Defaults to TRUE, i.e. the labels are shown.
#'
#' @param size.tip.label The size of tip labels; defaults to 2.
#'
#' @param fontface.tip.label The font face of text, defaults to 3 (italic), others
#' are 1 (plain), 2 (bold), 4 (bold.italic).
#'
#' @param fancy.tip.label A logical indicating whether to display tip labels with
#' distinct formatting for the taxon name and associated voucher or GenBank
#' accession number. Defaults to FALSE, maintaining the labels in the default
#' format specified by \code{fontface.tip.label}.
#'
#' @param abbrev.tip.label A logical indicating whether to show the tip labels
#' on the phylogeny. Defaults to FALSE, i.e. the labels are not abbreviated.
#'
#' @param highlight.taxa Define a vector with specific tip labels or the name of
#' any taxa (species or genus name) present in the tree to highlight them with
#' specific colors; defaults to "black".
#'
#' @param highlight.color Color vector to highlight specific tip labels in the
#' tree.
#'
#' @param understate.taxa Define a vector with specific tip labels or the name of
#' any taxa (species or genus name) present in the tree to understate them with
#' a gray color. All other tips will be in black, unless you have set the
#' the arguments \code{highlight.taxa} \code{highlight.color} that define taxa
#' to be highlight with an specific color.
#'
#' @param replace.taxa A vector of tip labels to be substituted with specific
#' alternative names. This is provided in the form of, for instance,
#' c("Harpalyce_formosa" = "Harpalyce_riparia",
#' "Harpalyce_cf_brasiliana" = "Harpalyce_magnibracteata"),
#' where the initial name represents the current tip label to be exchanged with
#' the corresponding alternative name.
#'
#' @param prune.taxa A vector with specific tip labels to be dropped from the tree.
#'
#' @param xlim.tree Set x axis limits specially for tree panel.
#'
#' @param hexpand Expand x axis limits by ratio of x axis range.
#'
#' @param gene.label A title for the tree; usually a name of an specific gene
#' (or set of genes) from which the tree was reconstructed.
#'
#' @param size.gene.label The size of gene text labels; defaults to 12.
#'
#' @param ylim.gene.label Set y axis limits for gene labels.
#'
#' @param phylogram.side If \code{TRUE}, a small phylogram with branch lengths
#' will be plotted on the upper left.
#'
#' @param phylogram.supports If \code{TRUE}, the side phylogram will include the
#' same symbols of node support as those plotted in the main edited tree.
#'
#' @param phylogram.height A number to adjust the size of the side phylogram.
#'
#' @param save Logical, if \code{TRUE}, the edited tree will be saved on disk.
#'
#' @param dpi One number in the range of 72-4000 referring to the image
#' resolution in the format of dots per inch in the output file. Default is to
#' create an output with 600 dpi.
#'
#' @param dir The path to the directory where the edited tree file will be saved
#' provided that the argument \code{save} is set up in \code{TRUE}. The default
#' is to create a directory named **RESULTS_edited_tree** and the tree will be
#' saved within a subfolder named after the current date.
#'
#' @param filename Name of the output file to be saved. The default is to
#' create a file entitled **edited_tree** and the specified \code{layout}. In the
#' case of using the default layout, the saved file will be named
#' **edited_tree_rectangular**
#'
#' @param format A character vector related to the file format of the global
#' map to be saved. The default is 'jpg' to save the output in Joint
#' Photographic Experts Group (.jpg), but you can also choose 'pdf' to save in
#' Portable Document Format (.pdf), 'tiff' to save in Tag Image File Format
#' (.tiff) or "png" to save in Portable Network Graphics (.png).
#'
#' @param ... Further arguments to be passed to \code{ggtree}.
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
#'           layout = "rectangular",
#'           branch.width = 0.5,
#'           branch.supports = TRUE,
#'           add.raxml.tree = Harpalyce_raxml_tree,
#'           add.parsi.tree = Harpalyce_parsimony_tree,
#'           highlight.clade = Harpalyce_clade,
#'           fill.gradient = "#D53E4F",
#'           show.tip.label = TRUE,
#'           size.tip.label = 4,
#'           fontface.tip.label = "italic",
#'           understate.taxa = outgroup_taxa,
#'           gene.label = c("ITS/5.8S", "ETS", "matK", "trnL intron"),
#'           size.gene.label = 12,
#'           ylim.gene.label = NULL,
#'           phylogram.side = TRUE,
#'           phylogram.supports = TRUE,
#'           phylogram.height = 25,
#'           save = TRUE,
#'           dir = "results_edited_phylogeny",
#'           filename = "Harpalyce_edited_tree",
#'           format = "pdf")
#'}
#'
#' @importFrom ggtree ggtree MRCA geom_tiplab geom_hilight xlim_tree geom_text2 geom_point2 geom_treescale geom_tree hexpand td_filter
#' @importFrom ggtext geom_richtext
#' @importFrom treeio isTip
#' @importFrom ggplot2 annotate annotation_custom aes scale_colour_manual ggplotGrob element_rect theme
#' @importFrom cowplot save_plot
#' @importFrom phangorn Descendants
#' @importFrom phytools getDescendants
#' @importFrom ape Ntip
#' @importFrom tibble tibble
#' @importFrom stringr str_extract
#' @importFrom dplyr mutate
#' @importFrom glue glue
#'
#' @export
#'

plotPhylo <- function(tree = NULL,
                      layout = "rectangular",
                      branch.width = 0.5,
                      branch.supports = TRUE,
                      add.raxml.tree = NULL,
                      add.parsi.tree = NULL,
                      highlight.clade = NULL,
                      fill.gradient = NULL,
                      show.tip.label = TRUE,
                      size.tip.label = 4,
                      fontface.tip.label = "italic",
                      fancy.tip.label = FALSE,
                      abbrev.tip.label = FALSE,
                      highlight.taxa = NULL,
                      highlight.color = NULL,
                      understate.taxa = NULL,
                      replace.taxa = NULL,
                      prune.taxa = NULL,
                      xlim.tree = NULL,
                      hexpand = NULL,
                      gene.label = NULL,
                      size.gene.label = 12,
                      ylim.gene.label = NULL,
                      phylogram.side = FALSE,
                      phylogram.supports = FALSE,
                      phylogram.height = NULL,
                      save = FALSE,
                      dpi = 600,
                      dir = "RESULTS_edited_tree",
                      filename = NULL,
                      format = "pdf",
                      ...) {

  #-----------------------------------------------------------------------------
  # Create a new directory to save the plot
  # If there is no directory... make one!

  if (save) {
    if (is.null(filename)) {
      filename <- paste0("edited_tree_", layout)
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

  if (!is.null(prune.taxa)) {
    tree <- treeio::drop.tip(tree, prune.taxa)
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

  if (!is.null(add.raxml.tree)) {
    raxml_tree <- add.raxml.tree
  } else {
    raxml_tree <- NULL
  }
  if (!is.null(add.parsi.tree)) {
    parsi_tree <- add.parsi.tree
  } else {
    parsi_tree <- NULL
  }

  if (branch.supports & layout != "circular") {
    if (!is.null(add.raxml.tree) |
        !is.null(add.parsi.tree)) {
      # Add support values from different phylogenetic estimation methods
      tree <- .phyloSupports(tree, raxml_tree, parsi_tree)
    }
  }

  if (show.tip.label) {
    # Update specific tip names
    if (!is.null(replace.taxa)) {
      for (i in seq_along(replace.taxa)) {
        tree@phylo$tip.label <- gsub(names(replace.taxa)[i],
                                     replace.taxa[i], tree@phylo$tip.label)
      }
    }
  }

  if (!any(names(tree@data) %in% "prob")) {
    names(tree@data)[1] = "prob"
  }

  #-----------------------------------------------------------------------------
  # Plotting the figure

  tree@phylo$tip.label <- gsub("_", " ", tree@phylo$tip.label)

  tiplabels <- tree@phylo$tip.label

  tree_plot <- ggtree(tree, branch.length = "none", layout = layout,
                      ladderize = TRUE, size = branch.width, ...)

  if (is.null(xlim.tree)) {
    xlim.tree <- max(tree_plot$data$x)+5
  }

  tree_plot <- tree_plot +
    xlim_tree(xlim.tree)

  if (!is.null(hexpand)) {
    tree_plot <- tree_plot +
      ggtree::hexpand(hexpand)
  }

  # Highlight clades
  if (!is.null(highlight.clade)) {
    if (inherits(highlight.clade, "character")) {
      highlight.clade <- list(highlight.clade)
    }
    for (i in seq_along(highlight.clade)) {

      highlight.clade[[i]] <- gsub("_", " ", highlight.clade[[i]])

      tree_plot <- tree_plot +
        geom_hilight(node = ggtree::MRCA(tree, highlight.clade[[i]]),
                     fill = fill.gradient[i], gradient = TRUE,
                     gradient.direction = "tr", alpha = 0.8) +
        geom_tree(layout = layout, size = branch.width)
    }
  }

  # Add and or change tip labels
  if (show.tip.label) {

    # Color highlight specific tips
    if (fancy.tip.label) {
      # No color highlight specific tip labels and use different sizes for taxa
      # and associated voucher information
      tipdata <- .fancy.tiplabels(tree_plot = tree_plot,
                                  size.tip.label = size.tip.label,
                                  highlight.taxa = highlight.taxa,
                                  highlight.color = highlight.color,
                                  understate.taxa = understate.taxa,
                                  abbrev.tip.label = abbrev.tip.label)
      if (layout == "circular") {
        tree_plot <- tree_plot %<+% tipdata +
          ggtext::geom_richtext(data = ggtree::td_filter(isTip),
                                aes(angle = angle, label = lab),
                                fill = NA, label.color = NA,
                                hjust = 0, nudge_x = -0.08)
      } else {
        tree_plot <- tree_plot %<+% tipdata +
          ggtext::geom_richtext(data = ggtree::td_filter(isTip),
                                aes(label = lab),
                                fill = NA, label.color = NA,
                                hjust = 0, nudge_x = -0.08)
      }

    } else if (is.null(understate.taxa) & is.null(highlight.taxa)) {
      # No color highlight for tip labels
      tree_plot <- tree_plot +
        geom_tiplab(offset = 0.02, size = size.tip.label,
                    fontface = fontface.tip.label,
                    show.legend = FALSE, ...)

    } else if (!is.null(highlight.taxa) & is.null(understate.taxa)) {
      # Color highlight only specific tip labels
      tipdata <- .col.tiplabels(tiplabels = tiplabels,
                                highlight.taxa = highlight.taxa,
                                highlight.color = highlight.color)
      cols <- c("black", highlight.color)

      tree_plot <- .plot_with_tiplabels(tree_plot, layout, tipdata,
                                        size.tip.label,
                                        fontface.tip.label,
                                        cols,
                                        ...)

    } else if (is.null(highlight.taxa) & !is.null(understate.taxa)) {
      # Color understate specific tip labels
      tipdata <- .col.tiplabels(tiplabels = tiplabels,
                                understate.taxa = understate.taxa)
      cols <- c("black", "gray60")

      tree_plot <- .plot_with_tiplabels(tree_plot, layout, tipdata,
                                        size.tip.label,
                                        fontface.tip.label,
                                        cols,
                                        ...)

    } else if (!is.null(highlight.taxa) & !is.null(understate.taxa)) {
      # Color highlight and understate specific tip labels
      tipdata <- .col.tiplabels(tiplabels = tiplabels,
                                highlight.taxa = highlight.taxa,
                                highlight.color = highlight.color,
                                understate.taxa = understate.taxa)
      cols <- c("black", highlight.color, "gray60")

      tree_plot <- .plot_with_tiplabels(tree_plot, layout, tipdata,
                                        size.tip.label,
                                        fontface.tip.label,
                                        cols,
                                        ...)
    }
  }

  # Add support values below branches
  if (branch.supports) {

    # Add supports as numbers on the phylogeny
    if (layout != "circular") {
      if (!is.null(add.raxml.tree) |
          !is.null(add.parsi.tree)) {
        if (!any(intree@data$prob > 1)) {
          tree_plot <- tree_plot +
            geom_text2(aes(subset = !isTip & gsub("[/].*", "", prob) >= 0.5,
                           label = prob),
                       size = 3, color = "gray40", hjust = 1.1, vjust = 1.5)
        } else {
          tree_plot <- tree_plot +
            geom_text2(aes(subset = !isTip & gsub("[/].*", "", prob) >= 50,
                           label = prob),
                       size = 3, color = "gray40", hjust = 1.1, vjust = 1.5)
        }

      } else {
        if (!any(intree@data$prob > 1)) {
          tree_plot <- tree_plot +
            geom_text2(aes(subset = !isTip & prob>=0.5 & prob<1,
                           label = stringr::str_extract(prob, "[[:digit:]][.][[:digit:]][[:digit:]]")),
                       size = 3, color = "gray40", hjust = 1.5, vjust = 1.5)
        } else {
          tree_plot <- tree_plot +
            geom_text2(aes(subset = !isTip & prob>=50 & prob<100,
                           label = prob),
                       size = 3, color = "gray40", hjust = 1.5, vjust = 1.5)
        }
      }
    }

    # Add supports as symbols on the phylogeny
    if (!any(intree@data$prob > 1)) {
      tree_plot <- tree_plot +
        # Add point when support for node is greater than 0.95 pp
        geom_point2(aes(subset = prob>=1 & isTip == FALSE), size = 3,
                    shape = ifelse(layout != "circular", 22, 21),
                    fill = "black", alpha = 0.8, stroke = 0.05) +
        geom_point2(aes(subset = prob>=0.9 & prob<1), size = 3,
                    shape = ifelse(layout != "circular", 22, 21),
                    fill = "#0072b2", alpha = 0.8, stroke = 0.05) +
        geom_point2(aes(subset = prob>=0.8 & prob<0.9), size = 3,
                    shape = ifelse(layout != "circular", 22, 21),
                    fill = "#e69f00", alpha = 0.8, stroke = 0.05) +
        geom_point2(aes(subset = prob>=0.7 & prob<0.8), size = 3,
                    shape = ifelse(layout != "circular", 22, 21),
                    fill = "#009e73", alpha = 0.8, stroke = 0.05)
    } else {
      tree_plot <- tree_plot +
        # Add point when support for node is greater than 0.95 pp
        geom_point2(aes(subset = prob>=100 & isTip == FALSE), size = 3,
                    shape = ifelse(layout != "circular", 22, 21),
                    fill = "black", alpha = 0.8, stroke = 0.05) +
        geom_point2(aes(subset = prob>=90 & prob<100), size = 3,
                    shape = ifelse(layout != "circular", 22, 21),
                    fill = "#0072b2", alpha = 0.8, stroke = 0.05) +
        geom_point2(aes(subset = prob>=80 & prob<90), size = 3,
                    shape = ifelse(layout != "circular", 22, 21),
                    fill = "#e69f00", alpha = 0.8, stroke = 0.05) +
        geom_point2(aes(subset = prob>=70 & prob<80), size = 3,
                    shape = ifelse(layout != "circular", 22, 21),
                    fill = "#009e73", alpha = 0.8, stroke = 0.05)
    }

  }

  # Add gene labels
  if (!is.null(gene.label) & layout != "circular") {

    gene.label <- paste0(gene.label, "\n", collapse = "")

    # Control for position of gene labels
    if (is.null(ylim.gene.label)) {
      ylim.gene.label <- (max(tree_plot$data$y)-(max(tree_plot$data$y)*25/100))-10
    }

    tree_plot <- tree_plot +
      annotate(geom = "text", x = 1, y = ylim.gene.label, label = gene.label,
               color = "gray60", size = size.gene.label, fontface = "italic")
  }

  # Create phylogram
  if (phylogram.side & layout != "circular") {
    phylogram <- ggtree(tree, layout = layout, ladderize=TRUE, size = 0.1) +
      theme(
        panel.background = element_rect(fill = "transparent", color = NA),
        plot.background = element_rect(fill = "transparent", color = NA),
      )

    if (!is.null(highlight.clade)) {
      for (i in seq_along(highlight.clade)) {
        phylogram <- phylogram +
          geom_hilight(node = ggtree::MRCA(tree, highlight.clade[[i]]),
                       fill = fill.gradient[i], gradient = TRUE,
                       gradient.direction = "tr", alpha = 0.8,
                       size = 0.1) +
          geom_tree(layout = layout, size = 0.1)
      }
    }

    # Add supports as symbols on the phylogram
    if (phylogram.supports) {
      if (!any(intree@data$prob > 1)) {
        phylogram <- phylogram +
          geom_point2(aes(subset = prob>=1 & isTip == FALSE), size = 1.5, shape = 22,
                      fill = "black", alpha = 0.8, stroke = 0.05) +
          geom_point2(aes(subset = prob>=0.9 & prob<1), size = 1.5, shape = 22,
                      fill = "#0072b2", alpha=0.8, stroke = 0.05) +
          geom_point2(aes(subset = prob>=0.8 & prob<0.9), size = 1.5, shape = 22,
                      fill = "#e69f00", alpha=0.8, stroke = 0.05) +
          geom_point2(aes(subset = prob>=0.7 & prob<0.8), size = 1.5, shape = 22,
                      fill = "#009e73", alpha = 0.8, stroke = 0.05)
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

    phylogram <- phylogram +
      geom_treescale(x = 0, y = max(phylogram$data$y)/1.1, color = "gray60",
                     fontsize = 3, linesize = 0.5, offset = 1)

    if (is.null(phylogram.height)) phylogram.height = max(tree_plot$data$y)*30/100

    if (abbrev.tip.label & fancy.tip.label == FALSE) {
      temp <- abbrevGen(tiplabels = fulltree$data$label[fulltree$data$isTip])
      fulltree$data$label[fulltree$data$isTip] <- temp$abbrev_tiplabels
    }

    # Remove duplicated tree layers
    tree_plot <- .remove_tree_layers(tree_plot)
    phylogram <- .remove_tree_layers(phylogram)

    fulltree <- tree_plot + ggplot2::annotation_custom(
      ggplot2::ggplotGrob(phylogram),
      xmin = -1,
      xmax = max(tree_plot$data$x)-(max(tree_plot$data$x)*70/100),
      ymin = max(tree_plot$data$y)-(max(tree_plot$data$y)*phylogram.height/100),
      ymax = max(tree_plot$data$y)-(max(tree_plot$data$y)*1/100))

  } else {
    # Remove duplicated tree layers
    tree_plot <- .remove_tree_layers(tree_plot)
    fulltree <- tree_plot
    if (abbrev.tip.label) {
      temp <- abbrevGen(tiplabels = fulltree$data$label[fulltree$data$isTip])
      fulltree$data$label[fulltree$data$isTip] <- temp$abbrev_tiplabels
    }
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


#-------------------------------------------------------------------------------
# Auxiliary function to remove duplicated tree layers

.remove_tree_layers <- function(tree_plot_obj) {
  tf <- grepl("stat_tree$|stat_tree[.][.][.]", names(tree_plot_obj[["layers"]]))
  if (length(which(tf)) >= 4) {
    layers <- names(tree_plot_obj[["layers"]])[tf]
    n_layers <- length(layers)
    layers_to_del <- head(layers, n_layers-2)
    for (i in seq_along(layers_to_del)) {
      tree_plot_obj[["layers"]][[layers_to_del[i]]] <- NULL
    }
  }
  return(tree_plot_obj)
}


#-------------------------------------------------------------------------------
# Auxiliary function to get color highlights for different tip labels
.col.tiplabels <- function(tiplabels = NULL,
                           highlight.taxa = NULL,
                           highlight.color = NULL,
                           understate.taxa = NULL) {

  if (!is.null(highlight.taxa)) {
    highlight.taxa <- gsub("_", " ", highlight.taxa)
    tf <- lapply(highlight.taxa, function(x) grepl(paste0(x, collapse = "|"),
                                                   tiplabels))
    tipdata_high <- data.frame(tiplabels = tiplabels,
                               tocolor = NA)
    for (i in seq_along(tf)) {
      tipdata_high$tocolor[tf[[i]]] <- highlight.color[i]
    }
    tipdata_high$tocolor[is.na(tipdata_high$tocolor)] <- "black"
  }

  if (!is.null(understate.taxa)) {
    understate.taxa <- gsub("_", " ", understate.taxa)
    tf <- grepl(paste0(understate.taxa, collapse = "|"), tiplabels)
    tipdata_under <- data.frame(tiplabels = tiplabels,
                                tocolor = NA)
    tipdata_under$tocolor[tf] <- "gray60"
    tipdata_under$tocolor[is.na(tipdata_under$tocolor)] <- "black"
  }

  if (!is.null(highlight.taxa) & is.null(understate.taxa)) {
    tipdata <- tipdata_high
  } else if (is.null(highlight.taxa) & !is.null(understate.taxa)) {
    tipdata <- tipdata_under
  } else if (!is.null(highlight.taxa) & !is.null(understate.taxa)) {
    tipdata_high$tocolor[tf] <- "gray60"
    tipdata <- tipdata_high
  }

  return(tipdata)
}


#-------------------------------------------------------------------------------
# Auxiliary function to plot tree with size and colors of tip labels
.plot_with_tiplabels <- function(tree_plot,
                                 layout,
                                 tipdata,
                                 size.tip.label,
                                 fontface.tip.label,
                                 cols,
                                 ...) {

  if (layout == "circular") {
    tree_plot <- tree_plot %<+% tipdata +
      geom_tiplab(offset = 0.02, size = size.tip.label,
                  fontface = fontface.tip.label,
                  aes(angle = angle, color = tocolor),
                  show.legend = FALSE, ...) +
      scale_colour_manual(values = cols,
                          breaks = cols, guide = "none")
  } else {
    tree_plot <- tree_plot %<+% tipdata +
      geom_tiplab(offset = 0.02, size = 3,
                  fontface = fontface.tip.label, aes(color = tocolor),
                  show.legend = FALSE, ...) +
      scale_colour_manual(values = cols,
                          breaks = cols, guide = "none")
  }
  return(tree_plot)
}


#-------------------------------------------------------------------------------
# Auxiliary function to plot tip labels in Markdown or HTML format
.fancy.tiplabels <- function (tree_plot = NULL,
                              size.tip.label = NULL,
                              highlight.taxa = NULL,
                              highlight.color = NULL,
                              understate.taxa = NULL,
                              abbrev.tip.label = FALSE) {

  tiplabels <- tree_plot$data$label[tree_plot$data$isTip]

  df <- splitTips(tiplabels = tiplabels)

  df$size <- paste0(size.tip.label, "pt")
  df$sizevoucher <- paste0(size.tip.label*(75/100), "pt")
  df$color <- "black"

  if (!is.null(highlight.taxa) & is.null(understate.taxa)) {
  tipcols <- .col.tiplabels(tiplabels = tiplabels,
                            highlight.taxa = highlight.taxa,
                            highlight.color = highlight.color)
  df$color <- tipcols$tocolor
  } else if (is.null(highlight.taxa) & !is.null(understate.taxa)) {
    tipcols <- .col.tiplabels(tiplabels = tiplabels,
                              understate.taxa = understate.taxa)
    df$color <- tipcols$tocolor
  } else if (!is.null(highlight.taxa) & !is.null(understate.taxa)) {
    tipcols <- .col.tiplabels(tiplabels = tiplabels,
                              highlight.taxa = highlight.taxa,
                              highlight.color = highlight.color,
                              understate.taxa = understate.taxa)
    df$color <- tipcols$tocolor
  }

  if (abbrev.tip.label) {
    temp <- abbrevGen(tiplabels = df$genus)
    df$genus <- temp$abbreviation
  }

  if ("doubtID" %in% names(df)) {
    df$doubtID <- gsub("aff", "aff.", df$doubtID)
    df$doubtID <- gsub("cf", "cf.", df$doubtID)
  }
  tf <- df$species %in% "sp" | grepl("^sp[[:upper:]]$|^sp[[:digit:]]", df$species)
  if (any(tf)) {
    df$species[tf] <- gsub("^sp", "sp.", df$species[tf])
  }

  # https://yulab-smu.top/treedata-book/faq.html

  ntaxa <- names(df)[names(df) %in% c("genus", "doubtID", "species", "infrasp")]
  nvouc <- names(df)[names(df) %in% c("voucher", "genbank")]

  # Example of the Markdown or HTML to format text
  # glue("<i style='font-size:15pt'>**{genus} {species}**</i> <span style='font-size:10pt'>{voucher}</span>")

  ntaxa <- paste0("{", ntaxa, "}", collapse = " ")
  nvouc <- paste0("{", nvouc, "}", collapse = " ")
  ntaxa <- paste0("<span style='color:{color}'> <i style='font-size:{size}'>**",
                  ntaxa,
                  "**</i> <span", collapse = "")
  nvouc <- paste0("style='font-size:{sizevoucher}'>",
                  nvouc,
                  "</span> </span>", collapse = "")

  exp <- paste(ntaxa, nvouc)

  df2 <- dplyr::mutate(df,
                       lab = glue::glue(exp))

  df2$lab <- gsub("\\sNA\\s", " ", df2$lab)
  df2$lab <- gsub("[>]NA\\s", ">", df2$lab)
  df2$lab <- gsub(" NA[*][*]", "**", df2$lab)
  df2$lab <- gsub("\\sNA[<]", "<", df2$lab)

  tf <- grepl("[>]NA[<]", df2$lab)
  if (any(tf)) {
    df2$lab[tf] <- gsub("\\s[<]span\\s.*", "", df2$lab[tf])
  }

  return(df2)
}

# Another way to plot multiple fonts but without multiple sizes
# exp <- ifelse(!n %in% c("voucher", "accession"),
#               paste0("bolditalic({", n, "})"),
#               paste0("{", n, "}"))
# exp <- paste0(exp, collapse = "~")
# d2 <- dplyr::mutate(d,
#                     lab = glue::glue("bolditalic({species})~{accession}"))
# d2 <- dplyr::mutate(d,
#                     lab = glue::glue(exp))
# d2$lab <- gsub("[~]bolditalic[(]NA[)]|[~]NA$", "", d2$lab)

# p <- ggtree(Harpalyce_bayes_tree) %<+% d2 +
#   geom_tiplab(aes(label=lab), parse=T) +
#   hexpand(1)

