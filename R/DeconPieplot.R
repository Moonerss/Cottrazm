# source('R/spatial_decon_settings.R')
# TumorST <- readr::read_rds("YourPath/TumorBoundary/1.BoundaryDefine/CRC1/TumorSTBoundaryDefine.rds.gz")
# DeconData <- openxlsx::read.xlsx(system.file("extdata/DeconData.xlsx",package = "Cottrazm"))
# plot_col <- colnames(DeconData)[2:ncol(DeconData)]
# img_path = system.file("extdata/outs/spatial/tissue_lowres_image.png",package = "Cottrazm")
# pie_scale = 0.4
# scatterpie_alpha = 0.8
# border_color = "grey"

#' Title Pie plot of deconvolution result
#'
#' Pie plot of deconvolution result
#'
#' @param DeconData A data frame of deconvolution result of each spatial spot
#' @param TumorST A Seurat object with defined malignant (Mal) spots, boundary (Bdy) spots, and non-malignnat (nMal) spots
#' @param plot_col Names of cell types to plot
#' @param img_path Path to tissue_lowers_image.jpg in Spaceranger output file
#' @param pie_scale Size of spot in pie plot
#' @param scatterpie_alpha Transparency of spot in pie plot
#' @param border_color Color of border in pie plot
#'
#' @return A ggplot object each pie represented a ST spot, each part of the pie represent the percentage of cell types infiltrated in the ST spot
#' @export
#'
#' @examples
#' \dontrun{
#'   TumorST <- readr::read_rds("YourPath/TumorSTBoundaryDefine.rds.gz")
#'   DeconData <- openxlsx::read.xlsx(system.file("extdata/DeconData.xlsx", package = "Cottrazm"))
#'   plot_col <- colnames(DeconData)[2:ncol(DeconData)]
#'   img_path <- system.file("extdata/outs/spatial/tissue_lowres_image.png", package = "Cottrazm")
#'   DeconPieplot(DeconData = DeconData, TumorST = TumorST, plot_col = plot_col, img_path = img_path,
#'                pie_scale = 0.4, scatterpie_alpha = 0.8, border_color = "grey")
#' }
#'
DeconPieplot <- function(DeconData = DeconData,
                         TumorST = TumorST,
                         plot_col = plot_col,
                         img_path = img_path,
                         pie_scale = pie_scale,
                         scatterpie_alpha = scatterpie_alpha,
                         border_color = border_color) {

  check_pkgs(c('jpeg', 'png', 'grid', 'scatterpie', 'cowplot'))

  DeconData <- DeconData[DeconData$cell_ID %in% rownames(TumorST@meta.data), ]

  ## Preprocess data
  slice <- names(TumorST@images)[1]

  spatial_coord <- data.frame(TumorST@images[[slice]]@coordinates) %>%
    tibble::rownames_to_column("cell_ID") %>%
    dplyr::mutate(
      imagerow_scaled =
        .data$imagerow * TumorST@images[[slice]]@scale.factors$lowres,
      imagecol_scaled =
        .data$imagecol * TumorST@images[[slice]]@scale.factors$lowres
    ) %>%
    dplyr::inner_join(DeconData, by = "cell_ID")

  ### Load histological image into R
  img_path <- img_path # lowers image png(input dir)
  img_frmt <- tolower(stringr::str_sub(img_path, -4, -1))

  if (img_frmt %in% c(".jpg", "jpeg")) {
    img <- jpeg::readJPEG(img_path)
  } else if (img_frmt == ".png") {
    img <- png::readPNG(img_path)
  }

  # Convert image to grob object
  img_grob <- grid::rasterGrob(img,
    interpolate = FALSE,
    width = grid::unit(1, "npc"),
    height = grid::unit(1, "npc")
  )


  ## Plot spatial scatterpie plot
  scatterpie_plt <- ggplot() +
    annotation_custom(
      grob = img_grob,
      xmin = 0,
      xmax = ncol(img),
      ymin = 0,
      ymax = -nrow(img)
    ) +
    scatterpie::geom_scatterpie(
      data = spatial_coord,
      aes(
        x = .data$imagecol_scaled,
        y = .data$imagerow_scaled
      ),
      cols = plot_col,
      color = border_color,
      alpha = scatterpie_alpha,
      pie_scale = pie_scale,
      lwd=0.1
    ) +
    scale_y_reverse() +
    ylim(nrow(img), 0) +
    xlim(0, ncol(img)) +
    cowplot::theme_half_open(11, rel_small = 1) +
    theme_void() +
    coord_fixed(
      ratio = 1,
      xlim = NULL,
      ylim = NULL,
      expand = TRUE,
      clip = "on"
    )
  return(scatterpie_plt)
}
