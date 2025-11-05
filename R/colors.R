#' Color palette: GAP_Rainbow
#'
#' A color palette with 13 colors in rainbow used for continuous data.
#'
#' @export
GAP_Rainbow <- c('#00008C', '#004DEB', '#007359', '#66A600', '#B3E600',
                 '#CCF200', '#FFFF00', '#FFE600', '#FFB300', '#FF6600',
                 '#FF0000', '#D90000', '#A60000')


#' Color palette: GAP_Blue_White_Red
#'
#' A diverging color palette from blue to white to red, suitable for visualizing correlations.
#'
#' @export
GAP_Blue_White_Red <- c('#00006B', '#00009E', '#0000BF', '#0000E8', '#2121FF','#4A4AFF',
                        '#7A7AFF', '#ADADFF', '#E0E0FF', '#FAFAFF', '#FFFAFA', '#FFC7C7',
                        '#FF9494', '#FF6161', '#FF3838', '#FF0000', '#CF0000', '#AD0000',
                        '#850000', '#520000')


#' Color palette: GAP_d
#'
#' A 16-category discrete color palette.
#'
#' @export
GAP_d <- c(
  rgb(1.0, 0.0, 0.0), # Red
  rgb(0.0, 0.5, 0.0), # Green
  rgb(0.0, 0.0, 1.0), # Blue
  rgb(1.0, 1.0, 0.0), # Yellow
  rgb(0.0, 1.0, 1.0), # Cyan
  rgb(1.0, 0.0, 1.0), # Magenta
  rgb(0.5, 0.5, 0.5), # Gray
  rgb(1.0, 0.5, 0.0), # Orange
  rgb(0.5, 0.0, 0.5), # Purple
  rgb(0.5, 1.0, 0.0), # Light Green
  rgb(0.5, 0.25, 0.25), # Brownish
  rgb(0.5, 0.5, 0.0), # Olive
  rgb(0.25, 0.5, 0.5), # Teal
  rgb(1.0, 0.5, 0.5), # Light Red / Salmon
  rgb(0.0, 0.0, 0.5), # Dark Blue
  rgb(1.0, 1.0, 1.0) # White
)


#' Color palette: grayscale_palette
#'
#' A simple two-color grayscale palette from white to black, often used for binary data visualization.
#'
#' @export
grayscale_palette <- gray(seq(1, 0, length.out = 2)) # 0 - white，1 - black
