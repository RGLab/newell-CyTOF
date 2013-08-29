#' Copies channel names from a source flowFrame to another flowFrame
#'
#' The flowFrame objects must have the same number of channels.
#'
#' @param flow_frame the \code{flowFrame} object that will be updated
#' @param source_flowframe the \code{flowFrame} object from which the channel
#' names are borrowed
#' @return an updated version of \code{flow_frame}
copy_flowframe_channels <- function(flow_frame, source_flowframe) {
  if (ncol(flow_frame) != ncol(source_flowframe)) {
    stop("Cannot proceed. ",
         "The two flowFrames have different numbers of channels.")
  }

  fr_rownames <- rownames(parameters(flow_frame)@data)
  channel_idx <- paste0(fr_rownames, "N")

  for (j in seq_along(channel_idx)) {
    channel_j <- channel_idx[j]
    flow_frame@parameters@data$name[j] <- source_flowframe@parameters@data$name[j]
    keyword(flow_frame)[[channel_j]] <- keyword(source_flowframe)[[channel_j]]
    flow_frame@description[[channel_j]] <- source_flowframe@description[[channel_j]]
  }

  colnames(exprs(flow_frame)) <- colnames(flow_frame)
  
  flow_frame
}
