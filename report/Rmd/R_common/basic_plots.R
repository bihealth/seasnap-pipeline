############################################
## define functions for generic plots
############################################

multiplot <- function(..., plotlist=NULL, cols=1, layout=NULL, fn=NULL, height=7, width=7) {
    plots <- c(list(...), plotlist)
    numPlots <- length(plots)
    if (is.null(layout)) {
        layout <- matrix(seq(1, cols * ceiling(numPlots/cols)), ncol=cols, nrow=ceiling(numPlots/cols))
    }

    if (numPlots == 1) {
        if (!is.null(fn)) {
            pdf(fn, height=height, width=width)
            print(plots[[1]])
            invisible(dev.off())
        }
        print(plots[[1]])
    } else {
        if (!is.null(fn)) {
            pdf(fn, height=height, width=width)
            for (p in plots) print(p)
            invisible(dev.off())
        }
        grid::grid.newpage()
        grid::pushViewport(grid::viewport(layout=grid::grid.layout(nrow(layout), ncol(layout))))
        for (i in 1:numPlots) {
            matchidx <- as.data.frame(which(layout==i, arr.ind=TRUE))
            print(plots[[i]], vp=grid::viewport(layout.pos.row=matchidx$row, layout.pos.col=matchidx$col))
        }
    }
}
