# Tomado de: http://www.sthda.com/english/wiki/a-complete-guide-to-3d-visualization-device-system-in-r-r-software-and-data-visualization

rgl_init <- function(new.device = FALSE, bg = "white", width = 640) { 
  if( new.device | rgl.cur() == 0 ) {
    rgl.open()
    par3d(windowRect = 50 + c( 0, 0, width, width ) )
    rgl.bg(color = bg )
  }
  rgl.clear(type = c("shapes", "bboxdeco"))
  rgl.viewpoint(theta = 15, phi = 20, zoom = 0.7)
}

rgl_add_axes <- function(x, y, z, axis.col = "grey",
                         xlab = "", ylab="", zlab="", show.plane = TRUE, 
                         show.bbox = FALSE, bbox.col = c("#333377","black"))
{ 
  
  lim <- function(x){c(-max(abs(x)), max(abs(x))) * 1.1}
  # Add axes
  xlim <- lim(x); ylim <- lim(y); zlim <- lim(z)
  rgl.lines(xlim, c(0, 0), c(0, 0), color = axis.col)
  rgl.lines(c(0, 0), ylim, c(0, 0), color = axis.col)
  rgl.lines(c(0, 0), c(0, 0), zlim, color = axis.col)
  
  # Add a point at the end of each axes to specify the direction
  axes <- rbind(c(xlim[2], 0, 0), c(0, ylim[2], 0), 
                c(0, 0, zlim[2]))
  rgl.points(axes, color = axis.col, size = 3)
  
  # Add axis labels
  rgl.texts(axes, text = c(xlab, ylab, zlab), color = axis.col,
            adj = c(0.5, -0.8), size = 2)
  
  # Add plane
  if(show.plane) 
    xlim <- xlim/1.1; zlim <- zlim /1.1
  rgl.quads( x = rep(xlim, each = 2), y = c(0, 0, 0, 0),
             z = c(zlim[1], zlim[2], zlim[2], zlim[1]))
  
  # Add bounding box decoration
  if(show.bbox){
    rgl.bbox(color=c(bbox.col[1],bbox.col[2]), alpha = 0.5, 
             emission=bbox.col[1], specular=bbox.col[1], shininess=5, 
             xlen = 3, ylen = 3, zlen = 3) 
  }
}

pan3d <- function(button, dev = rgl.cur(), subscene = currentSubscene3d(dev)) {
  start <- list()
  
  begin <- function(x, y) {
    activeSubscene <- par3d("activeSubscene", dev = dev)
    start$listeners <<- par3d("listeners", dev = dev, subscene = activeSubscene)
    for (sub in start$listeners) {
      init <- par3d(c("userProjection","viewport"), dev = dev, subscene = sub)
      init$pos <- c(x/init$viewport[3], 1 - y/init$viewport[4], 0.5)
      start[[as.character(sub)]] <<- init
    }
  }
  
  update <- function(x, y) {
    for (sub in start$listeners) {
      init <- start[[as.character(sub)]]
      xlat <- 2*(c(x/init$viewport[3], 1 - y/init$viewport[4], 0.5) - init$pos)
      mouseMatrix <- translationMatrix(xlat[1], xlat[2], xlat[3])
      par3d(userProjection = mouseMatrix %*% init$userProjection, dev = dev, subscene = sub )
    }
  }
  rgl.setMouseCallbacks(button, begin, update, dev = dev, subscene = subscene)
  cat("Callbacks set on button", button, "of rgl device", dev, "in subscene", subscene, "\n")
}

rgl_init()
rgl.spheres(mc_orloci[,1], mc_orloci[,2], mc_orloci[,3], r = 0.05, color = "yellow")
rgl.bbox(color = "grey90")
rgl_add_axes(mc_orloci[,1], mc_orloci[,2], mc_orloci[,3], show.bbox = TRUE)
texts3d(mc_orloci[,1], mc_orloci[,2], mc_orloci[,3],
        pos = 3, offset = 0.5,
        texts = rownames(mc_orloci), cex = 2)


rgl_init()
rgl.spheres(mc_orloci_norm[,1], mc_orloci_norm[,2], mc_orloci_norm[,3], r = 0.05, color = "yellow")
rgl.bbox(color = "grey90")
rgl_add_axes(mc_orloci_norm[,1], mc_orloci_norm[,2], mc_orloci_norm[,3], show.bbox = TRUE)
texts3d(mc_orloci_norm[,1], mc_orloci_norm[,2], mc_orloci_norm[,3],
        pos = 3, offset = 0.5,
        texts = rownames(mc_orloci), cex = 2)
