library(VennDiagram)

# A more complicated diagram
venn.plot <- draw.triple.venn(
  area1 = 265,
  area2 = 322,
  area3 = 78,
  n12 = 22,
  n23 = 11,
  n13 = 20,
  n123 = 4,
  category = c("FO-GFP", "T-GFP", "MZ-GFP"),
  fill = c("blue", "yellow", "green"), # fill = c("#7f7fff", "#ffff7f", "#7fff7f"),
  col = rep("black", 3), cat.pos=c(-20, 20, 180), cat.dist=c(0.05, 0.05, 0.035),
  cat.col = rep("black", 3), 
  cex = 2.5, alpha = rep(0.5, 3),
  cat.cex = 2.5, fontfamily=rep("sans",7),
  cat.fontfamily=rep("sans",3),
  fontface=rep("bold",7),
  cat.fontface=rep("bold",3),
  ind = FALSE)
grid.newpage()
grid.draw(venn.plot)

# Writing to file
# tiff(filename = "venn.tiff")
# grid.draw(venn.plot)
# dev.off()

pdf(file = "./Triple_Venn_diagram.pdf")
grid.draw(venn.plot)
dev.off()
