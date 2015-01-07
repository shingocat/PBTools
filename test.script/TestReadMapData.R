###############################################################################
# TODO: Add comment
# 
# Author: MQin
# Date: Jan 6, 2015
# FileName: TestReadMapData.R
###############################################################################


map <- read.map.data("E:\\Data\\map_demo.csv", marker ="Marker", chr ="Chromosome", pos = "Position")

map2 <- read.map.data("E:\\Demo\\SelectionTools-examples\\il.demo.map", header = F)