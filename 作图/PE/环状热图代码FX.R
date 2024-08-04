library(data.table)
library(tidyverse)
library(circlize)
library(ComplexHeatmap)
library(dplyr)
setwd("E:\\作图\\正向\\PE") 
data=fread("圈图PE11.txt")
hhh=data %>% as.data.frame()

setwd("E:\\作图\\正向\\PE") 

rownames(hhh)=hhh[,1]

hhh=hhh[,-1]

col_fun = list(
  col_1 = colorRamp2(c(0.8, 1, 1.3), 
                     c("#41ab5b", "white", "#CD0000")),
  col_2 = colorRamp2(c(0.7, 1, 1.5),
                     c("#1874CD", "white", "#CD3700")),
  col_3 = colorRamp2(c(0.8, 1, 1.3),
                     c("#96CDCD", "white", "#E9499F")),
  col_4 = colorRamp2(c(0.00009, 0.001, 0.05),
                     c("#000000", "#575757", "#FFFFFF"))
)

#绘图

pdf("plot.pdf", height = 12.5, width = 12.5) #设置输出的pdf文件长宽大小
circos.par$gap.degree <- 55 #相邻扇区的角度度数
circos.par$start.degree <- 30 #扇区从圆的顶部开始的几度角度处开始
circos.par$track.margin <- c(0.0001, 0.0001) #设置环形图中轨道之间的边距（间距）

for (i in 1:4) { #有几个结果就需要把i循环的上限改成相应的数字加1
  data_tmp <- as.matrix(hhh[,i]) #提取相应的结果
  if (i == 1) {
    rownames(data_tmp) <- rownames(hhh) #第一圈外围的图例
  }
  colnames(data_tmp) <- colnames(hhh)[i] #列名
  
  if (i <= 4) {
    #开始绘制圆形热图
    circos.heatmap(data_tmp, 
                   col = col_fun[[i]], #提取上述设定的颜色
                   rownames.side = "outside", #设定名字在扇形的外面
                   cluster = T, #分类和分组，能够让后面的class与相应的traits对应
                   cell.border = "white", #单元格之间的边界颜色
                   track.height = 0.06) #热图的轨道高度
  } else {
    circos.heatmap(data_tmp, 
                   col = col_fun[[i]],
                   rownames.side = "outside", 
                   cluster = F,
                   track.height = 0.04) 
  }
}

#设置图例
lgd1 <- Legend(title = colnames(hhh)[1], #相应结果的标题名字
               border = "black", #设置图例边框的颜色为黑色
               grid_height = unit(3, "mm"), #每个图例项之间的垂直间距为 3 毫米
               legend_width = unit(20, "mm"), #参数设置图例的宽度为 20 毫米
               at = c(0.5, 1.5),  #显示图例中各个cutoff的数值
               title_position = "topcenter", #标题的位置
               col_fun = col_fun[[1]], #颜色的图例
               direction = "horizontal",title_gp = gpar(fontsize = 8, fontface = "bold")) #图例的方向

lgd2 <- Legend(title = colnames(hhh)[2], border = "black", grid_height = unit(3, "mm"),
               legend_width = unit(20, "mm"),
               at = c(0.5, 1.5), 
               title_position = "topcenter",
               col_fun = col_fun[[2]], direction = "horizontal",title_gp = gpar(fontsize = 8, fontface = "bold"))

lgd3 <- Legend(title = colnames(hhh)[3], border = "black", grid_height = unit(3, "mm"),
               legend_width = unit(20, "mm"),
               at =c(0.5, 1.5),  
               title_position = "topcenter",
               col_fun = col_fun[[3]], direction = "horizontal",title_gp = gpar(fontsize = 8, fontface = "bold"))

lgd4 <- Legend(title = colnames(hhh)[4], border = "black", 
               grid_height = unit(3, "mm"),
               legend_width = unit(20, "mm"),
               at = c(0.00009, 0.05),
               title_position = "topcenter",
               col_fun = col_fun[[4]], direction = "horizontal",title_gp = gpar(fontsize = 8, fontface = "bold"))





pd <- packLegend(lgd1, lgd2, lgd3, lgd4, row_gap = unit(1, "mm")) 
draw(pd, x = unit(0.6, "npc"), y = unit(0.75, "npc")) 

dev.off() 
circos.clear() 

