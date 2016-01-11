args <- commandArgs(trailingOnly = TRUE)
print(args)

printFile <- function(filename, output_filename)
{
  len <- 50;
  
  duplication_rate = read.table(filename, header=T, sep="\t") 
  
  plot_colors = c("black","green","red")
  line_types = c("solid","solid","solid");
  
  #cairo_pdf(output_filename,width=9,height=7.5)
png(output_filename, width=800, height=600)

  par(mar=c(5, 4, 4, 6) + 0.1)
  max_y = max(duplication_rate)
  plot(duplication_rate$rate[c(1:len)], duplication_rate$non_unique[c(1:len)]+duplication_rate$unique[c(1:len)], log="y",type="l", lty=line_types[1], col=plot_colors[1], axes=FALSE, ann=FALSE)
  axis(1, at=2*(1:len/2), labels=2*(1:len/2))
  
  aty <- axTicks(2,log=TRUE)
  labels <- sapply(log(aty,10),function(i)
    as.expression(bquote(10^.(i))))
  axis(2, at=axTicks(2,log=TRUE), labels=format(aty, scientific=TRUE))
  
  lines(duplication_rate$rate[c(1:len)],duplication_rate$unique[c(1:len)],type="l", lty=line_types[2], col=plot_colors[2])
  title(xlab= "Duplication Rate")
  title(ylab= "Number of Reads")
  box()
  
  
  pcr_artifacts_percentage <-  duplication_rate$non_unique / (duplication_rate$unique + duplication_rate$non_unique)
  pcr_artifacts_absolute <-  duplication_rate$non_unique
  
  
  max_y_difference = max(pcr_artifacts_percentage)
  par(new=TRUE)
  plot(duplication_rate$rate[c(1:len)],pcr_artifacts_percentage[c(1:len)], type="l", lty=line_types[3], col=plot_colors[3],ann=F, axes=F)
  mtext("Percentage PCR artifacts",side=4,col="black",line=4) 
  usr<-par("usr")
  
  axis(4, at=axTicks(4), labels=format(axTicks(4), scientific=FALSE),las=1)
  
  title(main="Reads per Duplication Rate", col.main="black", font.main=4)
  
  print(usr)
  
  text((usr[1]+usr[2])/1.6,usr[4]*0.99, paste("mapped reads: ",as.character(sum(duplication_rate$unique)+sum(duplication_rate$non_unique)),"\n",
                     "unique reads: ",as.character(sum(duplication_rate$unique)),"\n",
                     "non-unique reads: ",as.character(sum(duplication_rate$non_unique)),"\n",
                     "PCR artifacts rate: ",as.character(signif(sum(duplication_rate$non_unique)/(sum(duplication_rate$non_unique)+sum(duplication_rate$unique))),3)), adj=c(1,1))
  
  legend("topright", c("total duplicates","unique duplicates", "PCR artifacts rate"), cex=0.8, col=plot_colors, lty=line_types)
  dev.off()
  
}

printFile(args[1], args[2])

