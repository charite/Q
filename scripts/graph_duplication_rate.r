# Read car and truck values from tab-delimited autos.dat
#args <- commandArgs(trailingOnly = TRUE)
#print(args)

printFile <- function(filename, output_filename)
{
  len <- 50;
  
  duplication_rate = read.table(filename, header=T, sep="\t") 
  
  # Compute the largest y value used in the data (or we could
  # just use range again)
  
  # Define colors to be used for cars, trucks, suvs
  plot_colors = c("black","black","darkgray")
  line_types = c("dotted","solid","solid");
  
  # Start PNG device driver to save output to figure.png
  cairo_pdf(output_filename,width=9,height=7.5)
  #pdf("duplication_rate.pdf")
  #png("duplication_rate.png", width=800, height=600)
  
  # Graph autos using y axis that ranges from 0 to max_y.
  # Turn off axes and annotations (axis labels) so we can 
  # specify them ourself
  par(mar=c(5, 4, 4, 6) + 0.1)
  max_y = max(duplication_rate)
  plot(duplication_rate$rate[c(1:len)], duplication_rate$non_unique[c(1:len)]+duplication_rate$unique[c(1:len)], log="y",type="l", lty=line_types[1], col=plot_colors[1], axes=FALSE, ann=FALSE)
  axis(1, at=2*(1:len/2), labels=2*(1:len/2))
  
  #ylims <- c(0.2, max_y)
  #get_axp <- function(x) 10^c(ceiling(x[1]), floor(x[2]))
  #usr.i <- log10(ylims)
  #(aty <- axTicks(side = 2, usr = usr.i,
  #axp = c(get_axp(usr.i), n = 3), log = TRUE, nintLog = Inf))
  aty <- axTicks(2,log=TRUE)
  labels <- sapply(log(aty,10),function(i)
    as.expression(bquote(10^.(i))))
  #axis(2, aty, labels = labels)
  axis(2, at=axTicks(2,log=TRUE), labels=format(aty, scientific=TRUE))
  #axis(2, at=axTicks(2,log=TRUE), labels=labels)
  
  # Graph trucks with red dashed line and square points
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
  
  #aty_diff <- seq(0,1,by=0.2)
  #labels_diff <- sapply(aty_diff,function(i)
  #  as.expression(bquote(10^.(i))))
  axis(4, at=axTicks(4), labels=format(axTicks(4), scientific=FALSE),las=1)
  
  #axis(4, ylim=c(min(duplication_rate_difference),max_y_difference))
  title(main="Reads per Duplication Rate", col.main="black", font.main=4)
  
  # Label the x and y axes with dark green text

  print(usr)
  
  text((usr[1]+usr[2])/1.6,usr[4]*0.99, paste("unique reads: ",as.character(sum(duplication_rate$unique)),"\n",
                     "nonunique reads: ",as.character(sum(duplication_rate$non_unique))), adj=c(1,1))
  
  legend("topright", c("total duplicates","unique duplicates", "PCR artifacts rate"), cex=0.8, col=plot_colors, lty=line_types)
  dev.off()
  
}

setwd("P:/promotion/duplication_rate")

filename1 <- "P:\\Promotion\\duplication_rate\\SRR1175698_duplication_rate_reads.txt"
filename2 <- "P:\\Promotion\\duplication_rate\\SRR1175699_duplication_rate_reads.txt"
filename3 <- "P:\\Promotion\\duplication_rate\\SRR1509029_duplication_rate_reads.txt"
filename4 <- "P:\\Promotion\\duplication_rate\\SRR1509030_duplication_rate_reads.txt"

duplication_rate = read.table(filename1, header=T, sep="\t") 

printFile(filename1, "SRR1175698.pdf")
printFile(filename2, "SRR1175699.pdf")
printFile(filename3, "SRR1509029.pdf")
printFile(filename4, "SRR1509030.pdf")

system(paste("pdftops -eps ","SRR1175698.pdf"), FALSE)
system(paste("pdftops -eps ","SRR1175699.pdf"), FALSE)
system(paste("pdftops -eps ","SRR1509029.pdf"), FALSE)
system(paste("pdftops -eps ","SRR1509030.pdf"), FALSE)
# Turn off device driver (to flush output to png)
