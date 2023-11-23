pkgTest <- function(x) {
  if (!require(x, character.only = TRUE)) {
    install.packages(x, dep = TRUE)
    if (!require(x, character.only = TRUE)) {
      BiocManager::install(x, ask = F)
    }
    require(x)
  }
}

pkgTest("ShortRead")
pkgTest("tidyverse")
pkgTest("data.table")

sample_files <- readRDS("rds_data/raw_fastq/sample_files.rds")

################ QC for 1 file ############################

file_name <- sample_files$batch1[1]

reads <- readFastq(file_name)
qaSummary <- qa(file_name, type="fastq")

total_reads <- qaSummary[["readCounts"]]$read

base_quality <- qaSummary[["perCycle"]]$quality
head(base_quality)


means <- base_quality %>% group_by(Cycle) %>% summarise(means = sum(Score*Count)/sum(Count))

get_quant <- function(xx, yy, q) { xx[which(cumsum(yy)/sum(yy) >=q)][[1]] }

q25 <- base_quality %>% group_by(Cycle) %>% summarise(q25 = get_quant(Score, Count, 0.25))

stats_fastq_qual <-
  base_quality %>% group_by(Cycle) %>% summarise(
    means = sum(Score * Count) / sum(Count),
    q10 = get_quant(Score, Count, 0.1),
    q25 = get_quant(Score, Count, 0.25),
    q50 = get_quant(Score, Count, 0.5),
    q75 = get_quant(Score, Count, 0.75),
    q90 = get_quant(Score, Count, 0.9)
  )


ggplot()+
  geom_boxplot(data = stats_fastq_qual, aes(ymin = q10, x = Cycle,
                              lower = q25, middle = q50, upper = q75,
                              ymax = q90, group = Cycle), stat = "identity",
               fill = "gold")

rectangle_color <- list("bad" = list("ymin" = 0, "ymax" = 20),
                        "med" = list("ymin" = 20, "ymax" = 30),
                        "good" = list("ymin" = 30, "ymax" = max(stats_fastq_qual$q90)))
rectangle_color <- rbindlist(lapply(names(rectangle_color), function(x) {
  data.frame(name = x,
             ymin = rectangle_color[[x]]$ymin,
             ymax = rectangle_color[[x]]$ymax,
             xmin = 0,
             xmax = max(stats_fastq_qual$Cycle)*1.05)
}))

ggplot() +
  geom_rect(data = rectangle_color, aes(fill = name,
                                        xmin = xmin,
                                        ymin = ymin,
                                        xmax = xmax,
                                        ymax = ymax),
            alpha = 0.3) +
  geom_boxplot(data = stats_fastq_qual, aes(ymin = q10, x = Cycle,
                              lower = q25, middle = q50, upper = q75,
                              ymax = q90, group = Cycle), stat = "identity",
               fill = "gold")+
  scale_fill_manual(values = c("bad" = "firebrick",
                               "med" = "orange",
                               "good" = "darkgreen"))+
  theme_minimal()+
  guides(fill="none")+
  theme(strip.text = element_text(size = 10),
        panel.spacing=unit(1,"lines"))+
  scale_x_continuous(limits = c(0, max(stats_fastq_qual$Cycle)*1.05), 
                     expand = expansion(mult = c(0, 0),
                                        add = c(3, 0)))+
  geom_line(data = stats_fastq_qual, aes(x = Cycle, y = q50), color = "#800020", linewidth = 2)+
  labs(x = "Position", y = "Quality Score")




ggplot()+
  geom_segment(data = read_len, aes(x = read_length, y = 0, xend = read_length, yend = num_reads))+
  geom_point(data = read_len, aes(x = read_length, y = num_reads), color = "#9370DB")+
  scale_x_continuous(limits = c(0,max(read_len$read_length)*1.1), expand = c(0,0)) +
  scale_y_continuous(limits = c(0,max(read_len$num_reads)*1.1),expand = c(0,0))+
  theme_minimal()+
  guides(fill="none")+
  labs(x = "Read len", y = "Frequency")


