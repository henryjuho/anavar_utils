library(ggplot2)

sfs <- round(100/(1:19))
print(sfs)

sfs_data <- as.data.frame(cbind(1:19, sfs))
colnames(sfs_data) <- c('freq', 'count')

plot = ggplot(sfs_data, aes(x=as.factor(freq), y=count)) +
    geom_bar(stat='identity', position='dodge') +
    theme_bw() +
    xlab('Derived allele frequency')  +
    ylab("Number of variants")

png('sfs_example.png', width=6, height=3, units='in', res=320)

plot

dev.off()
