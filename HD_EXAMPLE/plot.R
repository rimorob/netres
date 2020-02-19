lin = read.csv('linear.csv')
au = read.csv('autoencoder.csv')

r2Au = cor.test(au$allCAG, au$allCAGpred, method='pearson')$estimate^2
r2Lin = cor.test(lin$allCAG, lin$allCAGpred, method='pearson')$estimate^2

mAu = lm(allCAGpred ~ allCAG, data = au)
mLin = lm(allCAGpred ~ allCAG, data = lin)

library(plugdensity)
cagDens = plugin.density(lin$allCAG, nout=2*nrow(lin))
df = data.frame(x=cagDens$x, y = cagDens$y)
library(ggplot2)
library(ggpubr)

jitter = data.frame(x = lin$allCAG, y = 0)
rr = range(lin$allCAG)
rr[1] = rr[1] - 0.5
rr[2] = rr[2] + 0.5

perf = data.frame(pred=c(lin$allCAGpred, au$allCAGpred),
                  CAG=c(lin$allCAG, lin$allCAg),
                  type=c(rep('linear', nrow(lin)),
                  rep('autoencoder', nrow(au))))

p1 = ggplot(lin, aes(x=lin$allCAG, fill=lin$allCAG)) + geom_density(fill='skyblue', alpha=0.1, bw='SJ-ste') + xlab("\n Log(CAG)") + ylab("Density \n") +
    geom_jitter(data = jitter, aes(x,y), height = 0.01) + ggpubr::theme_pubclean() + guides(fill=FALSE) + xlim(rr)

p2 = ggplot(perf, aes(x=CAG, y = pred, color=type)) + geom_point() + geom_smooth(method=lm, aes(fill=type)) + scale_colour_manual(values = c("skyblue", "orange"), labels = c('autoencoder', 'linear')) +
    scale_fill_manual(values = c("skyblue", "orange"), labels = c('autoencoder', 'linear')) + ggpubr::theme_pubclean() 

pdf('cagPlots.pdf', width=11)
print(p1)
print(p2)
dev.off()

