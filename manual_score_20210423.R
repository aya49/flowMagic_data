avgF1 <- c(.92,.76,.98,.81,.95,.90,.68,.81,.94,.85)
scat=c("01_CD66CD45_Leukocyte", "02_CD3CD19_mononuclear_", "03_CD14CD7_NKLinNeg_", "04_CD14CD16_lin-_", "05_CD4CD8_Tcell_", "06_CD4CD45RA_Tcell_", "07_FoxP3CD25_CD4Tcell", "08_TCRgdCD3_CD4Tcell", "09_CD8CD45RA_NotCD4CD8Tcell_", "10_TbetCD45RA_CD8Tcell_")

shots <- data.frame(
    F1=c(c(.80,.32,.90,.41,.53,.86,.38,.77,.75,.65), 
         c(.88,.77,.95,.68,.89,.93,.42,.80,.94,.81), avgF1, 
         c(.92,.74,.98,.83,.96,.92,.73,.81,.97,.88),
         c(.96,.77,.98,.80,.97,.92,.75,.88,.96,.89)),
    scatterplot=rep(scat,5),
    shot=c(rep(1,10), rep(5,10), rep(10,10), rep(15,10), rep(20,10))
)

baseline <- data.frame(
    F1=  c(.97,.71,.97,.84,.98,.94,.73,.80,.95,.80),
    scatterplot=scat
)

ggplot2::ggplot(shots, ggplot2::aes(x=shot, y=F1)) +
    ggplot2::geom_line(ggplot2::aes(colour=scatterplot))
