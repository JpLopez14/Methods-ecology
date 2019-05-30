########################################################
############## Estructura de la comunidad  #############
########################################################

## Se muestras las principales funciones para el analisis ecológico de comunidades ## 

# 1. Calculo de la riqueza, abundancia, diversidad de H', Equidad, Dominancia
# 2. Curvas de rarefacción
# 3. Curva de diversidad de Renyi
# 4. Calculo de diversidad Beta
# 5. Estamación de especies (Chao)
# 6. Curvas de Whittaker

# Authores: Javier Pérez-López *(1,2)*, Ana Wegier *(1)*. 
# 1. Laboratorio de Genética de la Conservación, Jardín Botánico, Instituto de Biología, Universidad Nacional Autónoma de México, Ciudad de México, México.
# 2. Posgrado en Ciencias Biológicas, UNAM.


#####################################################################################
#Variable "Tabla" debe tener el mismo numero de renglones que la variable "factor"
#Variable "factor" debe tener distinto nombre


# Diversidad #

DiversidadCC <- function(Tabla,factor){
  require(vegan)
  require(SpadeR)
  SpecNum <- specnumber(Tabla) # rowSums(BCI > 0) #  Species richness
  ShannonD <- round(diversity(Tabla), 2) # Shannon entropy
  Pielou <- round(ShannonD/log(SpecNum), 2) # Pielou's evenness
  Simp <- round(diversity(Tabla, "simpson"), 2) #Indice de dominacia de Simpson
  Abun<-rowSums(Tabla)
  TablaF <- data.frame(factor,SpecNum, Abun, ShannonD, Simp, Pielou)
  print("Indicadors de Diversidad")
  print(TablaF)
}

# Rarefaccion #
RarefraccionCC <- function(Tabla,factor){
  require(vegan)
  Tabla1 <- data.frame(Tabla, row.names = factor)
  raremax <- min(rowSums(Tabla1))
  col1 <- seq(1:nrow(Tabla1)) #Para poner color a las lineas
  lty1 <- c("solid","dashed","longdash","dotdash")
  rarecurve(Tabla1, sample = raremax, col = "black", lty = lty1, cex = 0.6)
  #Para calcular el numero de especies de acuerdo a rarefraccion
  UUU <- rarefy(Tabla1, raremax)
  print(UUU)
}

# Renyi #
RenyiCC <- function(Tabla, factor){
  require(vegan) #Paquete para la funcion "renyi"
  require(ggplot2)#Paquete para hacer la funcion "qplot"
  require(reshape)#Paquete para la funcion "melt"
  Tabla <- data.frame(Tabla, row.names = factor)
  mod <- renyi(Tabla)
  vec <- seq(1:11)
  mod1 <- data.frame(vec,t(mod))
  mod2 <- melt(mod1, id = c("vec"))
  mod2
  #mod2$variable <- as.numeric(mod2$variable)
  orange <- qplot(vec, value, data = mod2, colour = variable, geom = 	"line") + theme_classic()+ylab("Diversity (alpha)")+xlab("Alfa")+ theme(legend.title = element_blank(), legend.key = element_rect(fill = "white", colour = "white"))
  orange + scale_x_continuous(breaks = c(1,2,3,4,5,6,7,8,9,10,11), 		labels = c("0","0.25","0.5","1","2","4","8","16","32","64","Inf"))
}

# Diversidad Beta#
DivBeta <- function(Tabla, Factor, n1){
  require(vegan) # Paquete para la funcion "betadiver"
  Tabla <- data.frame(Tabla, row.names=Factor)
  Tabla
  Tabla <- decostand(Tabla, "pa")
  mod <- betadiver(Tabla, method=n1)
  mod
  hca <- hclust(mod, method = "ward.D2")
  plot(as.dendrogram(hca), horiz = T, main="Diversidad Beta")
  print(mod)
}

# Estimadcion de Chao # 

ChaoF <- function(Tabla, Factor){
  require(vegan)
  require(ggplot2)
  require(reshape)
  Tabla <- data.frame(Tabla, row.names = Factor)
  Tabla
  UNO <- estimateR(Tabla)[1:2,]
  UNO1 <- melt(UNO)
  names(UNO1)[1] <- c("Var")
  names(UNO1)[2] <- c("Veg")
  LL <- ggplot(UNO1, aes(x = Veg, y = value, fill = Var)) 
  LL <- LL + geom_bar(position = "stack", stat = "identity") 
  LL <- LL + facet_wrap( ~ Var) + theme_bw() + coord_flip()
  print(LL)
  print(UNO)
}

# Curvas de Whittaker #

Whitaker1 <- function(Tabla, Factor){
  require(vegan)
  Tabla <- data.frame(Tabla, row.names = Factor)
  Tabla
  mod <- radfit(Tabla)
  mod
  plot(mod, pch = 19)
}