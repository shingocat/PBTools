#----------------------------------------------------------------
# This function is lifted from the AMMI function of agricolae package which is authoured by Felipe de Mendiburu.
# It draws a polygon or a circumference around the center of the Biplot with a proportional radio at the longest 
# distance of the genotype.
#
# ARGUMENTS:
# model - Object
# distance - Circumference radius >0 and <=1
# shape - Numerical, relating to the shape of the polygon outline
#
# Modifications made by: Nellwyn L. Sales
#----------------------------------------------------------------

AMMI.contour <-
function(model,distance,shape,...)
{
G<- subset(model,model$type=="GEN")
x<-G$PC1
y<-G$PC2
d<-sqrt(x^2+y^2)
r <-distance*max(d)
x<-seq(-r,r,length=shape)
A<-cbind(x,y=sqrt(r^2-x^2))
B<-cbind(x,y=-sqrt(r^2-x^2))
lines(A,type="l",...)
lines(B,type="l",...)
Gin <- d<=r
Gout<- d >r
GEN.in<-row.names(G)[Gin]
GEN.out<-row.names(G)[Gout]
cat("\nLimit, radio:",r)
cat("\nGenotype  in:",length(GEN.in))
cat("\nGenotype out:",length(GEN.out),"\n\n")
distance<-data.frame(row.names=row.names(G),distance=d)
return(list("GENOTYPE IN"=GEN.in, "GENOTYPE OUT"=GEN.out,Distance=distance))
}

