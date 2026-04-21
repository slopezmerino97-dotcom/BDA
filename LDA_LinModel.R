A = read.csv('Clas2006.csv', header = FALSE)
A = as.matrix(A)
meanA = colMeans(A)
Am = sweep(A,2,meanA,FUN = "-") # Centering A
Am
grp = c(array(-1,dim = 10),array(1,dim = 10)) # centered grp
grp

## LDA
reslda <- lda(x = Am,grp)
d = reslda$scaling
d

## Linear model
b = lm(grp ~ Am -1)
b$coefficients

## compare LDA and linear model
d_div_b = d/b$coefficients
df = data.frame(d,b$coefficients,d_div_b)
df
plot(b$coefficients,d)
c = cor(b$coefficients,d)
c
