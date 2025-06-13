library("ggplot2")

F100_box <- cbind(FF_100,"n=100")
F200_box <- cbind(FF_200,"n=200")
F400_box <- cbind(FF_400,"n=400")

Fbox <- rbind(F100_box,F200_box,F400_box)
Fbox <- as.data.frame(Fbox)
Fbox[,1] <- as.numeric(as.character(Fbox[,1]))
Fbox[,2] <- as.factor(Fbox[,2])
colnames(Fbox) <- c("PE","label")

pF <- ggplot(Fbox,aes(x=label,y=PE, fill=label))+
  geom_boxplot(outlier.size=1,outlier.alpha=0.2)+ 
  xlab("")+
  ylab(expression(italic("MSE")[italic("f")]))+ scale_color_manual(values=c("tomato4","#56B4E9","#E69F00","darkgreen","gray10"))+ylim(0.0,0.07)+
  theme(legend.position="none",axis.title.y = element_text(size=20),axis.title.x = element_text(size=25),axis.text = element_text(size=20))
plot(pF)


P100_box <- cbind(PP_100,"n=100")
P200_box <- cbind(PP_200,"n=200")
P400_box <- cbind(PP_400,"n=400")

Pbox <- rbind(P100_box,P200_box,P400_box)
Pbox <- as.data.frame(Pbox)
Pbox[,1] <- as.numeric(as.character(Pbox[,1]))
Pbox[,2] <- as.factor(Pbox[,2])
colnames(Pbox) <- c("PE","label")

pP <- ggplot(Pbox,aes(x=label,y=PE, fill=label))+
  geom_boxplot(outlier.size=1,outlier.alpha=0.2)+ 
  xlab("")+
  ylab(expression(italic("MSE")[italic(phi)]))+ scale_color_manual(values=c("tomato4","#56B4E9","#E69F00","darkgreen","gray10"))+ylim(0.013,0.0245)+
  theme(legend.position="none",axis.title.y = element_text(size=20),axis.title.x = element_text(size=25),axis.text = element_text(size=20))
plot(pP)


S100_box <- cbind(SS_100,"n=100")
S200_box <- cbind(SS_200,"n=200")
S400_box <- cbind(SS_400,"n=400")

Sbox <- rbind(S100_box,S200_box,S400_box)
Sbox <- as.data.frame(Sbox)
Sbox[,1] <- as.numeric(as.character(Sbox[,1]))
Sbox[,2] <- as.factor(Sbox[,2])
colnames(Sbox) <- c("PE","label")

pS <- ggplot(Sbox,aes(x=label,y=PE, fill=label))+
  geom_boxplot(outlier.size=1,outlier.alpha=0.2)+ 
  xlab("")+
  ylab(expression(italic("MSE")[italic(zeta)]))+ scale_color_manual(values=c("tomato4","#56B4E9","#E69F00","darkgreen","gray10"))+ylim(0.0,0.044)+
  theme(legend.position="none",axis.title.y = element_text(size=20),axis.title.x = element_text(size=25),axis.text = element_text(size=20))
plot(pS)

