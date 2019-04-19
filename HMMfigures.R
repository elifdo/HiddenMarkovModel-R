install.packages("ggplot2")
install.packages("gridExtra")

library(ggplot2)
require(gridExtra)
library(tidyr)

####### effect of number of states

test_data <-
  data.frame(
    fkm = c(0.8121997,0.809192,0.8106626,0.8106626,0.8121332,0.8048029),
    fkmed = c(0.85575,0.8734848,0.8840477,0.8689414,0.9168449,0.8746089),
    GKfkm= c(0.8181938,0.8136257,0.8326549,0.8418569,0.8626886,0.8731809)  ,
    number_of_states = c(5,6,7,8,9,10)
  )

p1<-  test_data %>%
  gather(state_selection_method,accuracy, fkm, fkmed,GKfkm) %>%
  ggplot(aes(x=number_of_states, y=accuracy, linetype=state_selection_method)) +
  theme(legend.position = c(0.2, 0.85),legend.key.size = unit(3, "mm"))+
  geom_line()+ggtitle("(a) 60 features chosen with hierarchical method")+
  guides(linetype=guide_legend(NULL),shape = guide_legend(override.aes = list(size = 1.5)))+
  ylim(0.79, 0.92)+
  xlab("Number of states") +
  ylab("Accuracy")

test_data2 <-
  data.frame(
    fkm = c(0.8049372,0.814964,0.841968,0.8152307,0.8671888,0.8567404 ),
    fkmed = c(0.8433481,0.8537101,0.8687485,0.8879992,0.8641611,0.8776186),
    GKfkm= c(0.8181938,0.8136257,0.8326549,0.8418569,0.8626886,0.8731809) ,
    number_of_states = c(5,6,7,8,9,10)
  )

p2<-  test_data2 %>%
  gather(state_selection_method,accuracy, fkm, fkmed,GKfkm) %>%
  ggplot(aes(x=number_of_states, y=accuracy, linetype=state_selection_method)) +
  geom_line()+ggtitle("(b) 60 features chosen with k-means method")+
  theme(legend.position = c(0.2, 0.85),legend.key.size = unit(3, "mm"))+
  guides(linetype=guide_legend(NULL),shape = guide_legend(override.aes = list(size = 1.5)))+
  ylim(0.79, 0.92)+
  xlab("Number of states") +
  ylab("Accuracy")

test_data3 <-
  data.frame(
    fkm = c( 0.7973615, 0.8122675,0.8417884,0.8237442, 0.8402539, 0.8760123),
    fkmed = c( 0.8496082,0.8790193,0.8655406,0.8865286,0.894127,0.8790213),
    GKfkm= c(0.8181938,0.8136257,0.8326549,0.8418569,0.8626886,0.8731809) ,
    number_of_states = c(5,6,7,8,9,10)
  )

p3<-  test_data3 %>%
  gather(state_selection_method,accuracy, fkm, fkmed,GKfkm) %>%
  ggplot(aes(x=number_of_states, y=accuracy, linetype=state_selection_method)) +
  guides(linetype=guide_legend(NULL),shape = guide_legend(override.aes = list(size = 1.5)))+
  geom_line( )+ggtitle("(c) 60 features chosen with k-medoids method")+
  theme(legend.position = c(0.2, 0.85),legend.key.size = unit(3, "mm"))+
  ylim(0.79, 0.92)+
  xlab("Number of states") +
  ylab("Accuracy")

test_data4 <-
  data.frame(
    fkm = c(0.8121332,0.8181718,0.818261,0.8152307,0.8508354,0.8361528),
    fkmed = c(0.8807566,0.8566506,0.8718028, 0.9031487,0.8834996,0.8897145),
    GKfkm= c(0.8626886,0.8136257,0.8326549,0.8418569,0.8626886,0.8731809) ,
    number_of_states = c(5,6,7,8,9,10)
  )

p4<-  test_data4 %>%
  gather(state_selection_method,accuracy, fkm, fkmed,GKfkm) %>%
  ggplot(aes(x=number_of_states, y=accuracy, linetype=state_selection_method)) +
  guides(linetype=guide_legend(NULL),shape = guide_legend(override.aes = list(size = 1.5)))+
  geom_line()+ggtitle("(d) Without any feature selection")+
  theme(legend.position = c(0.2, 0.85),legend.key.size = unit(3, "mm"))+
  ylim(0.79, 0.92)+
  xlab("Number of states") +
  ylab("Accuracy")

grid.arrange(p1,p2,p3,p4,ncol=2)


