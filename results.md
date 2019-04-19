You can find referred figures and tables with the name "figures.pdf" in the folder.

For brevity, I will refer to fuzzy k-means method as fkm, Gustafson-
Kessel like fuzzy k-means as GKfkm and fuzzy k-medoids as fkmed. First, we will
compare average accuracy across 10-validation folds formed while model training.

1. Effect of the number of states when other hyper parameters are fixed: Number of
states doesn’t have a linear effect on the accuracy values when other parameters
are fixed. Figure 4 shows accuracy as a function of number of states and feature
selection methods, with number of features used as 60 (a-c) or no features used
(d). As seen in Fig. 4, when GKfkm is used, accuracy increases with the number
of states except from the case when the number of states changes from 5 to 6.
In that case there is a slight decrease. There are no common patterns when other
techniques are used. There are 3 feature selection techniques, 7 different number
of features, 3 state selection techniques, which makes a total of 3*7*3 = 63
possible cases. This number becomes 66 when we include cases when we do not
perform any feature selection. Out of these 66 cases, 40 of them give the best
accuracy when number of states is 10, 14 of them give the best accuracy when
the number of states are 9, 5 cases give when the number of states 8, followed by
number of states 5,6 and 7. Therefore, usually, large number of states give better
results.

2. Effect of the number of features on the accuracy: Figure 5, 6 and 7 show accuracy
as a function of number of features for different feature selection and state
selection methods. As seen in the figures, the effect of the number of features
highly depends on the methods used. Change in the number of features does not
have an effect on accuracy when GKfkm is used. Also, fkm is very robust to
the changes in the number of features only when hierarchical feature selection
method is used. There are no common pattern for the other methods.

3. Effect of feature selection methods on the accuracy: As seen in Fig. 8, when fkm
state selection method is used, k-medoids method gives the best results almost
all the time except a few cases. Hierarchical method gives poor results and it
does not change with the number of features. K-means method, however, doesn’t
follow a common pattern, being worse than k-medoids for most of the times,
but having higher accuracy for a few cases. As seen in Fig. 9, feature selection
methods are more robust, none of them are particularly better than each other
when fkmed method is implemented. Change in the feature selection method
does not have an effect on the accuracy when GKfkm is used. Only exception
is when number of states is 5, accuracy values for this method does not change
for our range of number of features, but changes when we do not implement any
feature selection. When hierarchical feature selection method is used with fkm
state selection, number of features affect the accuracy slightly and accuracy value
is very poor.

4. Effect of the state selection methods on accuracy: When we compare state selection
methods, we see that GKfkm is very robust to the number of features
and the feature selection methods, but accuracy changes when number of states
change. Fkm is also very robust to the changes in the number of features when
hierarchical feature selection method is used. As seen in Fig. 5, when hierarchical
feature selection method is used, fkmed always gives the best results, GKfkm
gives worse results and fkm gives the worst. This is expected since fkmed method
is more robust to outliers in the data and GKfkm method captures non-spherical
patterns unlike fkm. As seen in Fig. 6, a similar pattern appears when k-means
feature selection is implemented, except in some cases fkm surpasses fkmed and
GKfkm. Fkm works more efficiently when it is used with k-medoids feature selection 
method. As seen in Fig. 7, in some cases fkm gives better results than
fkmed and on many cases gives better results than GKfkm. Overall, there are
132 different cases when all other hyper parameters are fixed except state selection
methods. Fkmed method is the best among state selection methods 118
times out of 132 cases, followed by fkm which is the best 14 times and GKfkm
is never the best among other methods.

5. Effect of feature selection on the accuracy: When fkm is implemented for state
selection, only k-medoids and sometimes k-means feature selection is giving
higher accuracy compared to no feature selection. When GK-fkm state selection
method is used, there is no difference between feature selection and no feature
selection. When fkmed is used as state selection, using all the dataset without any
feature selection is giving either best or comparable results to the models with
feature selection.
At the end of training process, we chose hierarchical feature selection method with
60 features and fuzzy k-medoids state selection with 9 number of states, these states
can be seen in Table 1. In their paper, Zhang et al.[35] suggests a method called
multiple property grouping. We applied this method to our dataset and compared
the results on both 10-fold cross validation values on the training data and on test
data. Measures used to compare are smaller on test data than the training data since
the model hasn’t seen test data throughout the training process. As seen in Table
3, the proposed model gives better results than multiple property grouping on both
training and test data on almost all measures except a slightly smaller value of the
precision on test data.