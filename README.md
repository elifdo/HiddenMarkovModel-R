You can find referred figures and tables with the name "figures.pdf" in the folder.

In this project, I used the HIV-1 protease cleavage 746 dataset[1]. The data contains
lists of octamers (8 amino acids) and a flag depending on whether HIV-1 protease
will cleave in the central position (between amino acids 4 and 5). There are 401
cleaved and 345 non-cleaved octamers. I also use physicochemical properties of
aminoacids from AAIndex database[2]. In this database, there are 544 properties
taken as continuous variables for each aminoacid. I discard 14 of them since they
contain null values.

State Creation: 

In modeling data via HMM, there are 8 bit observation sequences where each
observation is from a set of 20 standard aminoacids, namely, A, R, N, D, C, Q,
E, G, H, I, L, K, M, F, P, S, T, W, Y, V. Each observation has a hidden state behind
it, which is formed by using physicochemical properties of aminoacids. Furthermore,
we accept that if we replace an aminoacid of a cleaved sequence with
another aminoacid having similar properties, it is more likely that the new sequence
will also be a cleaved sequence. Therefore, I grouped aminoacids according to the
similarities based on physicochemical properties and use that information as their
hidden states in the model.
After discarding features with null values, we have 530 features. Usually, when
a clustering algorithm is used with a large number of features, it performs poorly
due to outliers or highly correlated variables. Therefore, I implemented some feature
selection methods to decrease the number of the features. For this purpose,
initially, I grouped the features via the same AAIndex data by treating features as
instances. Here, I use k-means[3], k-medoids[4] and hierarchical clustering[5].
In k-means and hierarchical clustering, we form a subset of features by choosing a
variable randomly from each cluster. In k-medoids, I selected the cluster medoids
as cluster representatives. I tried different number of feature subsets that changes
from 30 to 60 with an increment of 5. I also constructed the model without doing
any feature selection and compare it with the models with feature selection.
By using these subsets of features, I grouped aminoacids to create the states.
An aminoacid can share many different properties with multiple groups. For this
purpose, I performed the fuzzy clustering[6] rather than classical clustering approaches.
By this way, an aminoacid can belong to more than one cluster with a
membership degree between 0 and 1, and sum of the membership degrees adding
to 1. We use fuzzy k-means[6], Gustafson and Kessel-like fuzzy k-means[7] and
fuzzy k-medoids[8]. Fuzzy k-means is similar to usual k-means, where Gustafson
Kessel like fuzzy k-means considers non-spherical clusters too. In fuzzy k-medoids,
medoids are being taken as cluster representatives instead of artificial means. We
assign aminoacids to a state if the membership degree is greater than 0.1, otherwise
not. Since there are 20 aminoacids to cluster, the number of clusters cannot be more
than 10, and we loose too much information and get poor results when it is less
than 5. Therefore, I tried and compared the number of states from 5 to 10. Also, I
standardized the features before clustering to avoid any bias caused by the variance
of the features. In Table 1, 9 states created by fuzzy k-medoids using 60
features chosen by hiererchical clustering is being presented.

Initialization of the EM Algorithm

On the other side, while doing the inference on emission, starting and transition
probabilities, I used the Baum-Welch EM algorithm, which can converge to a local
maximum instead of the global maximum. Therefore, I gave a clever starting
point to the algorithm in order to increase the probability of reaching the global
maximum in our calculation. Hence, I used data in hand to make a good prediction
in the following way:
1. Calculation of initial probabilities: To calculate the starting probability of a state,
I counted all the sequences in the training data which starts with aminoacids
that this state includes. Then, for all states, I divide them to the sum of these
counts to turn these counts into probabilities. For example, let’s say we have 15
sequences in the training data where 5 of them starting with A, 2 starting with
S and 8 starting with V. State 3 contains A, therefore its count will be 2, State
4 contains both A and S, therefore its count will be 2+5, State 6 contains V,
therefore its count will be 8. Probability of first state being State 3 will be 2=17,
first state being State 4 will be 7=17 and first state being State 6 will be 8=17
while all other values in the vector P being 0.
2. Calculation of emission probabilities: To estimate the probability of observing an
aminoacid given a state, I apply the following procedure: with 0.9 probability,
we observe one of the aminoacids that this state includes, and with 0.1 probability,
other aminoacids that this state doesn’t include. As an example, State 1
includes 9 aminoacids where probability 0.9 is equally distributed among them,
each having probability 0.9\9=0.1, and rest of the aminoacids have 0.1 probability
equally distributed among them, each having probability 0.1\11.
3. Calculation of transition probabilities: We know corresponding states for each
aminoacid of our sequences in the data. For example, Table 2 shows the corresponding
states for aminoacids of the sequence AIMALKMR. As seen in Table 2, 
there are 2 transitions from State 5 to State 5, there are 1 transition from State
5 to State 2, etc. given only the sequence AIMALKMR. This way, we count all
transitions coming from all sequences in the training set. Afterwards, for each
state, we sum all transitions from this state to all states including itself and divide
counts of all transitions from this state to this number, this way we turn it into
a probability distribution. In case of the sum being 0, probability of this row is
being equally distributed as 1\N for each state.

Figures 2 and 3 show the count matrix and transition matrix produced by using only
sequence AIMALKMR.

Modeling the data via HMM

I splitted cleaved and non-cleaved data into 90% of training and 10% of test data. I
only use the test data after finding the optimal model parameters through training.
Using initializations for the model parameters, I apply the Baum-Welch EM algorithm
with 1000 maximum number of iterations, and convergence criteria for the
change of log likelihood equal to 0.001. I utilized the aphid R package for this calculation.
Optimum values of the hyper-parameters were selected by using 10-fold cross validation on
the training data. Cross validation is used to reduce the bias stems from the random
selection of data. Accordingly, training data was divided into 10 folds, 9 of them are
used for training and the last one is used for the validation. We repeat this process
10 times until we utilize all 10 folds as the validation data.
To classify the sequence as cleaved or non-cleaved, two separate HMMs are trained
on the cleaved and non-cleaved datasets respectively. We declare these sequences
as cleaved if the likelihood of belonging to cleaved HMM is greater than the noncleaved
HMM and vice versa. This way, we calculate false positive(FP), false negative(
FN), true positive(TP) and true negative(TN) values. To measure the quality of
our classification, we compute precision(pre), recall(rec), accuracy(acc), Matthews
correlation coefficient(MCC) and F-measure(F).
Finally, after deciding the final model, to avoid over optimism caused by the overfitting,
we declare results by using the test data which the final model has not seen
yet.


1. Rognvaldsson, T., You, L., Garwicz, D.: State of the art prediction of HIV-1 protease cleavage
sites. Bioinformatics 31, 1204–1210 (2015)
2. Nakai, K., Kidera, A., Kanehisa, M.: Cluster analysis of amino acid indices for prediction of
protein structure and function. Protein Eng. 2, 93–100 (1988)
3. Hartigan, J. A., Wong, M. A.: Algorithm AS 136: A K-means clustering algorithm. J. Royal
Stat. Soc. Series C (Applied Statistics) 28, 100–108 (1979)
4. Park, H., Jun, C.: A simple and fast algorithm for k-medoids clustering. Expert Syst. Appl.
36, 3336–3341 (2009)
5. Murtagh, F.: Multidimensional Clustering Algorithms. Physica-Verlag (1985)
6. Bezdek J.: Pattern Recognition With Fuzzy Objective Function Algorithms. Plenum Press,
New York (1981)
7. Gustafson, D.E., Kessel, W.C.: Fuzzy clustering with a fuzzy covariance matrix. Proc. IEEE
CDC 761–766 (1978)
8. Krishnapuram, R., Joshi, A., Nasraoui, O., Yi, L.: Low-complexity fuzzy relational clustering
algorithms for Web mining. IEEE Trans. Fuzzy Syst. 9(4), 595–607 (2001)
    