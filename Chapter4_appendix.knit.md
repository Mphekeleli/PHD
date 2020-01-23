---
output: pdf_document
---



\section{Appendix B: R Code for Chapter 4}
\label{sec:Appendix B: R Code for Chapter 4}

\singlespace

Required: R Packages from CRAN (in addition to packages already found in Chapter 3)

\small

```r
if (!require(MuMIn)){
  install.packages("MuMIn")
  library(MuMIn)
}
if (!require(parallel)){
  install.packages("parallel")
  library(parallel)
}
if (!require(R2HTML)){
  install.packages("R2HTML")
  library(R2HTML)
}
if (!require(e1071)){
  install.packages("e1071")
  library(e1071)
}
if (!require(ROCR)){
  install.packages("ROCR")
  library(ROCR)
}
```
\normalsize

\subsection{Data partitioning of the pre-processed OpRisk dataset into Training/Validation/Testing proportions, in preparation for machine learning model building treatments. The original dataset is partitioned into three random subsets initiated by a random number sequence with a randomly selected seed.}
\label{ssec:Data Training/Validation/Testing}

The function \texttt{getmode} specifies the modal class as the reference level in the GLM from which the corresponding observations are estmated and weighted against, see chapter \ref{DATA EXPLORATION AND EXPOSURE VARIABLE ANALYSIS} section \ref{sec:Exploratory data analysis} on page \pageref{sec:Exploratory data analysis} 

\small

```r
building <- TRUE
scoring <- ! building

# Load data
fname <- "file:///C:/Users/Mphekeleli/Documents/R PROJECT/OpRiskPHDGitHub
/OpRisk_PHD_Thesis/Data/OPriskDataSet_exposure.csv"
crs$dataset <- read.csv(fname,
              sep=";",
              dec=",",
              na.strings=c(".", "NA", "", "?"),
              strip.white=TRUE, encoding="UTF-8")
exposure <- crs$dataset[,ncol(crs$dataset)]

# Select variables for loss incident model
crs$dataset <- as.data.frame(crs$dataset)

# The following varaible selections have been noted

crs$input <- crs$dataset %>%
  group_by(UpdatedDay,
           UpdatedTime,
           TradedDay,
           TradedTime,
           Desk,
           CapturedBy,
           TradeStatus,
           TraderId,
           Instrument,
           Reason,
           EventTypeCategoryLevel1,
           BusinessLineLevel1) %>%
  transmute(LossesIndicator = LossIndicator,
            Losses = Loss,
            Exposure = exposure)

getmode <- function(x){
  u <- unique(x)
  as.integer(u[which.max(tabulate(match(x,u)))])
}

for (i in 5:(ncol(crs$input) - 3)){
     crs$input[[i]] <- relevel(crs$input[[i]], getmode(crs$input[[i]]))
}

#==========================================================================================
# A predefined value is used to reset the random seed so that results are repeatable
crv$seed <- 42 

# Build the training/validation/testing datasets. Set parameter values
   

set.seed(crv$seed)    # set random seed to make your partition reproducible

crs$nobs <- nrow(crs$input)                 # nobs=2331

crs$train <- sample(crs$nobs, 0.7*crs$nobs) # proportion of training data = 1632 

crs$nobs %>%
  seq_len() %>%
  setdiff(crs$train) %>%
  sample(0.15*crs$nobs) ->                  # proportion of validation data = 350 
  crs$validate

crs$nobs %>%
  seq_len() %>%
  setdiff(crs$train) %>%
  setdiff(crs$validate) ->                  # proportion of testing data = 349 
  crs$test


crs$training <- as.data.frame(crs$input[crs$train,])
crs$validation <- as.data.frame(crs$input[crs$validate,])
crs$testing <- as.data.frame(crs$input[crs$test,])
```
\normalsize

\doublespacing

\section{Models}
\label{sec:Models}

\singlespace

\subsection{Estimation of some poisson regression models for OpRisk loss frequency distribution: To build the model we pass on to the model building function \texttt{glm} i.e., the formula that describes the model to build.}
\label{ssec:Estimation of some poisson regression models for OpRisk loss frequency distribution}

use "LossesIndicator" as the dependent variable, \texttt{TradedDay} and \texttt{Desk} variables are predictor variables. 

\small

\normalsize

<!-- \begin{figure} -->
<!-- \centering -->
<!-- \includegraphics[scale=1.0]{SUMfreqfit1.pdf} -->
<!-- \caption[Poisson GLM Summary statistics]{Estimation of some Poisson distribution for target variable \texttt{LossesIndicator} and the two explanatory variables \texttt{TradedDay} and \texttt{Desk}} -->
<!-- \label{Two variable poisson estimation} -->
<!-- \end{figure} -->

\emph{Basic model build summary}

\begin{verbatim}
Call:
glm(formula = LossesIndicator ~ TradedDay + Desk, family = poisson(link = "log"), 
    data = crs$training, offset = log(Exposure))

Deviance Residuals: 
    Min       1Q   Median       3Q      Max  
-2.8706  -0.5300  -0.2286  -0.0545   4.3750  

Coefficients:
                      Estimate Std. Error z value             Pr(>|z|)    
(Intercept)          -8.053221   0.193632 -41.590 < 0.0000000000000002 ***
TradedDay            -0.014087   0.006381  -2.208             0.027266 *  
DeskAfrica            1.457695   0.413087   3.529             0.000417 ***
DeskBonds/Repos       1.764230   0.254571   6.930      0.0000000000042 ***
DeskCommodities       0.924033   0.235749   3.920      0.0000887114575 ***
DeskDerivatives      -0.577626   0.344672  -1.676             0.093763 .  
DeskEquity            1.365152   0.225948   6.042      0.0000000015232 ***
DeskManagement/Other -1.410706   1.014052  -1.391             0.164177    
DeskMM                0.339561   0.241501   1.406             0.159711    
DeskPrime Services    2.129594   0.223030   9.548 < 0.0000000000000002 ***
DeskSND              -0.716361   0.283796  -2.524             0.011596 *  
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

(Dispersion parameter for poisson family taken to be 1)

    Null deviance: 1994.9  on 1630  degrees of freedom
Residual deviance: 1764.7  on 1620  degrees of freedom
AIC: 2310.7

Number of Fisher Scoring iterations: 8
\end{verbatim}

\subsection{Let us now fit a broader GLM to be our global model, still \texttt{LossesIndicator} as the dependent variable, while the rest of the predictive variables will be predictor variables.}
\label{ssec:GLM estimation results}

\small

```r
freqfit <- glm(formula = LossesIndicator ~ UpdatedDay + UpdatedTime + TradedDay + 
    TradedTime + Desk + CapturedBy + TradeStatus + TraderId + 
    Instrument + Reason + EventTypeCategoryLevel1 + BusinessLineLevel1, 
    family = poisson(link = "log"), data = crs$training, 
    offset = log(Exposure)))
```
\normalsize

\emph{Global model build summary}

\begin{verbatim}

Call:
glm(formula = LossesIndicator ~ UpdatedDay + UpdatedTime + TradedDay + 
    TradedTime + Desk + CapturedBy + TradeStatus + TraderId + 
    Instrument + Reason + EventTypeCategoryLevel1 + BusinessLineLevel1, 
    family = poisson(link = "log"), data = crs$training, 
    offset = log(Exposure))

Deviance Residuals: 
    Min       1Q   Median       3Q      Max  
-4.7605  -0.3575  -0.1105  -0.0315   3.9754  

Deviance Residuals: 
    Min       1Q   Median       3Q      Max  
-4.7605  -0.3575  -0.1105  -0.0315   3.9754  

Coefficients:
                       Estimate     Std. Error z value             Pr(>|z|)    
(Intercept)           -7.8113982    0.6701090 -11.657 < 0.0000000000000002 ***
UpdatedDay             0.0043298    0.0097374   0.445             0.656565    
UpdatedTime            0.5850748    0.7146037   0.819             0.412935    
TradedDay             -0.0001216    0.0082899  -0.015             0.988297    
TradedTime            -0.8201083    0.7856849  -1.044             0.296572    
DeskAfrica             1.5831846    0.5281939   2.997             0.002723 ** 
DeskBonds/Repos        2.3432611    0.3958404   5.920     0.00000000322505 ***
DeskCommodities        0.6910019    0.4627163   1.493             0.135343    
DeskDerivatives        0.3832286    0.5020354   0.763             0.445255    
DeskEquity             1.1778919    0.3969116   2.968             0.003001 ** 
DeskManagement/Other -14.7535497  367.5957685  -0.040             0.967985    
DeskMM                 1.2027479    0.5903580   2.037             0.041618 *  
DeskPrime Services    -0.1311216    1.3579826  -0.097             0.923079    
DeskSND                0.0518179    0.7105464   0.073             0.941864    
CapturedByMIDOFFICE    0.2618294    0.3319433   0.789             0.430242    
CapturedByPROD ACCOUN
TANT                   0.0934629    0.5393415   0.173             0.862423    
CapturedByPROD 
CONTROLLER            -0.4702595    0.3321428  -1.416             0.156824    
CapturedByUNAUTHORISED-0.9738648    0.5370115  -1.813             0.069756 .  
TradeStatusBO-BO 
Confirmed             -0.3379764    0.2363953  -1.430             0.152801    
TradeStatusTerminated  2.1847205    1.3465905   1.622             0.104716    
TradeStatusTerminated
/Void                 11.3395140 1443.2959040  -0.008             0.993731    
TraderIdAMBA           0.4943858    0.5830511   0.848             0.396478    
TraderIdANALYST       -0.2062667    0.2944408  -0.701             0.483592    
TraderIdASSOCIATE     -0.7844915    0.3507247  -2.237             0.025301 *  
TraderIdATS            2.7382867    0.6662928   4.110     0.00003961139187 ***
TraderIdMNGDIRECTOR    0.6153055    0.3144102   1.957             0.050346 .  
TraderIdVICE PRINCIPAL 0.2229076    0.3662304   0.609             0.542754    
InstrumentBill         0.7859508    0.6781491   1.159             0.246471    
InstrumentBond        -0.3433677    0.4080143  -0.842             0.400035    
InstrumentBuySellback  0.1566566    0.6720075   0.233             0.815670    
InstrumentCall Deposit-1.2603996    0.3973057  -3.172             0.001512 ** 
InstrumentCombination  0.6328848    1.2513624   0.506             0.613028    
InstrumentCredit
DefaultSwap           -1.9878550    0.6529647  -3.044             0.002332 ** 
InstrumentCurr        -1.0355350    0.5045561  -2.052             0.040134 *  
InstrumentCurrSwap    -0.6964712    0.5165406  -1.348             0.177550    
InstrumentDeposit      0.0584655    0.3945317   0.148             0.882193    
InstrumentEquityIndex -0.8243925    1.0794767  -0.764             0.445048    
InstrumentETF          1.0751791    1.1397155   0.943             0.345489    
InstrumentFRA         -1.6745203    0.9177750  -1.825             0.068070 .  
InstrumentFRN          0.2792480    0.5903548   0.473             0.636201    
InstrumentFuture/
Forward               -0.7303512    0.3917942  -1.864             0.062305 .  
InstrumentIndexLinked
Bond                  -0.6022619    0.7595690  -0.793             0.427836    
InstrumentIndexLinked
Swap                  -0.0154039    0.6024804  -0.026             0.979602    
InstrumentOption       0.0277106    0.5621959   0.049             0.960688    
InstrumentOther       -0.1830206    0.2891929  -0.633             0.526821    
InstrumentRepo/Reverse-2.0252367    0.7133321  -2.839             0.004524 ** 
InstrumentSecurityLoan-1.0447821    0.3461497  -3.018             0.002542 ** 
InstrumentStock       -0.6782734    1.4169911  -0.479             0.632172    
InstrumentTotalReturn
Swap                  -1.0219974    0.3990992  -2.561             0.010444 *  
ReasonAcquirer         4.2877314    1.2671128   3.384             0.000715 ***
ReasonBrokerage 
Related                1.7788236    1.0318752   1.724             0.084730 .  
ReasonCalendar Related 1.2768935    0.3888510   3.284             0.001024 ** 
ReasonCapture Errors    
Direct/Amount/Rate/CP  0.2418004    0.4150033   0.583             0.560131  
ReasonClient Request
to Amend Economics of
Deal		      -0.9290734    1.0584940  -0.878             0.380089   
ReasonCommodities
Early Delivery        -0.4631230    1.2048262  -0.384             0.700689    
ReasonCorrected 
NDEUSSA Reset Level   -2.8475291    0.8317484  -3.424             0.000618 ***
ReasonCPI Fixings     -4.2562507    0.9924649  -4.289     0.00001798306834 ***
ReasonEarly Terminati
on/Delivery/ Close out-2.6225044    0.7518160  -3.488             0.000486 ***
ReasonFees/Commissions
Related               -1.4712607    0.9402884  -1.565             0.117655    
ReasonOperations 
request to change 
Economics of Trade    -0.1702668    0.8638805  -0.197             0.843753    
ReasonPayments Related-0.5039834    0.7179092  -0.702             0.482669    
ReasonPortfolio Move
 / Restructure        -1.8408735    1.1353256  -1.621             0.104921    
ReasonSales Credits    2.1487308    1.1446783   1.877             0.060498 .  
ReasonSystem Update 
Call Accounts         -4.2209045    0.9534873  -4.427     0.00000956382022 ***
ReasonTrade 
Restructure           -0.6497247    2.2083716  -0.294             0.768598    
ReasonTri-Optima      -3.2769186    1.4560359  -2.251             0.024412 *  
ReasonValuation Group 
Request                0.7265138    0.6359340   1.142             0.253273    
EventTypeCategoryLevel
1EL1                   2.4984605    0.6443312   3.878             0.000105 ***
EventTypeCategoryLevel
1EL4                   1.3739551    0.6054219   2.269             0.023243 *  
EventTypeCategoryLevel
1EL6                   5.6915077    0.8154721   6.979     0.00000000000296 ***
EventTypeCategoryLevel
1EL8                  -4.0471414    1.1004892  -3.678             0.000235 ***
BusinessLineLevel
1BL1                   2.0456821    0.7730115   2.646             0.008136 ** 
BusinessLineLevel
1BL3                   0.1576673    0.6308477   0.250             0.802642    
BusinessLineLevel
1BL4                  -1.2651162    0.4937915  -2.562             0.010406 *  
BusinessLineLevel
1BL5                  -1.1905846    0.5243289  -2.271             0.023166 *  
BusinessLineLevel
1BL6                   1.7544025    1.4007670   1.252             0.210403    
BusinessLineLevel
1BL7                   2.2553305    2.2524125   1.001             0.316684    
BusinessLineLevel
1BL9                   0.2496745 1613.4575790   0.000             0.999877    
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

(Dispersion parameter for poisson family taken to be 1)

    Null deviance: 1898.7  on 1630  degrees of freedom
Residual deviance: 1228.1  on 1553  degrees of freedom
AIC: 1888.1

Number of Fisher Scoring iterations: 15
\end{verbatim}

\subsection{Estimation of some binomial regression modelS: The formula that describes the model to build.}

The target variable is "remapped" or transformed into a factor variable from a numerical variable by use of the code below:

\small

```r
crs$dataset[["TFC_LossIndicator"]] <- 
  as.factor(crs$dataset[["LossIndicator"]])

ol <- levels(crs$dataset[["TFC_LossIndicator"]])
lol <- length(ol)
nl <- c(sprintf("[%s,%s]", ol[1], ol[1]), 
        sprintf("(%s,%s]", ol[-lol], ol[-1]))
levels(crs$dataset[["TFC_LossIndicator"]]) <- nl
```

\emph{We will use "LossesIndicator" as the dependent variable.} 

\small

```r
 freqfit <- glm(LossesIndicator ~ UpdatedDay + UpdatedTime + TradedDay
        + TradedTime + Desk + CapturedBy + TradeStatus + TraderId +
   Instrument + Reason + EventTypeCategoryLevel1 + BusinessLineLevel1,
   data=crs$training, family=binomial(link="logit"), offset=log(Exposure))

 summary(freqfit)
```
\normalsize

\emph{Basic model build summary}

\begin{verbatim}
Call:
glm(formula = LossesIndicator ~ UpdatedDay + UpdatedTime + TradedDay + 
    TradedTime + Desk + CapturedBy + TradeStatus + TraderId + 
    Instrument + Reason + EventTypeCategoryLevel1 + BusinessLineLevel1, 
    family = binomial(link = "logit"), data = crs$training, offset = log(Exposure))

Deviance Residuals: 
    Min       1Q   Median       3Q      Max  
-3.2957  -0.4346  -0.0987  -0.0001   4.0153  

Coefficients:
                       Estimate      Std. Error z value           Pr(>|z|)    
(Intercept)           -7.9749408     1.0332350  -7.718 0.0000000000000118 ***
UpdatedDay             0.0158775     0.0144454   1.099           0.271707    
UpdatedTime           -0.5014732     0.9397316  -0.534           0.593594    
TradedDay             -0.0007336     0.0120301  -0.061           0.951373    
TradedTime             0.9590326     1.0673571   0.899           0.368913    
DeskAfrica             4.3345648     0.9075486   4.776 0.0000017870595853 ***
DeskBonds/Repos        5.0241306     0.6499654   7.730 0.0000000000000108 ***
DeskCommodities        2.0972543     0.6780211   3.093           0.001980 ** 
DeskDerivatives       -0.1864882     0.7413218  -0.252           0.801380    
DeskEquity             1.8915605     0.5960148   3.174           0.001505 ** 
DeskManagement/Other -18.2048936  1520.1898190  -0.012           0.990445    
DeskMM                 1.9981588     0.7940263   2.516           0.011853 *  
DeskPrime Services     0.5307160     1.3075756   0.406           0.684832    
DeskSND                2.3274858     1.0644292   2.187           0.028771 *  
CapturedByMIDOFFICE    0.8804690     0.4547261   1.936           0.052836 .  
CapturedByPROD ACCOUN
TANT                   0.4517261     0.6963223   0.649           0.516512    
CapturedByPROD 
CONTROLLER            -1.2989586     0.4743150  -2.739           0.006170 ** 
CapturedByUNAUTHORISED 0.3986375     0.8119582   0.491           0.623456    
TradeStatusBO-BO 
Confirmed             -1.5378456     0.3303144  -4.656 0.0000032287903266 ***
TradeStatusTerminated  1.5155682     1.3644806   1.111           0.266685    
TradeStatusTerminated
/Void                -15.5553055  4600.1942577  -0.003           0.997302    
TraderIdAMBA          -0.0025153     0.7912919  -0.003           0.997464    
TraderIdANALYST       -0.6145646     0.4385093  -1.401           0.161069    
TraderIdASSOCIATE     -0.9420685     0.5319834  -1.771           0.076584 .  
TraderIdATS            2.6800452     0.9148818   2.929           0.003396 ** 
TraderIdMNGDIRECTOR    1.3616232     0.5228248   2.604           0.009205 ** 
TraderIdVICE PRINCIPAL-0.7050650     0.5070382  -1.391           0.164360    
InstrumentBill        -0.8718713     0.9354385  -0.932           0.351313    
InstrumentBond         0.2366580     0.6067829   0.390           0.696521    
InstrumentBuySellback  1.3218806     1.0448559   1.265           0.205824    
InstrumentCall Deposit-2.1675836     0.5991862  -3.618           0.000297 ***
InstrumentCombination-16.7614537   746.7004950  -0.022           0.982091    
InstrumentCredit
DefaultSwap           -1.5055382     0.9616767  -1.566           0.117458    
InstrumentCurr        -2.8712232     0.7263157  -3.953 0.0000771343104836 ***
InstrumentCurrSwap    -2.0228703     0.7723046  -2.619           0.008812 ** 
InstrumentDeposit      0.2775796     0.5010409   0.554           0.579575    
InstrumentEquityIndex-20.9093522  4247.6245656  -0.005           0.996072    
InstrumentETF          0.6949521     2.5050338   0.277           0.781456    
InstrumentFRA        -18.4966479  2217.1713497  -0.008           0.993344    
InstrumentFRN          1.0640690     0.6934803   1.534           0.124934    
InstrumentFuture
/Forward              -1.3009786     0.5043702  -2.579           0.009897 ** 
InstrumentIndexLinked
Bond                 -18.5083821  1980.8539595  -0.009           0.992545    
InstrumentIndexLinked
Swap                  -1.3276177     0.8399563  -1.581           0.113974    
InstrumentOption      -0.3366249     0.7407203  -0.454           0.649501    
InstrumentOther       -0.7233858     0.4259873  -1.698           0.089481 .  
InstrumentRepo/Reverse-1.6068842     0.8756603  -1.835           0.066498 .  
InstrumentSecurityLoan-1.6209368     0.5052288  -3.208           0.001335 ** 
InstrumentStock      -15.5171850   917.6844311  -0.017           0.986509    
InstrumentTotalReturn
Swap                  -2.1376528     0.5110535  -4.183 0.0000287895194553 ***
ReasonAcquirer       -15.7038902 10686.9945316  -0.001           0.998828    
ReasonBrokerage Related4.1338623     2.1977014   1.881           0.059973 .  
ReasonCalendar Related 2.8178751     0.5313717   5.303 0.0000001139021850 ***
ReasonCapture Errors
 - Direction / Amount
 / Rate / CP           1.4848841     0.5279634   2.812           0.004916 ** 
ReasonClient Request to
 Amend Economics 
of Deal               -0.3720879     1.2198215  -0.305           0.760340    
ReasonCommodities
Early Delivery         2.4241223     1.2227799   1.982           0.047427 *  
ReasonCorrected 
NDEUSSA Reset Level    0.4751118     0.9464671   0.502           0.615678    
ReasonCPI Fixings     -16.6208128  1172.4035169  -0.014           0.988689    
ReasonEarly 
Termination / Delivery
 / Close out          -43.1745557  1552.0277206  -0.028           0.977807    
ReasonFees/Commissions
Related                -0.8431985     1.0953186  -0.770           0.441407    
ReasonOperations request
 to change Economics
 of Trade             -40.9898225  1552.0278612  -0.026           0.978930    
ReasonPayments Related  0.6979688     0.7737247   0.902           0.367009    
ReasonPortfolio Move
 / Restructure        -14.7140963   850.1374455  -0.017           0.986191    
ReasonSales Credits     2.9232638     1.3710077   2.132           0.032990 *  
ReasonSystem Update Call
 Accounts             -79.5591356  1906.0449509  -0.042           0.966706    
ReasonTrade 
Restructure           -55.6249278  2864.4577747  -0.019           0.984507    
ReasonTri-Optima       -1.9567842     1.4525216  -1.347           0.177928    
ReasonValuation Group
 Request                2.2511805     0.9782985   2.301           0.021385 *  
EventTypeCategoryLevel1
EL1                    44.5975989  1552.0274644   0.029           0.977076    
EventTypeCategoryLevel1
EL4                     1.2665117     0.6101433   2.076           0.037916 *  
EventTypeCategoryLevel1
EL6                    81.6802380  1906.0447914   0.043           0.965819    
EventTypeCategoryLevel1
EL8                    -5.4027759     1.3335661  -4.051 0.0000509176014626 ***
BusinessLineLevel1
BL1                     2.7388864     1.0366354   2.642           0.008240 ** 
BusinessLineLevel1
BL3                     -1.6356998     0.8999922  -1.817           0.069147 .  
BusinessLineLevel1
BL4                     -4.1488349     0.7241412  -5.729 0.0000000100835606 ***
BusinessLineLevel1
BL5                     -2.2516395     0.7017310  -3.209           0.001333 ** 
BusinessLineLevel1
BL6                      0.8510447     1.4106248   0.603           0.546302    
BusinessLineLevel1
BL7                    -16.5197686  2691.6096464  -0.006           0.995103    
BusinessLineLevel1
BL9                     24.2348398  1520.1905636   0.016           0.987281    
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

(Dispersion parameter for binomial family taken to be 1)

    Null deviance: 2380.8  on 1630  degrees of freedom
Residual deviance: 1270.6  on 1553  degrees of freedom
AIC: 1426.6

Number of Fisher Scoring iterations: 18
\end{verbatim}

\section{Model selection and multimodel inference: MuMIn}
\label{sec:Model selection and multimodel inference: MuMIn}

\singlespace

\subsection{"All possible models" are considered by subjectively ad iteratively searching the data for patterns and "significance".}
\label{ssec:Data_mining}

we use "dredge" function to generate models using combinations of the terms in the global model. The function will also calculate AICc values and rank models according to it.

\small

```r
 cl <- makeCluster(2) # Assign R cores to the job
 options(na.action=na.fail)
 freqfits <- dredge(freqfit)
 stopCluster(cl)
 freqfits
```
\normalsize

"MuMLn::dredge" returns a list of $4097$ models, below is the build summary

\begin{figure}
\centering
\includegraphics[height=20cm, width=15cm]{Dredge_bin.pdf}
\caption[Data dredging]{Model selection data mining exercise}
\label{Dredge}
\end{figure}

\subsection{We use \texttt{get.models} function to generate a list in which its objects are the fitted models. We will also use the "model.avg" function to do a model averaging based on AICc.}
\label{sec:Model averaging function}

\small

```r
cl <- makeCluster(2) # Assign R cores to the job 
Amodel <- model.avg(get.models(freqfits, subset = TRUE))
summary(Amodel)
stopCluster(cl) 
```
\normalsize

Now we have AICc values for our models and we have the average model (or mean model).

\begin{figure}
\centering
\includegraphics[scale=1.0]{Average_model_Training_1.pdf}
\caption[]{Estimation of some Poisson distribution for target variable \texttt{LossesIndicator} showing the 4096 component models and their respective AICc values for our models giving rise to the average (mean) model \texttt{Amodel}}
\label{AModel_Summary1}
\end{figure}

<!-- \begin{figure} -->
<!-- \centering -->
<!-- \includegraphics[height=20cm, width=15cm]{Get_models_bin1.eps} -->
<!-- \caption[Model averaging]{component models for computing the average (mean) model \texttt{AModel}} -->
<!-- \label{AModel_Summary1} -->
<!-- \end{figure} -->

\begin{figure}
\centering
\includegraphics[height=20cm, width=15cm]{Get_models_bin2.pdf}
\caption[Poisson GLM Amodel summary statistics]{Continued from \ref{AModel_Summary1}}
\label{AModel_Summary2}
\end{figure}

\begin{figure}
\centering
\includegraphics[height=20cm, width=15cm]{Get_models_bin3.pdf}
\caption[Model averaging]{Continued from \ref{AModel_Summary2}}
\label{Poisson GLM Amodel summary statistics continued}
\end{figure}

\subsection{Evaluate model performance on the test dataset}
\label{sec:Evaluate model performance on the test dataset}

Obtain the response from the Linear model.

\small

```r
 Av.PredTT <- predict(Amodel, crs$testing, type = "response")
 
 # Export into excel
 
HTMLStart(); HTML(data.frame(Av.PredTT)); w <- HTMLStop()
browseURL(w)
 
shell(paste("start excel", w))
Est <- "file:///C:/Users/Mphekeleli/Documents/R PROJECT/OpRiskPHDGitHub/
OpRisk_PHD_Thesis/Data/OPriskDataSet_GoF_Amodel_InCode.csv"
pred <- read.csv(Est,
                   sep=",",
                   dec=".",
                   na.strings=c(".", "NA", "", "?"),
                   strip.white=TRUE, encoding="UTF-8")
```
\normalsize

Generate the confusion matrix showing counts and generate an ROC curve for the GLM model on the appropriate dataset partition


```r
confusionMatrix(table(pred$response, crs[["testing"]][["LossesIndicator"]]))

crs$pr <- Av.PredTT

# Remove observations with missing target.

no.miss   <- na.omit(crs[["testing"]][["LossesIndicator"]])
miss.list <- attr(no.miss, "na.action")
attributes(no.miss) <- NULL

if (length(miss.list))
{
  predic <- prediction(crs$pr[-miss.list], no.miss)
} else
{
  predic <- prediction(crs$pr, no.miss)
}

pe <- performance(predic, "tpr", "fpr")
au <- performance(predic, "auc")@y.values[[1]]
pd <- data.frame(fpr=unlist(pe@x.values), tpr=unlist(pe@y.values))
p <- ggplot(pd, aes(x=fpr, y=tpr))
p <- p + geom_line(colour="red")
p <- p + xlab("False Positive Rate") + ylab("True Positive Rate")
p <- p + ggtitle("ROC Curve Linear OPriskDataSet_exposure.csv [validate] LossesIndicator")
p <- p + theme(plot.title=element_text(size=10))
p <- p + geom_line(data=data.frame(), aes(x=c(0,1), y=c(0,1)), colour="grey")
p <- p + annotate("text", x=0.50, y=0.00, hjust=0, vjust=0, size=5,
                  label=paste("AUC =", round(au, 2)))
print(p)

# Calculate the area under the curve for the plot.


# Remove observations with missing target.

no.miss   <- na.omit(crs[["testing"]][["LossesIndicator"]])
miss.list <- attr(no.miss, "na.action")
attributes(no.miss) <- NULL

if (length(miss.list))
{
  predic <- prediction(crs$pr[-miss.list], no.miss)
} else
{
  predic <- prediction(crs$pr, no.miss)
}
performance(predic, "auc")
```

\begin{figure}
\centering
\includegraphics[height=20cm, width=15cm]{ConfusionMatrix_Appendix.eps}
\caption[Confusion Matrices comparison]{Comparison of Poisson, Binomial GLM's confusion matrices for Training, Validation \& Testing datasets}
\label{ConfusionMatrixAll}
\end{figure}

\begin{figure}[t!] % "[t!]" placement specifier just for this example
\begin{subfigure}{0.48\textwidth}
\includegraphics[width=\linewidth]{POI_ROC_Training.pdf}
\caption{AUC for poisson parameter depicting model's performance on the training set} \label{POI_ROC_Train}
\end{subfigure}\hspace*{\fill}
\begin{subfigure}{0.48\textwidth}
\includegraphics[width=\linewidth]{BIN_ROC_Training.pdf}
\caption{AUC for binomial parameters depicting model's performance on the training set} \label{BIN_ROC_Train}
\end{subfigure}

\medskip
\begin{subfigure}{0.48\textwidth}
\includegraphics[width=\linewidth]{POI_ROC_Validation.pdf}
\caption{As in (a) but depicting model's performance on the validation set} \label{POI_ROC_Validate}
\end{subfigure}\hspace*{\fill}
\begin{subfigure}{0.48\textwidth}
\includegraphics[width=\linewidth]{BIN_ROC_Validation.pdf}
\caption{As in (b) but depicting model's performance on the validation set} \label{BIN_ROC_Validate}
\end{subfigure}

\medskip
\begin{subfigure}{0.48\textwidth}
\includegraphics[width=\linewidth]{POI_ROC_Testing.pdf}
\caption{As in (a) \& (c) depicting model's performance on the testing set} \label{POI_ROC_Test}
\end{subfigure}\hspace*{\fill}
\begin{subfigure}{0.48\textwidth}
\includegraphics[width=\linewidth]{BIN_ROC_Testing.pdf}
\caption{As in (b) \& (d) depicting model's performance on the testing set} \label{BIN_ROC_Test}
\end{subfigure}

\caption{Comparison of area under the curve (AUC) for the Poisson, Binomial GLM ROC curve for training, validation and testing datasets} \label{ROCCurveAll}
\end{figure}

\clearpage

\section{Data augmentation code: Extrapolation simulation in Matlab}
\label{sec:Data augmentation code: Extrapolation simulation in Matlab}

Notably, the code for the predict condition was run via the Matlab Terminal:

\small

```bash

% Updated Time
% generate the vector DD

DDD = 1:31;

% generate the vector VVV

VVV = 1:12;

% generate the vector UUU

UUU = 2013 : -1 : 2006;

% making the full thrity five days vector
% Years
Thirty_five_days = [UUU';UUU';UUU';UUU(1:end-1)'];
%Months
Thirty_five_days2 = [VVV';VVV';VVV(1:7)'];
% Days
Thirty_five_days3 = DDD';

% The updated time algorithm

% initializing the time matrix
for i = 1 : length(UUU)
    
    UUU_trans{i} = num2cell(zeros(1,12));
    
end

for i = 1 : length(UUU)
    for j = 1 : length(UUU_trans{1,1})
        
        UUU_TRANS{1,i}{1,j} = num2cell(zeros(31,4));
    end
    
end

% The number of random numbers
H = 1000;
UPD =.789155092592539;
format long
% filling in the time matrix
for i = 1 : length(UUU)
    for j = 1 : length(UUU_trans{1,1})
        
        UUU_TRANS{1,i}{1,j}(:,end) = num2cell((1:31)');
        UUU_TRANS{1,i}{1,j}(:,end-1) = num2cell(VVV(j));
        
        
        
    end
    
end

UUU = num2cell(UUU);
%         UUU = sortrows(UUU,2);


for i = 1 : length(UUU)
    for j = 1 : length(UUU_trans{1,1})
        for k = 1 : length(UUU_TRANS{1,5}{1,1})
            % PART 1
            UUU{i,j} =num2cell((((i)^(0)).*((j)^(0)).*rand(1,31)));
            UUU{i,j} = UUU{i,j}';
            UUU{i,j}(:,2) = num2cell(Thirty_five_days(:,1));
            UUU{i,j} = sortrows(UUU{i,j},1);
            
            %PART 2
            UUU2{i,j} =num2cell((((i)^(0)).*((j)^(0)).*rand(1,31)));
            UUU2{i,j} = UUU2{i,j}';
            UUU2{i,j}(:,2) = num2cell(Thirty_five_days2(:,1));
            UUU2{i,j} = sortrows(UUU2{i,j},1);
            % PART 3
            UUU3{i,j} =num2cell((((i)^(0)).*((j)^(0)).*rand(1,31)));
            UUU3{i,j} = UUU3{i,j}';
            UUU3{i,j}(:,2) = num2cell(Thirty_five_days3(:,1));
            UUU3{i,j} = sortrows(UUU3{i,j},1);
            % PART 1
            UUU_TRANS{1,i}{1,j}(k,end-3) = UUU{i,j}(k,2);
            
            % PART 2
            UUU_TRANS{1,i}{1,j}(k,end-2) = UUU2{i,j}(k,2);
            %PART3
            UUU_TRANS{1,i}{1,j}(k,end-1) = UUU3{i,j}(k,2);
            rH{i,j} = num2cell(((i)^(0).*(j)^(0)).*rand(H,1));
            yH{i,j} = rH{i,j}(cell2mat( rH{i,j}) <= UPD);
            gH{i,j} = num2cell(cell2mat( yH{i,j}(1:31)));
            UUU_TRANS{1,i}{1,j}(k,end) = gH{i,j}(k,1);
        end
    end
    
end

UPDATED_TIME = UUU_TRANS;


%% Traded time

% generate the vector DD

DDDT = 1:31;

% generate the vector VVV

VVVT = 1:12;

% generate the vector UUU

UUUT = 2013 : -1 : 2006;

% making the full thrity five days vector
% Years
Thirty_five_daysT = [UUUT';UUUT';UUUT';UUUT(1:end-1)'];
%Months
Thirty_five_days2T = [VVVT';VVVT';VVVT(1:7)'];
% Days
Thirty_five_days3T = DDDT';

% The updated time algorithm

% initializing the time matrix
for i = 1 : length(UUUT)
    
    UUU_transT{i} = num2cell(zeros(1,12));
    
end

for i = 1 : length(UUUT)
    for j = 1 : length(UUU_transT{1,1})
        
        UUU_TRANST{1,i}{1,j} = num2cell(zeros(31,4));
    end
    
end

% The number of random numbers
HT = 1000;
UPDT =.789155092592539;
format long
% filling in the time matrix
for i = 1 : length(UUUT)
    for j = 1 : length(UUU_transT{1,1})
        
        UUU_TRANST{1,i}{1,j}(:,end) = num2cell((1:31)');
        UUU_TRANST{1,i}{1,j}(:,end-1) = num2cell(VVVT(j));
        
        
        
    end
    
end

UUUT = num2cell(UUUT);
%         UUU = sortrows(UUU,2);


for i = 1 : length(UUUT)
    for j = 1 : length(UUU_transT{1,1})
        for k = 1 : length(UUU_TRANST{1,5}{1,1})
            % PART 1
            UUUT{i,j} =num2cell((((i)^(0)).*((j)^(0)).*rand(1,31)));
            UUUT{i,j} = UUUT{i,j}';
            UUUT{i,j}(:,2) = num2cell(Thirty_five_daysT(:,1));
            UUUT{i,j} = sortrows(UUUT{i,j},1);
            
            %PART 2
            UUU2T{i,j} =num2cell((((i)^(0)).*((j)^(0)).*rand(1,31)));
            UUU2T{i,j} = UUU2T{i,j}';
            UUU2T{i,j}(:,2) = num2cell(Thirty_five_days2T(:,1));
            UUU2T{i,j} = sortrows(UUU2T{i,j},1);
            % PART 3
            UUU3T{i,j} =num2cell((((i)^(0)).*((j)^(0)).*rand(1,31)));
            UUU3T{i,j} = UUU3T{i,j}';
            UUU3T{i,j}(:,2) = num2cell(Thirty_five_days3T(:,1));
            UUU3T{i,j} = sortrows(UUU3T{i,j},1);
            % PART 1
            UUU_TRANST{1,i}{1,j}(k,end-3) = UUUT{i,j}(k,2);
            
            % PART 2
            UUU_TRANST{1,i}{1,j}(k,end-2) = UUU2T{i,j}(k,2);
            %PART3
            UUU_TRANST{1,i}{1,j}(k,end-1) = UUU3T{i,j}(k,2);
            rHT{i,j} = num2cell(((i)^(0).*(j)^(0)).*rand(HT,1));
            yHT{i,j} = rHT{i,j}(cell2mat( rHT{i,j}) <= UPDT);
            gHT{i,j} = num2cell(cell2mat( yHT{i,j}(1:31)));
            UUU_TRANST{1,i}{1,j}(k,end) = gHT{i,j}(k,1);
        end
    end
    
end

TRADED_TIME = UUU_TRANST;

%% RULE for correcting the traded time
SIZZZE = size(UUU_TRANS{1,3}{1,2});
for i = 1 : 8
    for j = 1 : length(UUU_transT{1,1})
        for k = 1 : length(UUU_TRANST{1,5}{1,1})
            
            
            if  TRADED_TIME{1,i}{1,j}{k,1} >= UPDATED_TIME{1,i}{1,j}{k,1}
                
                TRADED_TIME{1,i}{1,j}{k,1} = UPDATED_TIME{1,i}{1,j}{k,1};
            end
                
            if TRADED_TIME{1,i}{1,j}{k,1} >= UPDATED_TIME{1,i}{1,j}{k,1}...
                    && TRADED_TIME{1,i}{1,j}{k,2} >= UPDATED_TIME{1,i}{1,j}{k,2}
                
                
                TRADED_TIME{1,i}{1,j}{k,1} = UPDATED_TIME{1,i}{1,j}{k,1};
                TRADED_TIME{1,i}{1,j}{k,2} = UPDATED_TIME{1,i}{1,j}{k,2};
            end
                
                
                
            if TRADED_TIME{1,i}{1,j}{k,1} >= UPDATED_TIME{1,i}{1,j}{k,1}...
                    && TRADED_TIME{1,i}{1,j}{k,2} >= UPDATED_TIME{1,i}{1,j}{k,2}...
                    && TRADED_TIME{1,i}{1,j}{k,3} >= UPDATED_TIME{1,i}{1,j}{k,3}
                
                
                TRADED_TIME{1,i}{1,j}{k,1} = UPDATED_TIME{1,i}{1,j}{k,1};
                TRADED_TIME{1,i}{1,j}{k,2} = UPDATED_TIME{1,i}{1,j}{k,2};
                TRADED_TIME{1,i}{1,j}{k,3} = UPDATED_TIME{1,i}{1,j}{k,3};
            end
                
                
            if TRADED_TIME{1,i}{1,j}{k,1} >= UPDATED_TIME{1,i}{1,j}{k,1}...
                    && TRADED_TIME{1,i}{1,j}{k,2} >= UPDATED_TIME{1,i}{1,j}{k,2}...
                    && TRADED_TIME{1,i}{1,j}{k,3} >= UPDATED_TIME{1,i}{1,j}{k,3}...
                    && TRADED_TIME{1,i}{1,j}{k,4} >= UPDATED_TIME{1,i}{1,j}{k,4}
                
                
                
                TRADED_TIME{1,i}{1,j}{k,1} = UPDATED_TIME{1,i}{1,j}{k,1};
                TRADED_TIME{1,i}{1,j}{k,2} = UPDATED_TIME{1,i}{1,j}{k,2};
                TRADED_TIME{1,i}{1,j}{k,3} = UPDATED_TIME{1,i}{1,j}{k,3};
                TRADED_TIME{1,i}{1,j}{k,4} = UPDATED_TIME{1,i}{1,j}{k,4};
                
            end
        end
    end
    
end

%% The traded time table
% Fill the updated time
for i = 1 : 8
    
   TABLE_TRADED_TIME{1,i} = vertcat(TRADED_TIME{1,i}{1,1},...
       TRADED_TIME{1,i}{1,2}, TRADED_TIME{1,i}{1,3},...
       TRADED_TIME{1,i}{1,4}, TRADED_TIME{1,i}{1,5},...
       TRADED_TIME{1,i}{1,6}, TRADED_TIME{1,i}{1,7},...
       TRADED_TIME{1,i}{1,8}, TRADED_TIME{1,i}{1,9},...
       TRADED_TIME{1,i}{1,10}, TRADED_TIME{1,i}{1,11},...
       TRADED_TIME{1,i}{1,12});
    
end

% The final concatenation
FINAL_TABLE_TRADED_TIME = vertcat(TABLE_TRADED_TIME{1,1},...
    TABLE_TRADED_TIME{1,2},TABLE_TRADED_TIME{1,3},...
    TABLE_TRADED_TIME{1,4},TABLE_TRADED_TIME{1,5},...
    TABLE_TRADED_TIME{1,6},TABLE_TRADED_TIME{1,7},...
    TABLE_TRADED_TIME{1,8});


%% The updated time table
% Fill the updated time
for i = 1 : 8
    
   TABLE_UPDATED_TIME{1,i} = vertcat(UPDATED_TIME{1,i}{1,1},...
       UPDATED_TIME{1,i}{1,2}, UPDATED_TIME{1,i}{1,3},...
       UPDATED_TIME{1,i}{1,4}, UPDATED_TIME{1,i}{1,5},...
       UPDATED_TIME{1,i}{1,6}, UPDATED_TIME{1,i}{1,7},...
       UPDATED_TIME{1,i}{1,8}, UPDATED_TIME{1,i}{1,9},...
       UPDATED_TIME{1,i}{1,10}, UPDATED_TIME{1,i}{1,11},...
       UPDATED_TIME{1,i}{1,12});
    
end

% The final concatenation
FINAL_TABLE_UPDATED_TIME = vertcat(TABLE_UPDATED_TIME{1,1},...
    TABLE_UPDATED_TIME{1,2},TABLE_UPDATED_TIME{1,3},...
    TABLE_UPDATED_TIME{1,4},TABLE_UPDATED_TIME{1,5},...
    TABLE_UPDATED_TIME{1,6},TABLE_UPDATED_TIME{1,7},...
    TABLE_UPDATED_TIME{1,8});

% Find the unique traded times, Sort the traded times and assign
% unique trade number to each in ascending order
UNIQ = sortrows(unique(cell2mat(FINAL_TABLE_TRADED_TIME),'rows'),[1 2 3 4]);
UNIQO = sortrows(unique(cell2mat(FINAL_TABLE_TRADED_TIME),'rows'),[1 2 3 4]);

% generate a random number
raNDgen1 = (324434 : 26835144)';
raNDgen2 = rand(length(raNDgen1),1);

% merge the two vectors
raNDgen = [raNDgen1, raNDgen2];
% sort according to the second column
Sort_raNDgen = sortrows(raNDgen,2);

% Cut at 2976 and sort according to the first column
Sort_raNDgen1 = Sort_raNDgen(1: length(FINAL_TABLE_TRADED_TIME), :);
Sort_raNDgen2 = sortrows(Sort_raNDgen1,1);

% assign the computed trade numbers to the corresponding trade times
UNIQO(:,5) = Sort_raNDgen2(:,1);
% size of UNIQO
SASS = size(UNIQO);
% finding the position of the sorted traded times in the original times
for i = 1 : length(FINAL_TABLE_TRADED_TIME)
POS{i,1} = num2cell(find(ismember(cell2mat(FINAL_TABLE_TRADED_TIME(:,1:end)),UNIQO(i,1:4),'rows')~=0));
end


%% THE_FINAL_TRADED_TIME
THE_FINAL_TRADED_TIME = zeros(length(FINAL_TABLE_TRADED_TIME), SASS(2));

for i = 1 : length(FINAL_TABLE_TRADED_TIME)
    
THE_FINAL_TRADED_TIME(cell2mat(POS{i,1}),:) = UNIQO(i,:);   
       
end

 Headers = OPriskDataSetexposure(1,:);

MATRIX = zeros(length(FINAL_TABLE_TRADED_TIME),length(Headers));

%% The traded time

% converting time into usal time formats
% there are 24 hours in the a day , to fins the hour
rt = 24.*THE_FINAL_TRADED_TIME(:,end-1);
hh = round(rt);

% the minutes
rr = 60.*abs(rt - hh);

mm = round(rr);

% the seconds
rg = 60.*abs(rr - mm);

ss = round(rg);

% Updated time as a vectors
Vec_tedTime = [THE_FINAL_TRADED_TIME(:,1:end-2), hh, mm, ss];

% converting the date back to string
formatOut = 'yyyy-mm-dd HH:MM:SS PM';
Vec_tradedTimeSTRING = datestr(Vec_tedTime(:, 1:end),formatOut);

%% The updated time
% converting time into usal time formats
% there are 24 hours in the a day , to fins the hour
rt = 24.*cell2mat(FINAL_TABLE_UPDATED_TIME(:,end));
hh = round(rt);

% the minutes
rr = 60.*abs(rt - hh);

mm = round(rr);

% the seconds
rg = 60.*abs(rr - mm);

ss = round(rg);

% Updated time as a vectors
Vec_updatedTime = [cell2mat(FINAL_TABLE_UPDATED_TIME(:,1:end-1)), hh, mm, ss];

% converting the date back to string
formatOut = 'yyyy-mm-dd HH:MM:SS PM';
Vec_updatedTimeSTRING = datestr(Vec_updatedTime(:, 1:end),formatOut);

%% generate the compatible columns
% capturedBy
UNI_STRINGS = unique(OPriskDataSetexposure(2:end,9));
% TraderID
UNI_STRINGS1 = unique(OPriskDataSetexposure(2:end,11));

for j = 1 : length(UNI_STRINGS)
    
    CapturedBy{j} = OPriskDataSetexposure(strcmp(OPriskDataSetexposure(:,9),...
        UNI_STRINGS(j))==1,9);
    % percentage proportion
    LEngC(j) = length(CapturedBy{j})./length(OPriskDataSetexposure);
    format long
   N_STR(j) = ceil(LEngC(j).* length(THE_FINAL_TRADED_TIME));
   
   N_STRRR{j} = num2cell(zeros(N_STR(j),1));
   
end


for j = 1 : length(UNI_STRINGS)
   N_STRRR{j}(:,1) = (UNI_STRINGS(j,1));
end
%
CAPTUREDBY_TOTAL = vertcat(N_STRRR{1,1},N_STRRR{1,2},...
    N_STRRR{1,3},N_STRRR{1,4},N_STRRR{1,5});

CAPTUREDBY_TOTAL = CAPTUREDBY_TOTAL(1:length(THE_FINAL_TRADED_TIME));
CAPTUREDBY_TOTAL = [CAPTUREDBY_TOTAL, num2cell(rand(length(CAPTUREDBY_TOTAL),1))];
CAPTUREDBY_TOTAL = sortrows(CAPTUREDBY_TOTAL,2);
%%


%% generate the compatible columns
% TraderID
Tra_UNI_STRINGS = unique(OPriskDataSetexposure(2:end,11));
% TraderID
Tra_UNI_STRINGS1 = unique(OPriskDataSetexposure(2:end,11));

for j = 1 : length(Tra_UNI_STRINGS)
    
    TraderID{j} = OPriskDataSetexposure(strcmp(OPriskDataSetexposure(:,11),...
        Tra_UNI_STRINGS(j))==1,11);
    % percentage proportion
    LEngC(j) = length(TraderID{j})./length(OPriskDataSetexposure);
    format long
   TRAN_STR(j) = ceil(LEngC(j).* length(THE_FINAL_TRADED_TIME));
   
   TRAN_STRRR{j} = num2cell(zeros(TRAN_STR(j),1));
   
end


for j = 1 : length(Tra_UNI_STRINGS)
   TRAN_STRRR{j}(:,1) = (Tra_UNI_STRINGS(j,1));
end

TRADERID_TOTAL = vertcat(TRAN_STRRR{1,1},TRAN_STRRR{1,2},...
    TRAN_STRRR{1,3},TRAN_STRRR{1,4},TRAN_STRRR{1,5},...
    TRAN_STRRR{1,6},TRAN_STRRR{1,7});

 TRADERID_TOTAL = TRADERID_TOTAL(1:length(THE_FINAL_TRADED_TIME));
TRADERID_TOTAL = [TRADERID_TOTAL, num2cell(rand(length(TRADERID_TOTAL),1))];
TRADERID_TOTAL = sortrows(TRADERID_TOTAL,2);
%% Business lines

BL1 = [cellstr('BL1') cellstr('BL1') ;...
    cellstr('Credit Derivatives') cellstr('Investment Banking')]';

BL2 = [cellstr('BL2') cellstr('BL2') cellstr('BL2') cellstr('BL2')...
    cellstr('BL2') cellstr('BL2') cellstr('BL2') cellstr('BL2')...
    cellstr('BL2') cellstr('BL2') cellstr('BL2') cellstr('BL2')...
    cellstr('BL2');...
    cellstr('Rates') cellstr('MM') cellstr('Equity')...
    cellstr('Commodities') cellstr('Africa')...
    cellstr('Options') cellstr('Bonds/Repos')...
    cellstr('Forex') cellstr('Prime Services')...
    cellstr('Credit Derivatives') cellstr('Management')...
    cellstr('Group Treasury') cellstr('SND')]';

BL3 = [cellstr('BL3') cellstr('BL3') cellstr('BL3') ;...
    cellstr('Africa') cellstr('MM') cellstr('SND')]';

BL4 = [cellstr('BL4') cellstr('BL4') cellstr('BL4')...
     cellstr('BL4') cellstr('BL4') cellstr('BL4');...
    cellstr('ACBB') cellstr('Credit Derivatives') cellstr('Funding')...
     cellstr('MM') cellstr('Portfolio Management') cellstr('SND')]';
 
 BL5 = [cellstr('BL5') cellstr('BL5') ;...
    cellstr('Credit Derivatives') cellstr('MM')]';

BL6 = [cellstr('BL6') cellstr('BL6') ;...
    cellstr('Management') cellstr('Prime Services')]';

BL7 = [cellstr('BL7') cellstr('BL7') ;...
    cellstr('Portfolio Management') cellstr('SND')]';

BL9 = [cellstr('BL9') cellstr('Portfolio Management')];

%% generate the compatible columns
% Business line
BUB_UNI_STRINGS = unique(OPriskDataSetexposure(2:end,22));
%
for j = 1 : length(BUB_UNI_STRINGS)
    
    Bus{j} = OPriskDataSetexposure(strcmp(OPriskDataSetexposure(:,22),...
        BUB_UNI_STRINGS(j))==1,22);
    % percentage proportion
    LEngB(j) = length(Bus{j})./length(OPriskDataSetexposure);
    format long
    Bu_STR(j) = ceil(LEngB(j).* length(THE_FINAL_TRADED_TIME));
    
    BUs_STRRR{j} = num2cell(zeros(Bu_STR(j),1));
    
end


for j = 1 : length(BUB_UNI_STRINGS)
    BUs_STRRR{j}(:,1) = (BUB_UNI_STRINGS(j,1));
end
%
BUS_TOTAL = vertcat(BUs_STRRR{1,1},BUs_STRRR{1,2},...
    BUs_STRRR{1,3},BUs_STRRR{1,4},BUs_STRRR{1,5},...
    BUs_STRRR{1,6}, BUs_STRRR{1,7}, BUs_STRRR{1,8});

BUS_TOTAL = BUS_TOTAL(1:length(THE_FINAL_TRADED_TIME));
BUS_TOTAL = [BUS_TOTAL, num2cell(rand(length(BUS_TOTAL),1))];
BUS_TOTAL = sortrows(BUS_TOTAL,2);

%% BUSINESS LINES

BUSINESS_LINESL = [BL1;BL2;BL3;BL4;BL5;BL6;BL7;BL9];

for j = 1 : length(THE_FINAL_TRADED_TIME)
    
    BUSINESS_LINES{j} = num2cell((j.^0).*zeros(length(BUSINESS_LINESL),3));
    BUSINESS_LINES{j}(:,3) = num2cell((j.^0).*rand(length(BUSINESS_LINESL),1));
    BUSINESS_LINES{j}(:,1:2) = BUSINESS_LINESL(:,1:2);
    
    POSITION{j} = find(strcmp(BUSINESS_LINES{j}(:,1),...
        BUS_TOTAL(j,1))==1, 1, 'last' );
    % the desk
    DESK(j,1) = BUSINESS_LINESL(POSITION{j},2);
end
%%
% fill in the matrix
MATRIX(:,1) = (THE_FINAL_TRADED_TIME(:,end));
MATRIX(:,7) = (THE_FINAL_TRADED_TIME(:,end-1));
MATRIX(:,6) = (THE_FINAL_TRADED_TIME(:,end-2));
MATRIX(:,4) = cell2mat(FINAL_TABLE_UPDATED_TIME(:,end));
MATRIX(:,3) = cell2mat(FINAL_TABLE_UPDATED_TIME(:,end-1));


MATRIX = num2cell(MATRIX);

MATRIX(:,8) = DESK(:,1);
MATRIX(:,22) = BUS_TOTAL(:,1);
MATRIX(:,9) = CAPTUREDBY_TOTAL(:,1);
MATRIX(:,11) = TRADERID_TOTAL(:,1);
% fill in the updated time and the traded time
MATRIX(:,2) = cellstr(Vec_updatedTimeSTRING);
MATRIX(:,5) = cellstr(Vec_tradedTimeSTRING);


% % Exporting the results to Excel
% filename = 'Hoohlo.xlsx';
% writetable(cell2table(MATRIX ,...
%     'VariableNames',Headers),...
%     filename,'Sheet',1,'Range','A1')

```
\normalsize

\subsection{Deploying the R Model.}

Import extrapolated data (1 year data), use the estimated model \texttt{Amodel} to predict forecasted estimates

\small

```r
 newfname <- "file:///C:/Users/Mphekeleli/Documents/R PROJECT/OpRiskPHDGitHub
/OpRisk_PHD_Thesis/Data/Extrap_Data_Model.csv"

 crs$newdataset <- read.csv(newfname,
                sep=",",
                dec=".",
                na.strings=c(".", "NA", "", "?"),
                strip.white=TRUE, encoding="UTF-8",header=T)

 head(crs$newdataset)

 Forecast <- predict(Amodel, crs$newdataset, type = "link")
 Forecast
```
\normalsize

\clearpage

\section{Appendix C: R Code for GAMLSS}
\label{sec:Appendix B: R Code for GAMLSS}

\singlespace

Required: R Packages from CRAN (in addition to packages already found in Chapter 3)

\small

```r
if (!require(gamlss)){
  install.packages("gamlss")
  library(gamlss)
}
```
\normalsize

\subsection{Data exploration of OpRisk loss severity dataset in preparation for GAMLSS machine learning treatment. The raw data and the pre-processed datasets are initiated and analysed.}
\label{ssec:Data exploration GAMLSS}

Plots exploring loss severity characteristics, see chapter \ref{EXPOSURE-BASED OPERATIONAL RISK ANALYSIS} section \ref{sec:The estimation of some  generalised additive models for location scale and shape (GAMLSS) for severity loss estimation} on page \pageref{sec:The estimation of some  generalised additive models for location scale and shape (GAMLSS) for severity loss estimation} 

\small

```r
options(scipen = 999)
file_loc <- "C:/Users/Mphekeleli/Documents/R PROJECT/OpRiskPHDGitHub/
OpRisk_PHD_Thesis/Data"
setwd(file_loc)
list.files(file_loc)

frequency <- openxlsx::read.xlsx("Raw_Formatted_Data.xlsx",
                                 check.names = TRUE, sheet = "Frequency")
severity <- openxlsx::read.xlsx("Raw_Formatted_Data.xlsx",
                                check.names = TRUE, sheet = "Severity")
projdata <- openxlsx::read.xlsx("OPriskDataSet_GAMLSS.xlsx",
                                check.names = TRUE, sheet = "CleanedData")

names(projdata) <- sub("\\.", "", names(projdata))
dput(names(projdata))

par(mar=c(1,1,1,1))
PPP <- par(mfrow=c(2,2))
plot(Loss ~ UpdateTime, data = projdata, col=gray(0.7), pch=15, cex=0.5)
plot(Loss ~ UpdatedDay, data = projdata, col=gray(0.7), pch=15, cex=0.5)
par(PPP)
```
\normalsize

\subsection{Data partitioning of the pre-processed OpRisk dataset into Training/Validation/Testing proportions, in preparation for machine learning model building treatments. The original dataset is partitioned into three random subsets initiated by a random number sequence with a randomly selected seed.}
\label{ssec:GAMLSS Data Training/Validation/Testing}

\small

```r
# Load packages
library(rattle, quietly = TRUE)
library(magrittr, quietly = TRUE) # Utilize the %>% and %>% pipeline operators
library(Hmisc, quietly = TRUE)
library(chron, quietly = TRUE)
library(dplyr, quietly = TRUE)
library(ggplot2)
library(caTools)
library(caret)
library(gamlss)

building <- TRUE
scoring <- ! building

# A predefined value is used to reset the random seed so that results are repeatable

crv$seed <- 42 # set random seed to make your partition reproducible
# Load the dataset OPriskDataSet_exposure
#===================================================================================
fname <- "file:///C:/Users/Mphekeleli/Documents/R PROJECT/OpRiskPHDGitHub/
OpRisk_PHD_Thesis/Data/OPriskDataSet_exposure_severity.csv"
crs$dataset <- read.csv(fname,
                        sep=",",
                        dec=".",
                        na.strings=c(".", "NA", "", "?"),
                        strip.white=TRUE, encoding="UTF-8")

exposure <- crs$dataset[,ncol(crs$dataset)]
crs$dataset <- as.data.frame(crs$dataset)

# The following varaible selections have been noted
crs$input <- crs$dataset %>%
  group_by(UpdatedDay,
           UpdatedTime,
           TradedDay,
           TradedTime,
           Desk,
           CapturedBy,
           TradeStatus,
           TraderId,
           Instrument,
           Reason,
           EventTypeCategoryLevel1,
           BusinessLineLevel1) %>%
  transmute(LossesIndicator = LossIndicator,
            Losses = Loss,
            Exposure = exposure)

# Create function "getmode" which finds the modal class in the categorical variables
getmode <- function(x){
  u <- unique(x)
  as.integer(u[which.max(tabulate(match(x,u)))])
}
# Reorder the categorical variables so that the modal class
# is specified as the reference level
for (i in 5:(ncol(crs$input) - 3)){
  crs$input[[i]] <- relevel(crs$input[[i]], getmode(crs$input[[i]]))
}

# Build the training/validation/testing datasets
# nobs=2331 training=1632 validation=350 testing=349

set.seed(crv$seed)

crs$nobs <- nrow(crs$input)

crs$train <- sample(crs$nobs, 0.7*crs$nobs)

crs$nobs %>%
  seq_len() %>%
  setdiff(crs$train) %>%
  sample(0.15*crs$nobs) ->
  crs$validate

crs$nobs %>%
  seq_len() %>%
  setdiff(crs$train) %>%
  setdiff(crs$validate) ->
  crs$test


crs$training <- as.data.frame(crs$input[crs$train,])
crs$validation <- as.data.frame(crs$input[crs$validate,])
crs$testing <- as.data.frame(crs$input[crs$test,])
```
\normalsize 

\subsection{Selection of the distribution.}
\label{sec:Selection of the distribution.}

Begin to model the response Losses (y) using the three parameter Zero adjusted gamma (ZAGA) distribution, followed by the Zero adjusted inverse gamma (ZAIG) distribution, and lastly using the four parameter Generalized beta type 2 ($\mbox{GB}2$) distribution. Start by fitting the full linear model $\mu$ including all explanatory variables. 

\small

```r
mod1 <- gamlss(Losses ~ UpdatedDay + Desk + CapturedBy + TradeStatus + TraderId 
     + Instrument + Reason + EventTypeCategoryLevel1 + BusinessLineLevel1 + Exposure,
     mu.start = NULL,  sigma.start = NULL, nu.start = NULL, tau.start = NULL,
     sigma.fo = ~1, nu.fo = ~1, data=crs$training,  family = ZAGA, n.cyc=80)
mod1

mod2 <- gamlss(Losses ~ UpdatedDay + Desk + CapturedBy + TradeStatus + TraderId
     + Instrument + Reason + EventTypeCategoryLevel1+ BusinessLineLevel1 + Exposure,
     mu.start = NULL,  sigma.start = NULL, nu.start = NULL, tau.start = NULL,
     sigma.fo = ~1, nu.fo = ~1, data=crs$training,  family = ZAIG, n.cyc=80)
mod2

mod3 <- gamlss(Losses ~ UpdatedDay + Desk + CapturedBy + TradeStatus + TraderId 
     + Instrument + Reason + EventTypeCategoryLevel1 + BusinessLineLevel1 + Exposure,
     mu.start = NULL,  sigma.start = NULL, nu.start = NULL, tau.start = NULL,
     sigma.fo = ~1, nu.fo = ~1, tau.fo = ~1, data=crs$training,  family = GB2, n.cyc=1700)
mod3

GAIC(mod1, mod2, mod3)
```
\normalsize

\subsection{Selection of terms.}
\label{sec:Selection of terms.}

Now we use \emph{\texttt{drop1()}} to check whether any linear terms can be dropped and then using \emph{\texttt{add1()}}, we consider adding a two-way interaction term into the linear model mod1 

\small

```r
drop1(mod3)
add1(mod3, scope=~(UpdatedDay + Desk + CapturedBy + TradeStatus + TraderId +
Instrument + Reason + EventTypeCategoryLevel1 + BusinessLineLevel1 + Exposure)^2)
```
\normalsize

We use \emph{\texttt{FORM}} as an upper bound for \emph{\texttt{scope}}, starting from all explanatory terms sos that all interactions are considered

\small

```r
FORM <- as.formula("~(UpdatedDay + Desk + CapturedBy + TradeStatus + TraderId
+ Instrument + Reason + EventTypeCategoryLevel1 + BusinessLineLevel1 + Exposure)^2
+pb(UpdatedDay) + pb(Desk) + pb(CapturedBy) + pb(TradeStatus)+ pb(TradeStatus)
+ pb(TraderId) + pb(Instrument) + pb(Reason)+ pb(EventTypeCategoryLevel1) 
+ pb(BusinessLineLevel1) + pb(Exposure)")
```
\normalsize

\subsection{The estimation of some GAMLSS for OpRisk loss severity distribution: To build the model we pass on to the model building function \texttt{gamlss} i.e., the formula that describes the model to build.}
\label{sec:The estimation of some GAMLSS for OpRisk loss severity distribution}

\small

```r
mod14 <- stepGAICAll.A(mod1, scope=list(lower=~1, upper=FORM), k=log(371))

GAIC(mod14,mod24,mod34)
mod14
plot(mod14)
mod14$anova
```
\normalsize

\begin{figure}
\centering
\begin{subfigure}[b]{0.5\textwidth}
   \includegraphics[width=1\linewidth]{GAIC_statistics.eps}
   \caption{}
   \label{Residuals summary statistics} 
\end{subfigure}

\begin{subfigure}[b]{0.5\textwidth}
   \includegraphics[width=1\linewidth]{ROC_4Plot_GB2.eps}
   \caption{}
   \label{Residuals_GB2}
\end{subfigure}

\caption[Normalized quantile residuals from model GB2]{(a)Summary statistics of quantile residuals from models GB2, ZAIG \& ZAGA (b)Displays (normalized quantile) residuals from model $\mbox{GB}2(\mu,\sigma,\nu,\tau)$.}
\end{figure}

\begin{figure}
\centering
\begin{subfigure}[b]{0.5\textwidth}
   \includegraphics[width=1\linewidth]{ROC_4Plot_ZAIG.eps}
   \caption{}
   \label{Residuals_ZAIG} 
\end{subfigure}

\begin{subfigure}[b]{0.5\textwidth}
   \includegraphics[width=1\linewidth]{ROC_4Plot_ZAGA.eps}
   \caption{}
   \label{Residuals_ZAGA}
\end{subfigure}

\caption[Normalized quantile residuals from models ZAIG \& ZAGA]{(a) Display of (normalized quantile) residuals from model $\mbox{ZAIG}(\mu,\sigma,\nu)$ (b) As for (a) but on the model $\mbox{ZAGA}(\mu,\sigma,\nu)$}
\end{figure}
