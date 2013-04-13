READ ME:

I. Intro:
The Generalized Treatment Effect Toolbox(GTET) is an on going effort to facilitate program evaluation in Python environment.
The vision of the GTET project is to provide an array of user friendly treatment effect estimators that have sophisticated underlying logic.

II. What's in Your Box?

The current version is a replica of Heckman, Urzua and Vytlacil (2006), whose detailed documentation can be found at "http://jenni.uchicago.edu/underiv/".

As a pre-requisite, you need to have the following libraries:
(i) numpy
(ii)Parallel Python(pp)
(iii)Python Spatial Analysis Library(pysal)  - [needed for Probit model, which I am too lazy to write myself]
(iv) os - needed to run my example

The current version offers the following four estimator for Generalized Roy Model:
(1) Linear in outcome equation, normal distribution for the error term
(2) Linear in outcome equation, polynomial estimation for MTE
(3) Nonparametric in outcome equation, polynomial estimation for MTE
(4) Nonparametric in outcome equation, Local Linear Regression estimation for MTE (LIV)

The current version has the following limitation:
(a) The covariates in outcome equation has be specified separately from the covariates in the choice equation.
(b) The weights is calculated by population average, rather than conditioning on specific X0
(c) Method (3) and (4) is computational intensive and requires further optimization.
(d) The propensity score is estimated with a probit for all 4 methods.

III. How to Use it?
It is as simple as it gets. For all four flavors, see the estimateExampleScript.py

Step 0: build the obj
Actual Code : " testEstObj = ei.treatmentEffectEst() "


Step 1:
You need to supply the toolbox with D,Z for the choice equation and X,D,Y for the outcome equation., 
where D is the treatment variable, Y is the outcome variable, X and Z are covariates.

Actual Code : " testEstObj.assignData(Y, D, X, Z) "


Step 2: 
You need to supply the toolbox with model speficiation. Currently you need to specify 4 parameters

				Method (1) Method(2) method(3) method(4)
(1) Linear			1			1		0		0
(2) Normal			1			0		0		0
(3) LIVMTE			*			*		0		1
(4) PolyOrder			*			k		k		*

where "*" means any integer value and "k" means the integer order of polynomial

Actual Code : " testEstObj.assignParam(linear, normal, LLRMTE, polyOrder) "

Step 3:
Call the relevant method of the estimator object with the x0 to be conditioned on. The three treatment effect methods available are:
(1) ate
(2) att
(3) atut

Notice that you need to call the lower case method.

Actual Code : " testEstObj.ate(x0) "

IV. What to Expect?

The 2nd release is going to be an expansion for matching estimator. I plan to release 
(1) Inverse Probability Weighting 
(2) Kernel Matching
(3) Nearest K neighborhood matching

all of which are matched on propensity score.

Therefore, the 3rd release will add more propensity score estimator other than probit model. I plan to release
(1) Logit/Probit
(2) Polynomial logit
(3) Polynomial logit on steroid (RKHS version)
(4) Normal Mixture (MLE - EM Algorithm)

V. Join the Party!
If you want to contribute to the project, please do not hesistate to contact me.
If you find a bug in the toolbox, please do not hesistate to contact me.
If you have specific request for the toolbox, please do not hesistate to contact me.

My email is Junchen@uchicago.edu

Reference:
James J. Heckman, Sergio Urzua and Edward Vytlacil, "Understanding Instrumental Variables in Models with Essential Heterogeneity",The Review of Economics and Statistics, Vol. 88, No. 3 (Aug., 2006), pp. 389-432

