\name{girtest}
\alias{girtest}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Confidence Intervals for the Strength of Identification
%%  ~~function to do ... ~~
}
\description{This function computes confidence intervals for the concentration parameter, the bias of the 2SLS estimator and for the size distortion of the associated Wald statistic (5\%nominal level). Proposition 1 in Ganics, Inoue and Rossi (2018).
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
girtest(Fstat, K_2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Fstat}{F-statistic of the first-stage regression on the instruments (excluded exogenous variables).
%%     ~~Describe \code{x} here~~
}
  \item{K_2}{Number of instruments (excluded exogenous variables) used in the first-stage regression.}
}
\details{Note it is only applicable when there is only one endogenous regressor, and under homoskedasticity and serial uncorrelation.
%%  ~~ If necessary, more details than the description above ~~
}
\value{Three 95\% asymptotic level confidence intervals: for the concentration parameter; for the bias of the 2SLS estimator; and for the size distortion of the associated Wald statistic. The critical values are computed using the tables in the Appendix D in Ganics, Inoue and Rossi (2018), Table D.7 (bias) and Table D.10 (size distortion).
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{Ganics, Inoue and Rossi (2018)
%% ~put references to the literature/web site here ~
}
\author{Oriol Gonzalez and Marko Irisarri
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (x)
{ Example of use: Replication of Angrist and Krueger (1991)/Bound et al. (1995)
  Angrist and Krueger (1991): 28 instruments and F-stat = 1.61 --> girtest(Fstat = 1.61, K_2 = 28)
  Bound et al. (1995): 3 instruments and F-stat = 13.49        --> girtest(Fstat = 13.49, K_2 = 3)

  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
