\name{switchnpreg}
\alias{switchnpreg}

\title{
  Fit a switching nonparametric regression model 
}

\description{
  Estimate the parameters of a switching nonparametric regression model
  using the EM algorithm as proposed by De Souza and Heckman (2013). The
  package allows two different estimation approaches (Bayesian and
  penalized log-likelihood) and two different types of hidden states
  (iid and Markov). The smoothing parameters are chosen by
  cross-validation. Standard errors for the estimates of the parameters
  governing the distribution of the state process are also provided.
}

\usage{
switchnpreg(x, y, f, alpha, sigma2, lambda, ...,
            method = c("pl", "bayes"), var.equal = TRUE,
            z.indep = TRUE, eps.cv, eps.em, maxit.cv, maxit.em)
}

\arguments{
  \item{x}{
    The sequence of covariates \eqn{x_1, \dots, x_n}.
  }
  \item{y}{
    The sequence of response variables \eqn{y_1, \dots, y_n}.
  }
  \item{f}{
    The \eqn{n \times J} matrix of initial values for the functions,
    where column \eqn{j} corresponds to the function \eqn{f_j}.
  }
  \item{alpha}{
    The initial values for the parameters of the latent state process.
    If the latent states are iid alpha is a vector containing the initial
    mixing proportions \eqn{p_j} for \eqn{j=1,\dots,J}. If the latent
    states follow a Markov structure then alpha is a list of two
    components: A and PI, where A is the initial \eqn{J \times J} matrix
    of transition probabilities A and PI is the initial \eqn{J}-vector
    of initial probabilities.
  }
  \item{sigma2}{
    The initial \eqn{J}-vector of regression error variances.
  }
  \item{lambda}{
   The initial \eqn{J}-vector of smoothing parameters.
  }
  \item{\dots}{
    Optional arguments to parameter update functions.
  }
  \item{method}{
    Character string 'pl' or 'bayes' to choose whether the model is
    fitted using the penalized log-likelihood approach or the Bayesian
    approach, respectively. 
  }
  \item{var.equal}{
    Logical indicating whether \eqn{\sigma^2_j} are equal for all \eqn{j}.
  }
  \item{z.indep}{
    Logical indicating whether the hidden states \eqn{z_i, \dots, z_n}
    are considered iid or Markovian.
  }
  \item{eps.cv}{
    Convergence value for the cross-validation procedure.
  }
  \item{eps.em}{
    Convergence value for the EM algorithm.
  }
  \item{maxit.cv}{
    Maximum number of iterations of the EM+CV procedure.
  }
  \item{maxit.em}{
    Maximum number of iterations of each EM loop.
  }
}

\value{
  A list with following elements:
  \item{current}{The final estimate of \eqn{\theta}, represented as a
    list with the elements named after the respective model parameter:
    \describe{
      \item{f}{The final function estimates.}
      \item{sigma2}{The final variance estimates.}
      \item{alpha}{The final estimates for the parameters of the latent
        state process.}
      \item{pij}{The matrix of size \eqn{n \times\ J} with \eqn{ij}-th
        element giving the final estimate of \eqn{p(z_i=j|y,\theta)}.}
    }
  }
  
  \item{lambda}{Chosen smoothing parameters.}
  \item{iter.cv}{Number of iterations of the EM+CV procedure.}
  \item{stderr}{Standard errors for the parameter estimates of the
    latent state process.}
}

\references{
  de Souza and Heckman (2013), \dQuote{Switching nonparametric
  regression models and the motorcycle data revisited}, submitted for
  peer review. Available at
  \href{http://arxiv.org/abs/1305.2227v2}{arXiv.org}, article-id:
  arXiv:1305.2227v2.
}

\author{
  Camila de Souza \email{camila@stat.ubc.ca} and Davor Cubranic
  \email{cubranic@stat.ubc.ca}.
}

\seealso{
  \code{demo(simulated_data_indep_example)},
  \code{demo(simulated_data_Markov_example)}
}

\examples{
  
  ## The motorcycle data set revisited ##
  
  x <- MASS::mcycle$times
  set.seed(30)
  x[duplicated(x)] <- round(jitter(x[duplicated(x)]),3)
  
  y <- MASS::mcycle$accel
  
  n <- length(y)
      
  spline_fit <- smooth.spline(x, y)
  
  ## set up the initial functions
  f.initial <- t(apply(as.matrix(spline_fit$y), 1,
                       `+`, c(30, 0, -30)))
  J <- ncol(f.initial)
  sig2 <- rep((sum((y-predict(spline_fit, x)$y)^2) / (n - spline_fit$df))/J, J)

  ## B and R parameters for penalized log-likelihood method
  basis <- create.bspline.basis(range(x), nbasis = 40)
  B <- getbasismatrix(x, basis)
  R <- getbasispenalty(basis)
    
  estimates <- switchnpreg(x = x, y = y,
                           f = f.initial,
                           alpha = rep(1, J) / J,
                           sigma2 = sig2,
                           lambda = rep(.5, J),
                           B = B, R = R,
                           var.equal = FALSE,
                           interval = log(c(1E-4, 1E3)),
  
                           eps.cv = rep(1E-1, J),
                           eps.em = rep(c(1E-1, 1E-2, 1E-3), each = J),
                           maxit.cv = 10,
                           maxit.em = 100)
  
  plot(x, y, ylim = c(-150,90),
       ylab = 'Head acceleration',
       xlab = 'Time')
  matlines(x, estimates$current$f, type='l', lty = 1, col = 1:J)
  matlines(sort(x), f.initial, lty = 2, col = 'gray')
}

\keyword{ nonparametric }
\keyword{ regression }
\keyword{ switching }
\keyword{ latent }
