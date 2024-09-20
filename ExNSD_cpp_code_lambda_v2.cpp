// A fully Bayesian implementation of nonparametric statistical downscaling (incorporating thinning and 
// using an efficient algorithm for the Cholesky decompositions):

#include <cmath>
#include <RcppArmadillo.h>
#include <iostream>

//[[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

// [[Rcpp::export]]
arma::mat calcDistsArma(arma::mat coordsArma){
  NumericMatrix coords = wrap(coordsArma);
  int nRows = coords.nrow();
  double d;
  NumericMatrix dists(nRows,nRows);
  int i,j;
  
  for(i=0; i<nRows; i++){
    for(j=i+1; j<nRows; j++){
      NumericVector v1 = coords.row(i);
      NumericVector v2 = coords.row(j);
      NumericVector v3 = v1-v2;
      d = sqrt(sum(pow(v3,2)));
      dists(j,i)=d;
      dists(i,j)=d;
    }
  }
  arma::mat distsArma(dists.begin(),dists.nrow(),dists.ncol(),false);
  return distsArma;
}


// [[Rcpp::export]]
arma::mat H0functionArma(double phi, arma::mat dists){
  arma::mat Ha = -phi*dists;
  arma::mat H0 = arma::exp(Ha);
  return H0;
}


// [[Rcpp::export]]
arma::mat choleskyEfficient(arma::mat Amat){
  int n = Amat.n_rows;
  arma::vec z(n);
  for(int l=0; l<n; l++){
    z(l) = arma::min(arma::find(Amat.col(l)));
  }
  int dimA = n*n;
  NumericVector L(dimA);
  NumericVector A = wrap(arma::vectorise(Amat));
  for (int i = 0; i < n; i++){
    for (int j = z(i); j < (i+1); j++){
      double s = 0;
      for (int k = 0; k < j; k++){
        s += L(i * n + k) * L(j * n + k);
      }  
      if(i == j){
        L(i * n + j) = sqrt(A(i * n + i) - s);
      }else{
        L(i * n + j) = (1.0 / L(j * n + j) * (A(i * n + j) - s));
      }
    }
  }
  arma::vec LArma(L.begin(),L.size());
  arma::mat LArmaMat(LArma);
  arma::mat LArmaMat2 = arma::reshape(LArmaMat,n,n);
  return LArmaMat2;
}

// [[Rcpp::export]]
double normalPDF(double x, double mean, double variance) {
  double exponent = -0.5 * std::pow((x - mean) / std::sqrt(variance), 2);
  double coefficient = 1.0 / (std::sqrt(2 * M_PI * variance));
  return coefficient * std::exp(exponent);
}

// [[Rcpp::export]]
double mvnPDF(const arma::vec& x, const arma::vec& mean, const arma::mat& covariance) {
  // Check if dimensions match
  if (x.n_elem != mean.n_elem || x.n_elem != covariance.n_rows || covariance.n_rows != covariance.n_cols) {
    std::cerr << "Error: Input dimensions do not match!" << std::endl;
    return 0.0;
  }
  
  // Constants
  const double pi = arma::datum::pi;
  
  // Calculate the exponent term in the PDF
  arma::vec diff = x - mean;
  arma::mat invCovariance = inv(covariance); // Inverse of covariance matrix
  double exponent = -0.5 * as_scalar(diff.t() * invCovariance * diff);
  
  // Calculate the normalization term
  double detCovariance = det(covariance);
  double normalization = 1.0 / (std::pow(2.0 * pi, x.n_elem / 2.0) * std::sqrt(detCovariance));
  
  // Calculate the PDF
  double pdf = normalization * std::exp(exponent);
  
  return log(pdf);
  
}

// [[Rcpp::export]]
double sumWithoutNaN(const arma::vec& v)
{
  double sum = 0.0;
  for (size_t i = 0; i < v.n_elem; ++i)
  {
    if (!std::isnan(v(i))) // Check if element is not NaN
    {
      sum += v(i);
    }
  }
  return sum;
}

// [[Rcpp::export]]
double invlogit(double x)
{
  // Compute the inverse logit transformation
  double expNegX = std::exp(-x);
  double invLogit = 1.0 / (1.0 + expNegX);
  if(invLogit==1.0){
    invLogit=0.999;
  }
  if(invLogit==0.0){
    invLogit=0.001;
  }
  return invLogit;
}

// [[Rcpp::export]]
double diracgpdPDF(double x_value, double xi, double sigma, double p) {
  double pdf_value;
  double z = x_value/sigma;
  if(x_value == 0) {
    pdf_value = (1-p);
  } else {
    if(xi >= 0) {
      if(xi == 0 || abs(xi) < 1e-10){
        pdf_value = p*1/sigma * std::exp(-1.0 * z);
      } else{
        pdf_value = p*1/sigma * std::pow(1.0 + (xi * z), -1.0 / xi - 1.0);
      }
    } else {
      if(x_value <= -sigma/xi) {
        pdf_value = p*1/sigma * std::pow(1.0 + (xi * z), -1.0 / xi - 1.0);
      } else{
        pdf_value = 1e-200;
      }
    }
  }
  return pdf_value;
}

// [[Rcpp::export]]
double roundToSig(double x, int n) {
  if (x == 0.0) {
    return 0.0; // Early return for zero to avoid log of zero
  }
  // Calculate the factor needed to shift all significant digits desired into the integer part
  double factor = std::pow(10.0, n - 1 - std::floor(std::log10(std::fabs(x))));
  // Scale the number, round it, and then scale back
  return std::round(x * factor) / factor;
}

// [[Rcpp::export]]
double sample_diracDelta(double sigma, double xi, double p) {
  double r = arma::randu();  // Generate a uniform random number
  if (r < p) {
    double rs = arma::randu();  // Generate another uniform random number for GPD
    // Implementing the quantile function of the Generalized Pareto Distribution
    if (xi == 0) {
      return sigma * log(1 / (1 - rs));  // Limiting case for xi -> 0
    } else {
      return sigma * (pow(1 - rs, -xi) - 1) / xi;
    }
  } else {
    return 0.0;  // Return 0 with probability 1-p
  }
}

// [[Rcpp::export]]
double laplace_trans_pdf(double x, double b, double mu) {
  double output = 0.0;
  // Calculate bottom part of the PDF
  double bottom = 0.5 * (2 - exp((-0.5 - mu) / b) - exp((mu - 0.5) / b));
  // Check if x is within the range [-0.5, 0.5]
  if (x >= -0.5 && x <= 0.5) {
    // Calculate top part of the PDF
    double top = (1 / (2 * b)) * exp(-abs(x - mu) / b);
    // Compute the PDF value
    output = log(top / bottom);
  } else {
    // Return a very small value instead of 0 to avoid numerical issues
    output = log(1e-200);
  }
  
  return output;
}

// [[Rcpp::export]]
arma::vec sampleFromNormal(int n, double mean, double sd) {
  arma::vec samples = arma::randn(n);  // Generate n standard normal random numbers
  return samples * sd + mean;          // Scale and shift to desired mean and sd
}

// [[Rcpp::export]]
int sampleFromRange(int nData) {
  // Seed for random number generation
  std::random_device rd;  
  std::mt19937 gen(rd()); 
  
  // Create a uniform integer distribution from 1 to nData
  std::uniform_int_distribution<> dis(0, nData-1);
  
  // Sample a single integer from the distribution
  return dis(gen);
}


// [[Rcpp::export]]
double binomPDF(int x, double p) {
  if (x == 1) {
    return p;
  } else if (x == 0) {
    return 1 - p;
  } else {
    Rcpp::stop("x must be 0 or 1");
    return 0; // This line won't be reached, just to suppress compiler warning
  }
}
// [[Rcpp::export]]
NumericMatrix ExNSDmodel(int nIter, int nThin, List yData, List xData, List xPred, NumericMatrix coordsData, NumericMatrix coordsPred, List By, List Bx, 
                         NumericMatrix ByPred, List BxPred, double phiAlpha, double phiBeta,
                         double aAlpha, double bAlpha, double aBeta, double bBeta, //double aY, double bY, 
                         double aC, double bC, NumericVector muD, NumericMatrix SigmaD,//double aX, double bX, 
                         double sigmaAlphaPrecInit, double sigmaBetaPrecInit,double sigmaCPrecInit,
                         NumericMatrix cInit,  NumericMatrix alphaInit, NumericMatrix betaInit, NumericMatrix dInit,//double sigmaYPrecInit, double sigmaXPrecInit,
                         double xiXInit, double xiYInit, 
                         int uniqueXi,  
                         double sigmaATune, double sigmaBTune, double sigmaCTune, double sigmaDTune, double xiXTune, double xiYTune, 
                         double xiXb, double xiYb, double xiXMu, double xiYMu, 
                         List xData_cov, List xPredData_cov, //NumericMatrix lambds_run, 
                         NumericMatrix lambdaInit, NumericVector lambdaTune, NumericVector lambdaMu, NumericVector lambdaSig, List yData_ind,
                         List xBasis_ind, List yBasis_ind, List xBasisPred_ind
){
  int nToSave = nIter/nThin;
  NumericMatrix Bx1 = Bx[1];
  int m = Bx1.ncol();
  int nData = coordsData.nrow();
  int nPred = coordsPred.nrow();
  // NumericVector q = sapply(By,nrow)
  // int p = Bx.nrow();
  int r = ByPred.nrow();
  arma::vec exes = xData[1];
  int n_times = exes.size();
  int nlam =4;
  
  
  arma::mat coordsDataArma(coordsData.begin(),coordsData.nrow(),coordsData.ncol(),false);
  arma::mat coordsPredArma(coordsPred.begin(),coordsPred.nrow(),coordsPred.ncol(),false);
  arma::mat coordsAllArma = arma::join_cols(coordsPredArma,coordsDataArma);
  arma::mat distsAllArma = calcDistsArma(coordsAllArma);
  
  arma::mat H0AlphaAll = H0functionArma(phiAlpha,distsAllArma);
  arma::mat H0BetaAll = H0functionArma(phiBeta,distsAllArma);
  arma::mat H0Alpha11 = H0AlphaAll.submat(0,0,(nPred-1),(nPred-1));
  arma::mat H0Alpha12 = H0AlphaAll.submat(0,nPred,(nPred-1),(nData+nPred-1));
  arma::mat H0Alpha21 = H0AlphaAll.submat(nPred,0,(nData+nPred-1),(nPred-1));
  arma::mat H0Alpha22 = H0AlphaAll.submat(nPred,nPred,(nData+nPred-1),(nData+nPred-1));
  arma::mat H0Alpha22Inv = arma::inv_sympd(H0Alpha22);
  arma::mat H0Beta11 = H0BetaAll.submat(0,0,(nPred-1),(nPred-1));
  arma::mat H0Beta12 = H0BetaAll.submat(0,nPred,(nPred-1),(nData+nPred-1));
  arma::mat H0Beta21 = H0BetaAll.submat(nPred,0,(nData+nPred-1),(nPred-1));
  arma::mat H0Beta22 = H0BetaAll.submat(nPred,nPred,(nData+nPred-1),(nData+nPred-1));
  arma::mat H0Beta22Inv = arma::inv_sympd(H0Beta22);
  
  arma::mat SigmaDArma(SigmaD.begin(),SigmaD.nrow(),SigmaD.ncol(),false);
  arma::mat SigmaDArmaInv = arma::inv_sympd(SigmaDArma);
  arma::vec muDArma(muD.begin(),muD.size());
  
  arma::mat ByPredArma(ByPred.begin(),ByPred.nrow(),ByPred.ncol(),false);
  arma::mat SigmaAlphaPredA = H0Alpha11-H0Alpha12*H0Alpha22Inv*H0Alpha21;
  arma::mat SigmaBetaPredA = H0Beta11-H0Beta12*H0Beta22Inv*H0Beta21;
  arma::mat SigmaAlphaPredB = H0Alpha12*H0Alpha22Inv;
  arma::mat SigmaBetaPredB = H0Beta12*H0Beta22Inv;
  
  arma::mat CholSigmaAlphaPredA = choleskyEfficient(SigmaAlphaPredA);
  arma::mat CholSigmaBetaPredA = choleskyEfficient(SigmaBetaPredA);
  
  int nVar2 = 3 + 4 * m * nData + 4 * m * nPred + r * nPred + 2 + 2 + 4 * nData + 4 * nPred;
  
  arma::mat outputArma(nToSave+1,nVar2);
  outputArma(0,0) = sigmaAlphaPrecInit;
  outputArma(0,1) = sigmaBetaPrecInit;
  outputArma(0,2) = sigmaCPrecInit;
  arma::mat alphaInitArma(alphaInit.begin(),alphaInit.nrow(),alphaInit.ncol(),false);
  arma::mat betaInitArma(betaInit.begin(),betaInit.nrow(),betaInit.ncol(),false);
  arma::mat cInitArma(cInit.begin(),cInit.nrow(),cInit.ncol(),false);
  arma::mat dInitArma(dInit.begin(),dInit.nrow(),dInit.ncol(),false);
  // Additions
  arma::mat lambdaInitArma(lambdaInit.begin(),lambdaInit.nrow(),lambdaInit.ncol(),false);
  double xiXs = xiXInit;
  double xiYs = xiYInit;
  arma::vec lambdaTunes(lambdaTune.begin(), lambdaTune.size(), false);
  arma::vec lambdaSigs(lambdaSig.begin(), lambdaSig.size(), false);
  arma::vec lambdaMus(lambdaMu.begin(), lambdaMu.size(), false);
  // arma::mat lambdaInitArma(nData, nlam);
  // for(int k = 0; k < nData; ++k){
  //   lambdaInitArma.row(k) = lambdsArma(0,arma::span(k*nlam,(k+1)*nlam-1));
  // }
  arma::vec alphaInitArmaVec = arma::vectorise(alphaInitArma);
  arma::vec betaInitArmaVec = arma::vectorise(betaInitArma);
  arma::vec cInitArmaVec = arma::vectorise(cInitArma);
  arma::vec dInitArmaVec = arma::vectorise(dInitArma);
  // Alterations
  arma::vec lambdaInitArmaVec = arma::vectorise(lambdaInitArma);
  arma::vec cPredInitArmaVec = arma::vectorise(cInitArma);
  arma::vec dPredInitArmaVec = arma::vectorise(dInitArma);
  
  outputArma(0,arma::span(3,2+m*nData)) = alphaInitArmaVec.t();
  outputArma(0,arma::span(3+m*nData,2+2*m*nData)) = betaInitArmaVec.t();
  outputArma(0,arma::span(3+2*m*nData,2+3*m*nData)) = cInitArmaVec.t();
  // outputArma(0,3*m*nData+4) = sigmaXPrecInit;
  outputArma(0,arma::span(3+3*m*nData,2+4*m*nData)) = dInitArmaVec.t();
  // Alterations
  outputArma(0,3+4*m*nData) = xiXs;
  outputArma(0,4+4*m*nData) = xiYs;
  outputArma(0,arma::span(5+4*m*nData,4+4*m*nData+4*nData)) = lambdaInitArmaVec.t();
  outputArma(0,arma::span(5+4*m*nData+2*m*nPred + 4*nData,4+4*m*nData+3*m*nPred+4*nData)) = arma::zeros(m*nPred).t();//dInitArma.col(1).t();// //dpred
  outputArma(0,arma::span(5+4*m*nData+3*m*nPred + 4*nData,4+4*m*nData+4*m*nPred+4*nData)) = arma::zeros(m*nPred).t();//cInitArma.col(1).t();//cpred
  
  // NumericMatrix output = wrap(outputArma);
  // return output;
  int l_count = 0;
  for(int i=1; i<(nToSave+1); i++){
    arma::mat subchainOutputArma(nThin+1,nVar2);
    subchainOutputArma.row(0) = outputArma.row(i-1);
    for(int j=1; j<(nThin+1); j++){
      arma::mat alphaPrev(m,nData);
      arma::mat betaPrev(m,nData);
      arma::mat cPrev(m,nData);
      arma::mat dPrev(m,nData);
      // Alterations
      arma::mat lambdaPrev(nData,nlam);
      arma::mat dPredPrev(m,nPred);
      arma::mat cPredPrev(m,nPred);
      //
      // // Rcout << "Did we reach here? Yes! F : " << nData << "\n";
      for(int k=0; k<nData; k++){
        alphaPrev.col(k) = subchainOutputArma(j-1,arma::span(3+k*m,2+(k+1)*m)).t();
        betaPrev.col(k) = subchainOutputArma(j-1,arma::span(3+nData*m+k*m,2+nData*m+(k+1)*m)).t();
        cPrev.col(k) = subchainOutputArma(j-1,arma::span(3+2*nData*m+k*m,2+2*nData*m+(k+1)*m)).t();
        dPrev.col(k) = subchainOutputArma(j-1,arma::span(3+3*nData*m+k*m,2+3*nData*m+(k+1)*m)).t();
        lambdaPrev.row(k) = subchainOutputArma(j-1,arma::span(5+4*m*nData+k*nlam,4+4*nData*m+(k+1)*nlam));
      }
      // Alterations
      double xiXPrev = subchainOutputArma(j-1,3+4*m*nData);
      double xiYPrev = subchainOutputArma(j-1,4+4*m*nData);
      //
      // // Rcout << "before : " << xiYPrev << "\n";
      for(int k=0; k<nPred; k++){
        dPredPrev.col(k) = subchainOutputArma(j-1,arma::span(5+4*m*nData+4*nData + 2*m*nPred + k*m,4+4*m*nData+4*nData + 2*m*nPred + (k+1)*m)).t();
        cPredPrev.col(k) = subchainOutputArma(j-1,arma::span(5+4*m*nData+4*nData + 3*m*nPred + k*m,4+4*m*nData+4*nData + 3*m*nPred + (k+1)*m)).t();
      }
      // // Rcout << "after : " << xiYPrev << "\n";
      /////////////// Gibb's sampling ///////////////
      // Rcout << "Did we reach here? Yes! F 2 : " << nData << "\n";
      // Update precision of Alpha (sigmaAlphaPrecUpdate)
      double aAlphaPrecUpdate = aAlpha+m*nData/2;
      double bAlphaPrecUpdate = bAlpha+0.5*arma::trace(H0Alpha22Inv*alphaPrev.t()*alphaPrev);
      double sigmaAlphaPrecUpdate = R::rgamma(aAlphaPrecUpdate,1/bAlphaPrecUpdate);
      // Update precision of Beta (sigmaBetaPrecUpdate)
      double aBetaPrecUpdate = aBeta+m*nData/2;
      double bBetaPrecUpdate = bBeta+0.5*arma::trace(H0Beta22Inv*(betaPrev-arma::ones(m,nData)).t()*(betaPrev-arma::ones(m,nData)));
      double sigmaBetaPrecUpdate = R::rgamma(aBetaPrecUpdate,1/bBetaPrecUpdate);
      
      // Update precision of C (sigmaCPrecUpdate)
      double aCPrecUpdate = aC+0.5*m*nData;
      arma::mat dMABG = cPrev-(alphaPrev+betaPrev%dPrev);
      double bCPrecUpdate = bC+0.5*arma::trace(dMABG.t()*dMABG);
      double sigmaCPrecUpdate = R::rgamma(aCPrecUpdate,1/bCPrecUpdate);
      // Rcout << "Did we reach here? Yes! F 3 : " << nData << "\n";
      /////////////// Metropolis Hastings ///////////////
      // Initialise variables
      arma::mat alphaUpdate(m,nData);
      arma::mat betaUpdate(m,nData);
      arma::mat cUpdate(m,nData);
      arma::mat dUpdate(m,nData);
      arma::mat lambdaUpdate(nData, nlam);
      arma::mat ps(nData,n_times);
      arma::mat alpha_var = (1/sigmaAlphaPrecUpdate) * arma::inv_sympd(H0Alpha22Inv);
      arma::mat beta_var = (1/sigmaBetaPrecUpdate) * arma::inv_sympd(H0Beta22Inv);
      
      // Update Alpha
      for(int k=0; k<m; k++){
        // Define old values
        arma::vec alphaPrevJ = alphaPrev.row(k).t();
        arma::vec betaPrevJ = betaPrev.row(k).t();
        arma::vec cPrevJ = cPrev.row(k).t();
        arma::vec dPrevJ = dPrev.row(k).t();
        arma::vec muCold = alphaPrevJ + betaPrevJ%dPrevJ;
        // Rcout << "Did we reach here? Yes! F 4 : " << nData << "\n";
        // Define new values
        arma::vec newAlpha(nData);
        for(int l=0;l<nData;l++){
          newAlpha(l) = alphaPrevJ(l) + R::rnorm(0,sigmaATune);
        }
        arma::vec muCnew = newAlpha + betaPrevJ%dPrevJ;
        
        // Calculate Alpha Priors
        double oldAlphaPrior = mvnPDF(alphaPrevJ, arma::zeros(nData), alpha_var);
        double newAlphaPrior = mvnPDF(newAlpha, arma::zeros(nData), alpha_var);
        
        // Calculate likelihood of C
        arma::vec coldlik(nData);
        arma::vec cnewlik(nData);
        for(int l=0;l<nData;l++){
          coldlik(l) = log(normalPDF(cPrevJ(l), muCold(l), 1/sigmaCPrecUpdate));
          cnewlik(l) = log(normalPDF(cPrevJ(l), muCnew(l), 1/sigmaCPrecUpdate));
        }
        
        // Calculate Acceptance ratio
        double pNew = sum(cnewlik) + newAlphaPrior;
        double pOld = sum(coldlik) + oldAlphaPrior;
        double paccept = pNew - pOld;
        
        // Update alpha if accepted
        if (std::isnan(paccept)) {
          paccept = -1000.0;
        }
        if (paccept >= std::log(R::runif(0.0, 1.0))) {
          alphaUpdate.row(k) = newAlpha.t();
        } else {
          alphaUpdate.row(k) = alphaPrevJ.t();
        }
      }
      
      // Update Beta
      for(int k=0; k<m; k++){
        // Define old values
        arma::vec alphaPrevJ = alphaUpdate.row(k).t();
        arma::vec betaPrevJ = betaPrev.row(k).t();
        arma::vec cPrevJ = cPrev.row(k).t();
        arma::vec dPrevJ = dPrev.row(k).t();
        arma::vec muCold = alphaPrevJ + betaPrevJ%dPrevJ;
        
        // Define new values
        arma::vec newBeta(nData);
        for(int l=0;l<nData;l++){
          newBeta(l) = betaPrevJ(l) + R::rnorm(0,sigmaBTune);
        }
        arma::vec muCnew = alphaPrevJ + newBeta%dPrevJ;
        
        // Calculate Beta Priors
        // arma::vec onesVector(nData, arma::fill::ones);
        double oldBetaPrior = mvnPDF(betaPrevJ, arma::ones(nData), beta_var);
        double newBetaPrior = mvnPDF(newBeta, arma::ones(nData), beta_var);
        
        // Calculate likelihood of C
        arma::vec coldlik(nData);
        arma::vec cnewlik(nData);
        for(int l=0;l<nData;l++){
          coldlik(l) = log(normalPDF(cPrevJ(l), muCold(l), 1/sigmaCPrecUpdate));
          cnewlik(l) = log(normalPDF(cPrevJ(l), muCnew(l), 1/sigmaCPrecUpdate));
        }
        
        // Calculate Acceptance ratio
        double pNew = sum(cnewlik) + newBetaPrior;
        double pOld = sum(coldlik) + oldBetaPrior;
        double paccept = pNew - pOld;
        
        // Update alpha if accepted
        if (std::isnan(paccept)) {
          paccept = -1000.0;
        }
        if (paccept >= std::log(R::runif(0.0, 1.0))) {
          betaUpdate.row(k) = newBeta.t();
        } else {
          betaUpdate.row(k) = betaPrevJ.t();
        }
      }
      
      // Update ps
      for(int k=0; k<nData; k++){
        //lambdaUpdate.row(k) = lambdsArma(l_count+1,arma::span(k*nlam,(k+1)*nlam-1));
        NumericMatrix xcovs = xData_cov[k];
        arma::mat xcovKJ(xcovs.begin(),xcovs.nrow(),xcovs.ncol(),false);
        NumericVector yData_indJ = yData_ind[k];
        arma::vec yData_indArma(yData_indJ.begin(),yData_indJ.size(),false);
        // arma::mat ByJArma(ByJ.begin(),ByJ.nrow(),ByJ.ncol(),false);
        
        int xcov_rows = xcovKJ.n_rows;
        
        // Old values
        arma::vec lambdaPrevJ = lambdaPrev.row(k).t();
        arma::vec oldp(xcov_rows);
        
        for(int l=0;l<xcov_rows;l++){
          arma::vec tmp1 = (xcovKJ.row(l)*lambdaPrevJ);
          double tmp2 = sumWithoutNaN(tmp1);
          oldp(l) = invlogit(tmp2);
        }
        
        // New values
        arma::vec lambdaNew(nlam);
        for(int l=0; l<nlam; l++){
          lambdaNew(l) = lambdaPrevJ(l) + sampleFromNormal(1,0.0,lambdaTunes(l))(0);
        }
        arma::vec newp(xcov_rows);
        // arma::vec tmp(xcov_rows);
        for(int l=0;l<xcov_rows;l++){
          arma::vec tmp1 = (xcovKJ.row(l)*lambdaNew);
          double tmp2 = sumWithoutNaN(tmp1);
          newp(l) = invlogit(tmp2);
        }
        
        // Compute Priors
        arma::mat LamD = arma::diagmat(lambdaSigs);
        double oldPrior = mvnPDF(lambdaPrevJ, lambdaMus, LamD);
        double newPrior = mvnPDF(lambdaNew, lambdaMus, LamD);
        
        // Compute Likelihood
        arma::vec oldLik(xcov_rows);
        arma::vec newLik(xcov_rows);
        
        for(int l=0; l<xcov_rows; l++){
          if (!std::isnan(yData_indArma(l))) { 
            oldLik(l) = log(binomPDF(yData_indArma(l),oldp(l)));
            newLik(l) = log(binomPDF(yData_indArma(l),newp(l)));
          }
          oldLik(l) = NA_REAL;
          newLik(l) = NA_REAL;
        }
        
        // Accept or reject
        double pOld = sumWithoutNaN(oldLik) + oldPrior;
        double pNew = sumWithoutNaN(newLik) + newPrior;
        double paccept = pNew - pOld;
        
        if (paccept >= std::log(R::runif(0.0, 1.0))) {
          lambdaUpdate.row(k) = lambdaNew.t();
          ps.row(k) = newp.t();
        } else {
          lambdaUpdate.row(k) = lambdaPrevJ.t();
          ps.row(k) = oldp.t();
        }
        // lambdaUpdate.row(k) = lambdaPrevJ.t();
        // ps.row(k) = oldp.t();
      }
      
      // Update Xi (both unique=T and unique=F)
      double xiXUpdate;
      double xiYUpdate;
      if(uniqueXi==1){
        // Initialise variables
        arma::vec pxiXold_v(nData);
        arma::vec pxiXnew_v(nData);
        
        // Define old values
        double oldXi = xiXPrev;
        double newXi = oldXi + R::rnorm(0,xiXTune);
        while (newXi > 0.5 || newXi < -0.5){
          newXi = oldXi + R::rnorm(0,xiXTune);
        }
        newXi = roundToSig(newXi,5);
        
        for(int k=0; k<nData; k++){
          // Define previous parameter values
          arma::vec cPrevJ = cPrev.col(k).t();
          arma::vec dPrevJ = dPrev.col(k).t();
          NumericVector xJ = xData[k];
          arma::vec xJArma(xJ.begin(),xJ.size());
          NumericMatrix BxJ = Bx[k];
          arma::mat BxJArma(BxJ.begin(),BxJ.nrow(),BxJ.ncol(),false);
          NumericVector yJ = yData[k];
          arma::vec yJArma(yJ.begin(),yJ.size());
          NumericMatrix ByJ = By[k];
          arma::mat ByJArma(ByJ.begin(),ByJ.nrow(),ByJ.ncol(),false);
          
          // 
          int xn = xJ.size();
          arma::vec px(xn);
          for(int l=0; l<xn; l++){
            if(xJArma(l) > 0){
              px(l) = 1;
            } else {
              px(l) = 0;
            }
          }
          arma::vec py = ps.row(k).t();
          
          // Calculate scale parameters
          arma::vec sigxiX = arma::exp(BxJArma * dPrevJ);
          arma::vec sigxiY = arma::exp(ByJArma * cPrevJ);
          
          // Obtain likelihood
          arma::vec tmp_xiXlikOld(xn);
          for(int l=0; l<xn; l++){
            tmp_xiXlikOld(l) = diracgpdPDF(xJArma(l), oldXi, sigxiX(l), px(l));
          }
          double xiX_likOld = sum(log(tmp_xiXlikOld));
          
          int yn = yJ.size();
          arma::vec tmp_xiYlikOld(yn);
          for(int l=0; l<yn; l++){
            tmp_xiYlikOld(l) = diracgpdPDF(yJArma(l), oldXi, sigxiY(l), py(l));
          }
          double xiY_likOld = sum(log(tmp_xiYlikOld));
          
          arma::vec tmp_xiXlikNew(xn);
          for(int l=0; l<xn; l++){
            tmp_xiXlikNew(l) = diracgpdPDF(xJArma(l), newXi, sigxiX(l), px(l));
          }
          double xiX_likNew = sum(log(tmp_xiXlikNew));
          
          arma::vec tmp_xiYlikNew(yn);
          for(int l=0; l<yn; l++){
            tmp_xiYlikNew(l) = diracgpdPDF(yJArma(l), newXi, sigxiY(l), py(l));
          }
          double xiY_likNew = sum(log(tmp_xiYlikNew));
          
          // Compute acceptance probabilities
          pxiXold_v(k) = xiX_likOld + xiY_likOld;
          pxiXnew_v(k) = xiX_likNew + xiY_likNew ;
        }
        // Xi Priors
        double xiXPriorOld = laplace_trans_pdf(oldXi, xiXb, xiXMu);
        double xiXPriorNew = laplace_trans_pdf(newXi, xiXb, xiXMu);
        double ncOld = normalPDF(oldXi, 0, xiXTune);
        double ncNew = normalPDF(newXi, 0, xiXTune);
        double pxiold = sum(pxiXold_v) + xiXPriorOld + log(ncOld);
        double pxinew = sum(pxiXnew_v) + xiXPriorNew + log(ncNew);
        double pxiaccept = pxinew - pxiold ;
        
        if (pxiaccept >= std::log(R::runif(0.0, 1.0))) {
          xiXUpdate = newXi;
        } else {
          xiXUpdate = oldXi;
        }
        
      } else {
        // Initialise variables
        arma::vec pxiXold_v(nData);
        arma::vec pxiXnew_v(nData);
        arma::vec pxiYold_v(nData);
        arma::vec pxiYnew_v(nData);
        
        // Update XiX
        // Define old values
        double oldXiX = xiXPrev;
        double newXiX = oldXiX + R::rnorm(0,xiXTune);
        while (newXiX > 0.5 || newXiX < -0.5){
          newXiX = oldXiX + R::rnorm(0,xiXTune);
        }
        // if(std::abs(newXiX)<1e-10){
        //   // Rcout << "newXiX small is : " << newXiX << "\n";
        //   newXiX = 0.0;
        // }
        // newXiX = roundToSig(newXiX,5);
        
        for(int k=0; k<nData; k++){
          // Define previous parameter values
          arma::vec dPrevJ = dPrev.col(k);
          NumericVector xJ = xData[k];
          arma::vec xJArma(xJ.begin(),xJ.size());
          NumericMatrix BxJ = Bx[k];
          arma::mat BxJArma(BxJ.begin(),BxJ.nrow(),BxJ.ncol(),false);
          
          int xn = xJ.size();
          arma::vec px(xn);
          for(int l=0; l<xn; l++){
            if(xJArma(l) > 0){
              px(l) = 1.0;
            } else {
              px(l) = 0.0;
            }
          }
          
          // Calculate scale parameters
          arma::vec sigxiX = arma::exp(BxJArma * dPrevJ);
          
          // Obtain likelihood
          arma::vec tmp_xiXlikOld(xn);
          for(int l=0; l<xn; l++){
            tmp_xiXlikOld(l) = diracgpdPDF(xJArma(l), oldXiX, sigxiX(l), px(l));
          }
          double xiX_likOld = sumWithoutNaN(log(tmp_xiXlikOld));
          
          arma::vec tmp_xiXlikNew(xn);
          for(int l=0; l<xn; l++){
            tmp_xiXlikNew(l) = diracgpdPDF(xJArma(l), newXiX, sigxiX(l), px(l));
          }
          double xiX_likNew = sumWithoutNaN(log(tmp_xiXlikNew));
          
          // Compute acceptance probabilities
          pxiXold_v(k) = xiX_likOld;
          pxiXnew_v(k) = xiX_likNew;
        }
        // Rcout << "Did we reach here? Yes! D : " << nData << "\n";
        // XiX Priors
        double xiXPriorOld = laplace_trans_pdf(oldXiX, xiXb, xiXMu);
        double xiXPriorNew = laplace_trans_pdf(newXiX, xiXb, xiXMu);
        double ncOldX = R::pnorm(oldXiX,0,xiXTune,TRUE,TRUE);  //normalPDF(oldXiX, 0, xiXTune);
        double ncNewX = R::pnorm(newXiX,0,xiXTune,TRUE,TRUE); //normalPDF(newXiX, 0, xiXTune);
        // Calculate acceptance ratios
        double pxiXold = sum(pxiXold_v) + xiXPriorOld + ncNewX;
        double pxiXnew = sum(pxiXnew_v) + xiXPriorNew + ncOldX;
        double pxiXaccept = pxiXnew - pxiXold;
        if (pxiXaccept >= std::log(R::runif(0.0, 1.0))) {
          xiXUpdate = newXiX;
        } else {
          xiXUpdate = oldXiX;
        }
        
        // Update XiY
        // Define old values
        double oldXiY = xiYPrev;
        double newXiY = oldXiY + R::rnorm(0,xiYTune);
        while (newXiY > 0.5 || newXiY < -0.5){
          newXiY = oldXiY + R::rnorm(0,xiYTune);
        }
        
        for(int k=0; k<nData; k++){
          // Define previous parameter values
          arma::vec cPrevJ = cPrev.col(k);
          NumericVector yJ = yData[k];
          arma::vec yJArma(yJ.begin(),yJ.size());
          NumericMatrix ByJ = By[k];
          arma::mat ByJArma(ByJ.begin(),ByJ.nrow(),ByJ.ncol(),false);
          NumericVector yData_indJ = yData_ind[k];
          arma::vec yData_indArma(yData_indJ.begin(),yData_indJ.size(),false);
          int yn = yJ.size();
          
          arma::vec psJ = ps.row(k).t();
          arma::vec py(n_times);
          int count = 0;
          for (int l=0; l < n_times; ++l) {
            if (!std::isnan(yData_indJ(l))) {
              py(count++) = psJ(l);
            }
          }
          py.resize(count);
          
          // Calculate scale parameters
          arma::vec sigxiY = arma::exp(ByJArma * cPrevJ);
          // Rcout << "sigxiY : " << sigxiY.n_elem << "\n";
          // // Rcout << "py length : " << py.n_elem << "\n";
          
          // Obtain likelihood
          arma::vec tmp_xiYlikOld(yn);
          for(int l=0; l<yn; l++){
            tmp_xiYlikOld(l) = diracgpdPDF(yJArma(l), oldXiY, sigxiY(l), py(l));
          }
          double xiY_likOld = sumWithoutNaN(log(tmp_xiYlikOld));
          
          arma::vec tmp_xiYlikNew(yn);
          for(int l=0; l<yn; l++){
            tmp_xiYlikNew(l) = diracgpdPDF(yJArma(l), newXiY, sigxiY(l), py(l));
          }
          double xiY_likNew = sumWithoutNaN(log(tmp_xiYlikNew));
          
          // Compute acceptance probabilities
          pxiYold_v(k) = xiY_likOld;
          pxiYnew_v(k) = xiY_likNew;
        } 
        
        // XiY Priors
        double xiYPriorOld = laplace_trans_pdf(oldXiY, xiYb, xiYMu);
        double xiYPriorNew = laplace_trans_pdf(newXiY, xiYb, xiYMu);
        double ncOldY = R::pnorm(oldXiY,0,xiYTune,TRUE,TRUE);  //normalPDF(oldXiX, 0, xiXTune);
        double ncNewY = R::pnorm(newXiY,0,xiYTune,TRUE,TRUE); //normalPDF(newXiX, 0, xiXTune);
        double pxiYold = sum(pxiYold_v) + xiYPriorOld + ncNewY;
        double pxiYnew = sum(pxiYnew_v) + xiYPriorNew + ncOldY;
        double pxiYaccept = pxiYnew - pxiYold;
        // // Rcout << "sum(pYiYold_v) : " << sum(pxiYold_v)<< "\n";
        // // Rcout << "sum(pxiYnew_v) : " << sum(pxiYnew_v) << "\n";
        // // Rcout << "pxiYnew: " << pxiYnew<< "\n";
        // // Rcout << "pxiYold : " << pxiYold << "\n";
        // // Rcout << "xiYPriorOld : " << xiYPriorOld << "\n";
        // // Rcout << "xiYPriorNew : " << xiYPriorNew << "\n";
        // // Rcout << "log(ncOld) : " << log(ncOld) << "\n";
        // // Rcout << "log(ncNew): " << log(ncNew)<< "\n";
        if (pxiYaccept >= std::log(R::runif(0.0, 1.0))) {
          xiYUpdate = newXiY;
        } else {
          xiYUpdate = oldXiY; 
        }
        
      }
      // Rcout << "Did we reach here? Yes! C : " << xiYPrev << "\n";
      // Update D (dUpdate)
      for (int k = 0; k < nData; ++k) {
        NumericVector xJ = xData[k];
        arma::vec xJArma(xJ.begin(),xJ.size());
        NumericMatrix BxJ = Bx[k];
        arma::mat BxJArma(BxJ.begin(),BxJ.nrow(),BxJ.ncol(),false);
        arma::vec alphaJ = alphaUpdate.col(k);
        arma::vec betaJ = betaUpdate.col(k);
        arma::vec cJ = cPrev.col(k);
        List xBasisind = xBasis_ind[k];
        // List xBasisPredind = xBasisPred_ind[k];
        
        int xn = xJ.size();
        arma::vec px(xn);
        for(int l=0; l<xn; l++){
          if(xJArma(l) > 0){
            px(l) = 1;
          } else {
            px(l) = 0;
          }
        }
        
        arma::vec dJold = dPrev.col(k);
        arma::vec dJnew = dPrev.col(k) + sampleFromNormal(m,0.0,sigmaDTune);
        double xiXJ = xiXUpdate;
        
        // Compute scale parameters
        arma::vec dxsigmaold = arma::exp(BxJArma * dJold);
        arma::vec dxsigmanew = arma::exp(BxJArma * dJnew);
        
        // Example of computing cMuOld and cMuNew directly in Armadillo
        arma::vec cMuOld = alphaJ + betaJ%dJold;
        arma::vec cMuNew = alphaJ + betaJ%dJnew;
        
        // Obtain the prior of c for the old and the new values
        arma::vec cPrioriOldd(m);
        arma::vec cPrioriNewd(m);
        for(int l=0;l<m;l++){
          cPrioriOldd(l) = log(normalPDF(cJ(l), cMuOld(l), 1/sigmaCPrecUpdate));
          cPrioriNewd(l) = log(normalPDF(cJ(l), cMuNew(l), 1/sigmaCPrecUpdate));
        }
        
        // Obtain the prior of d for the old and the new values
        double dPrioriOldd = mvnPDF(dJold, muDArma, SigmaDArma);
        double dPrioriNewd = mvnPDF(dJnew, muDArma, SigmaDArma);
        
        // Obtain likelihood
        arma::vec xlikOld(m);
        arma::vec xlikNew(m);
        for(int l=0;l<m;l++){
          // Define the vectors for iteration
          NumericVector xBasisindVec = xBasisind[l];
          arma::vec xBasisArma(xBasisindVec.begin(),xBasisindVec.size());
          arma::vec tmp(xBasisArma.n_elem);
          arma::vec tmp2(xBasisArma.n_elem);
          for(int t=0; t<xBasisArma.n_elem;t++){
            int tmp_int = xBasisArma(t)-1;
            tmp(t) = diracgpdPDF(xJArma(tmp_int), xiXJ, dxsigmaold(tmp_int), px(tmp_int));
            tmp2(t) = diracgpdPDF(xJArma(tmp_int), xiXJ, dxsigmanew(tmp_int), px(tmp_int));
          }
          xlikOld(l) = sum(log(tmp));
          xlikNew(l) = sum(log(tmp2));
        }
        
        // Acceptance
        arma::vec pDold = dPrioriOldd + xlikOld + cPrioriOldd;
        arma::vec pDnew = dPrioriNewd + xlikNew + cPrioriNewd;
        arma::vec pDaccept = pDnew - pDold;
        for(int l=0;l<m;l++){
          if (pDaccept(l) >= std::log(R::runif(0.0, 1.0))) {
            dUpdate(l,k) = dJnew(l);
          } else {
            dUpdate(l,k) = dJold(l);
          }
        }
      }
      
      // Rcout << "Did we reach here? Yes! B : " << xiYPrev << "\n";
      // Update C (cUpdate)
      for(int k=0; k<nData; k++){
        // Define old values
        NumericVector yJ = yData[k];
        arma::vec yJArma(yJ.begin(),yJ.size());
        NumericMatrix ByJ = By[k];
        arma::mat ByJArma(ByJ.begin(),ByJ.nrow(),ByJ.ncol(),false);
        arma::vec alphaUpdateJ = alphaUpdate.col(k);
        arma::vec betaUpdateJ = betaUpdate.col(k);
        arma::vec dUpdateJ = dUpdate.col(k);
        NumericVector yData_indJ = yData_ind[k];
        arma::vec yData_indArma(yData_indJ.begin(),yData_indJ.size(),false);
        
        List yBasisind = yBasis_ind[k];
        int yn = yJ.size();
        arma::vec psJ = ps.row(k).t();
        arma::vec py(n_times);
        int count = 0;
        for (int l=0; l < n_times; ++l) {
          if (!std::isnan(yData_indJ(l))) {
            py(count++) = psJ(l);
          }
        }
        py.resize(count);
        // 
        // Initialize parameters
        arma::vec cPrioriOld(m);
        arma::vec cPrioriNew(m);
        arma::vec ylikOldj(yJ.size());
        arma::vec ylikNewj(yJ.size());
        // arma::vec py = ps.row(k).t();
        
        // Obtain old C and propose new C
        arma::vec cJold = cPrev.col(k);
        arma::vec cJnew = cPrev.col(k) + sampleFromNormal(m,0.0,sigmaCTune);
        
        // Define parameters
        arma::vec muTemp = alphaUpdateJ + betaUpdateJ % dUpdateJ;
        double xiYJ;
        if(uniqueXi==1){
          xiYJ = xiXUpdate;
        } else {
          xiYJ = xiYUpdate;
        }
        arma::vec cysigmaold = arma::exp(ByJArma * cJold);
        arma::vec cysigmanew = arma::exp(ByJArma * cJnew);
        
        // Obtain the prior of c for the old and the new values
        arma::vec cPrioriOldc(m);
        arma::vec cPrioriNewc(m);
        for(int l=0;l<m;l++){
          cPrioriOldc(l) = log(normalPDF(cJold(l), muTemp(l), 1/sigmaCPrecUpdate));
          cPrioriNewc(l) = log(normalPDF(cJnew(l), muTemp(l), 1/sigmaCPrecUpdate));
        }
        // Obtain likelihood
        arma::vec ylikOld(m);
        arma::vec ylikNew(m);
        for(int l=0;l<m;l++){
          // Define the vectors for iteration
          NumericVector yBasisindVec = yBasisind[l];
          arma::vec yBasisArma(yBasisindVec.begin(),yBasisindVec.size());
          arma::vec tmp(yBasisArma.n_elem);
          arma::vec tmp2(yBasisArma.n_elem);
          for(int t=0; t<yBasisArma.n_elem;t++){
            int tmp_int = yBasisArma(t)-1;
            tmp(t) = diracgpdPDF(yJArma(tmp_int), xiYJ, cysigmaold(tmp_int), py(tmp_int));
            tmp2(t) = diracgpdPDF(yJArma(tmp_int), xiYJ, cysigmanew(tmp_int), py(tmp_int));
            
          }
          
          ylikOld(l) = sum(log(tmp));
          ylikNew(l) = sum(log(tmp2));
        }
        
        // Calculate acceptance ratios
        arma::vec pCold = ylikOld + cPrioriOldc;
        arma::vec pCnew = ylikNew + cPrioriNewc;
        arma::vec pCaccept = pCnew - pCold;
        
        for(int l=0;l<m;l++){
          if(pCaccept(l) >= std::log(R::runif(0.0, 1.0))) {
            cUpdate(l,k) = cJnew(l);
          } else {
            cUpdate(l,k) = cJold(l);
          }
        }
      }
      // Rcout << "Did we reach here? Yes! A : " << xiYPrev << "\n";
      // Store updates
      subchainOutputArma(j,0) = sigmaAlphaPrecUpdate;
      subchainOutputArma(j,1) = sigmaBetaPrecUpdate;
      subchainOutputArma(j,2) = sigmaCPrecUpdate;
      arma::vec alphaUpdateVec = arma::vectorise(alphaUpdate);
      arma::vec betaUpdateVec = arma::vectorise(betaUpdate);
      arma::vec cUpdateVec = arma::vectorise(cUpdate);
      arma::vec dUpdateVec = arma::vectorise(dUpdate);
      //
      arma::vec lambdaUpdateVec = arma::vectorise(lambdaUpdate.t());
      //
      subchainOutputArma(j,arma::span(3,2+m*nData)) = alphaUpdateVec.t();
      subchainOutputArma(j,arma::span(3+m*nData,2+2*m*nData)) = betaUpdateVec.t();
      subchainOutputArma(j,arma::span(3+2*m*nData,2+3*m*nData)) = cUpdateVec.t();
      subchainOutputArma(j,arma::span(3+3*m*nData,2+4*m*nData)) = dUpdateVec.t();
      
      // Adapt and survive, little bear
      subchainOutputArma(j,3+4*m*nData) = xiXUpdate;
      subchainOutputArma(j,4+4*m*nData) = xiYUpdate;
      subchainOutputArma(j,arma::span(5+4*m*nData,4+4*m*nData+4*nData)) = lambdaUpdateVec.t();
      //
      arma::mat alphaPredUpdate(m,nPred);
      arma::mat betaPredUpdate(m,nPred);
      arma::mat dPredUpdate(m,nPred);
      arma::mat cPredUpdate(m,nPred);
      
      //
      arma::mat lambdaPredUpdate(nPred,nlam);
      arma::mat p_PredUpdate(nPred, r);
      // arma::mat dPredOld(m,nPred);
      // arma::mat cPredOld(m,nPred);
      //
      
      // Update Alpha and Beta Predict
      for(int k=0; k<m; k++){
        arma::vec alphaJ = alphaUpdate.row(k).t();
        arma::vec muAlphaPredJ = SigmaAlphaPredB*alphaJ;
        // arma::mat SigmaAlphaPredJ = sigmaAlphaPrecUpdate*SigmaAlphaPredA;
        alphaPredUpdate.row(k) = muAlphaPredJ.t() + arma::randn(1,nPred)*sqrt(1/sigmaAlphaPrecUpdate)*CholSigmaAlphaPredA;
        
        arma::vec betaJ = betaUpdate.row(k).t();
        arma::vec muBetaPredJ = arma::ones(nPred) + SigmaBetaPredB*(betaJ-arma::ones(nData));
        // arma::mat SigmaBetaPredJ = sigmaBetaPrecUpdate*SigmaBetaPredA;
        betaPredUpdate.row(k) = muBetaPredJ.t() + arma::randn(1,nPred)*sqrt(1/sigmaBetaPrecUpdate)*CholSigmaBetaPredA;
      }
      // Rcout << "Did we reach here? Yes! : " << xiYPrev << "\n";
      // Update D Predict
      for(int k=0; k<nPred; k++){
        NumericVector xPredJ = xPred[k];
        arma::vec xPredJArma(xPredJ.begin(),xPredJ.size());
        NumericMatrix BxPredJ = BxPred[k];
        arma::mat BxPredJArma(BxPredJ.begin(),BxPredJ.nrow(),BxPredJ.ncol(),false);
        arma::vec alphaJ = alphaPredUpdate.col(k);
        arma::vec betaJ = betaPredUpdate.col(k);
        List xBasisPredind = xBasisPred_ind[k];
        
        int xn = xPredJ.size();
        arma::vec pxPred(xn);
        for(int l=0; l<xn; l++){
          if(xPredJArma(l) > 0){
            pxPred(l) = 1.0;
          } else {
            pxPred(l) = 0.0;
          }
        }
        // // Rcout << "xn : " << xn << "\n";
        // // Rcout << "dJnew : " << dJnew << "\n";
        arma::vec dJold = dPredPrev.col(k);
        arma::vec cUpdateJ = cPredPrev.col(k);
        // if(i>1){
        //   dJold = dPredUpdate.col(k);
        //   cUpdateJ = cPredUpdate.col(k);
        // } else {
        //   double sam = sampleFromRange(nData);
        //   dJold = dUpdate.col(sam);
        //   cUpdateJ = cUpdate.col(sam);
        // }
        
        arma::vec dJnew = dJold + sampleFromNormal(dJold.size(),0.0,sigmaDTune);
        double xiXJ = xiXUpdate;
        
        
        // Compute scale parameters
        arma::vec dxsigmaold = arma::exp(BxPredJArma * dJold);
        arma::vec dxsigmanew = arma::exp(BxPredJArma * dJnew);
        
        // Example of computing cMuOld and cMuNew directly in Armadillo
        arma::vec cMuOld = alphaJ + betaJ % dJold;
        arma::vec cMuNew = alphaJ + betaJ % dJnew;
        
        // Obtain the prior of c for the old and the new values
        arma::vec cPrioriOldd(m);
        arma::vec cPrioriNewd(m);
        for(int l=0;l<m;l++){
          cPrioriOldd(l) = log(normalPDF(cUpdateJ(l), cMuOld(l), 1/sigmaCPrecUpdate));
          cPrioriNewd(l) = log(normalPDF(cUpdateJ(l), cMuNew(l), 1/sigmaCPrecUpdate));
        }
        
        // Obtain the prior of d for the old and the new values
        double dPrioriOldd = mvnPDF(dJold, muDArma, SigmaDArma);
        double dPrioriNewd = mvnPDF(dJnew, muDArma, SigmaDArma);
        // int xnPred = xPredJ.size();
        
        // Obtain likelihood
        arma::vec xlikOld(m);
        arma::vec xlikNew(m);
        for(int l=0;l<m;l++){
          // Define the vectors for iteration
          NumericVector xBasisPredindVec = xBasisPredind[l];
          arma::vec xBasisPredArma(xBasisPredindVec.begin(),xBasisPredindVec.size());
          arma::vec tmp(xBasisPredArma.n_elem);
          arma::vec tmp2(xBasisPredArma.n_elem);
          for(int t=0; t<xBasisPredArma.n_elem;t++){
            int tmp_int = xBasisPredArma(t)-1;
            tmp(t) = diracgpdPDF(xPredJArma(tmp_int), xiXJ, dxsigmaold(tmp_int), pxPred(tmp_int));
            tmp2(t) = diracgpdPDF(xPredJArma(tmp_int), xiXJ, dxsigmanew(tmp_int), pxPred(tmp_int));
          }
          xlikOld(l) = sum(log(tmp));
          xlikNew(l) = sum(log(tmp2));
        }
        
        // Acceptance
        arma::vec pDold = dPrioriOldd + xlikOld + cPrioriOldd;
        arma::vec pDnew = dPrioriNewd + xlikNew + cPrioriNewd;
        arma::vec pDaccept = pDnew - pDold;
        for(int l=0;l<m;l++){
          if (pDaccept(l) >= std::log(R::runif(0.0, 1.0))) {
            dPredUpdate(l,k) = dJnew(l);
          } else {
            dPredUpdate(l,k) = dJold(l);
          }
        }
        
        // Maybe I should replicate this here too
        // Obtain likelihood
        // Acceptance
        // double pDold = dPrioriOldd + xlikOldi + sum(cPrioriOldd);
        // double pDnew = dPrioriNewd + xlikNewi + sum(cPrioriNewd);
        // double pDaccept = pDnew - pDold;
        // if (pDaccept >= std::log(R::runif(0.0, 1.0))) {
        //   dPredUpdate.col(k) = dJnew;
        // } else {
        //   dPredUpdate.col(k) = dJold;
        // }
      }
      // dPredOld = dPredUpdate;
      
      // Update C Predict
      for(int k=0; k<m; k++){
        for(int l=0; l<nPred; l++){
          cPredUpdate(k,l) = R::rnorm(alphaPredUpdate(k,l)+betaPredUpdate(k,l)*dPredUpdate(k,l),sqrt(1/sigmaCPrecUpdate));
        }
      }
      // cPredOld = cPredUpdate;
      
      // Predict (random samples) of p
      for(int k=0; k<nPred; k++){
        NumericVector xPredJ = xPred[k];
        arma::vec xPredJArma(xPredJ.begin(),xPredJ.size());
        NumericMatrix xcovPredJ = xPredData_cov[k];
        arma::mat xcovPredJArma(xcovPredJ.begin(),xcovPredJ.nrow(),xcovPredJ.ncol(),false);
        arma::vec lams(nlam);
        for(int l=0; l<4; l++){
          lams(l) = lambdaUpdate(sampleFromRange(nData),l);
        }
        
        lambdaPredUpdate.row(k) = lams.t();
        
        arma::vec tmp(r);
        for(int l=0;l<r;l++){
          tmp(l) = invlogit(sumWithoutNaN((xcovPredJArma.row(l) % lambdaPredUpdate.row(k)).t()));
        }
        p_PredUpdate.row(k) = tmp.t();
      }
      
      
      arma::vec alphaPredUpdateVec = arma::vectorise(alphaPredUpdate);
      arma::vec betaPredUpdateVec = arma::vectorise(betaPredUpdate);
      arma::vec dPredUpdateVec = arma::vectorise(dPredUpdate);
      arma::vec cPredUpdateVec = arma::vectorise(cPredUpdate);
      //
      arma::vec lambdaPredUpdateVec = arma::vectorise(lambdaPredUpdate.t());
      //
      
      
      subchainOutputArma(j,arma::span(5+4*m*nData+4*nData,4+4*m*nData+m*nPred + 4*nData)) = alphaPredUpdateVec.t();
      subchainOutputArma(j,arma::span(5+4*m*nData+m*nPred + 4*nData,4+4*m*nData+2*m*nPred+4*nData)) = betaPredUpdateVec.t();
      subchainOutputArma(j,arma::span(5+4*m*nData+2*m*nPred + 4*nData,4+4*m*nData+3*m*nPred+4*nData)) = dPredUpdateVec.t();
      subchainOutputArma(j,arma::span(5+4*m*nData+3*m*nPred + 4*nData,4+4*m*nData+4*m*nPred+4*nData)) = cPredUpdateVec.t();
      //
      subchainOutputArma(j,5+4*m*nData+4*m*nPred+4*nData) = xiXUpdate;
      subchainOutputArma(j,6+4*m*nData+4*m*nPred+4*nData) = xiYUpdate;
      subchainOutputArma(j,arma::span(7+4*m*nData+4*m*nPred + 4*nData,6+4*m*nData+4*m*nPred+4*nData + 4*nPred)) = lambdaPredUpdateVec.t();
      // 
      
      // Predict Y
      arma::mat yPredUpdate(r,nPred);
      for(int k=0; k<nPred; k++){
        arma::vec cJ = cPredUpdate.col(k);
        arma::vec SigmaYPred = arma::exp(ByPredArma * cJ);
        arma::vec ps = p_PredUpdate.row(k).t();
        arma::vec yTemp(r);
        for(int l=0; l<r; l++){
          yTemp(l) = sample_diracDelta(SigmaYPred(l),xiYUpdate,ps(l));
        }
        yPredUpdate.col(k) = yTemp;
      }
      arma::vec yPredUpdateVec = arma::vectorise(yPredUpdate);
      subchainOutputArma(j,arma::span(7+4*m*nData+4*m*nPred+4*nData + 4*nPred,6+4*m*nData+4*m*nPred+4*nData+4*nPred+r*nPred)) = yPredUpdateVec.t();
    }
    outputArma.row(i) = subchainOutputArma.row(nThin);
  }
  NumericMatrix output = wrap(outputArma);
  return output;
}

