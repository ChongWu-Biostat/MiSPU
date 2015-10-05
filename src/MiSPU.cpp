# include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;


// [[Rcpp::export()]]
List GUniFracCpp (arma::mat cum,arma::mat brlen) {
    
    // inputs
    const int n = cum.n_cols ;
    
    // containers
    arma::mat d0(n,n) ;
    d0.fill(0) ; // initialize d0 to 0
    
    arma::mat d5(n,n) ;
    d5.fill(0) ; // initialize d5 to 0
    
    arma::mat d1(n,n) ;
    d1.fill(0) ; // initialize d1 to 0
    
    for (int i=1; i< n; i ++) {
        for (int j=0; j<=i-1; j++) {
            arma::mat  cum11 = cum.col(i);
            arma::mat cum21 = cum.col(j);
            arma::vec cum1 = cum11.elem(arma::find( (cum11+cum21)>0));
            arma::vec cum2 = cum21.elem(arma::find( (cum11+cum21)>0));
            arma::vec brlen2 =brlen.elem(arma::find((cum11+cum21)>0 ) );
            
            arma::mat diff = abs(cum1 - cum2) / (cum1 + cum2);
            
            arma::mat  w = brlen2 ;
            d0(i,j) = accu(diff % w) /accu(w);
            d0(j,i) = d0(i,j);
            
            w = brlen2 % arma::sqrt(cum1 + cum2);
            d5(i,j) = accu(diff % w) /accu(w);
            d5(j,i) = d5(i,j);
            
            w = brlen2 %(cum1 + cum2);
            d1(i,j) = accu(diff % w) /accu(w);
            d1(j,i) = d1(i,j);
            
        }
    }
    
    // returns
    List ret ;
    ret["d0"] = d0 ;
    ret["d5"] = d5 ;
    ret["d1"] = d1 ;
    return(ret) ;
}


// [[Rcpp::export()]]
arma::cube GUniFracCpp2 (arma::mat cum,arma::mat brlen, arma::mat powV) {
    
    // inputs
    const int n = cum.n_cols ;
    const int powN = powV.n_cols;
    
    // containers
    arma::cube d(n,n,powN) ;
    d.zeros() ; // initialize distance cube to 0
    
    for (int i=1; i< n; i ++) {
        for (int j=0; j<=i-1; j++) {
            arma::mat  cum11 = cum.col(i);
            arma::mat cum21 = cum.col(j);
            arma::vec cum1 = cum11.elem(arma::find( (cum11+cum21)>0));
            arma::vec cum2 = cum21.elem(arma::find( (cum11+cum21)>0));
            arma::vec brlen2 =brlen.elem(arma::find((cum11+cum21)>0 ) );
            
            arma::mat diff = abs(cum1 - cum2) / (cum1 + cum2);
            
            // calculate the distance for different weights
            for (int powi = 0; powi<powN; powi++) {
                arma::mat weight = pow(cum1+cum2,powV(0,powi));
                arma::mat  w = brlen2 % weight;
                d(i,j,powi) = accu(diff % w) /accu(w);
                d(j,i,powi) = d(i,j,powi);
            }
        }
    }
    
    // returns
    return(d) ;
}




