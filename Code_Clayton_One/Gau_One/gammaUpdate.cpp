#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

//' @keywords internal
// [[Rcpp::export]]
List gammaUpdate(const Rcpp::List& b_, const Rcpp::List& z_, const Rcpp::List& w_,
                 const Rcpp::List& pb_, const arma::vec& haz, const Rcpp::List& v_,
                 const Rcpp::List& u_,const Rcpp::List& h_, const int& K, const int& q,
                 const int& qq, const int& nev, const arma::vec& jcount) {
   
   // Newton-Raphson updates of gamma (E-step and M-step) using an exact observed
   // information calculation
   
   // declare score, E[delta x v*], and information matrix
   arma::mat Si = arma::zeros<arma::mat>(q+qq+K, w_.size()); //Score vector
   arma::vec S = arma::zeros<arma::vec>(q+qq+K);
   arma::mat Evstari = arma::zeros<arma::mat>(q+qq+K, w_.size());
   arma::mat I = arma::zeros<arma::mat>(q+qq+K, q+qq+K); //Information matrix
   arma::mat Gammaj = arma::zeros<arma::mat>(q+qq+K, nev);
   arma::mat Gamma = arma::zeros<arma::mat>(q+qq+K, q+qq+K);
   
   // loop over subjects
   for (int i=0; i<w_.size(); i++) {
      
      Rcpp::checkUserInterrupt();
      
      // extract matrices from lists for subject i
      arma::mat b = Rcpp::as<arma::mat>(b_[i]);
      arma::mat z = Rcpp::as<arma::mat>(z_[i]);
      arma::mat w = Rcpp::as<arma::mat>(w_[i]);
      arma::vec pb = Rcpp::as<arma::vec>(pb_[i]);
      arma::vec v = Rcpp::as<arma::vec>(v_[i]);
      arma::vec u = Rcpp::as<arma::vec>(u_[i]); //this is a column vector
      Rcpp::DataFrame h = Rcpp::as<Rcpp::DataFrame>(h_[i]);
      
      // subjects who are censored before the first failure time
      // do not contribute towards \gamma estimation
      int tj_ind = h["tj.ind"];
      if (tj_ind == 0) continue;
      
      int nj = w.n_cols;  // number of failure times upto time T_i
      int delta = h["delta"]; // delta_i
      
      arma::mat Ii_int = arma::zeros<arma::mat>(q+qq+K, q+qq+K); // information matrix (uninitialized) for subject i
      arma::mat bzt = b * z; // b x t(z)
      arma::mat bztev = bzt % repmat(w, 1, K); 
      // b x t(Z) . exp{W_1 %*% gamma_{11} + Z^*_11(t) %*% gamma_{12} + Z^*_11(t) %*% b * gamma_{12}} 
      arma::mat uev = repmat(trans(u), b.n_rows, 1) % repmat(w, 1, K); //1*nj
      // Z^*_{11}(t) .* exp{W_1 %*% gamma_{11} + Z^*_11(t) %*% gamma_{12} + Z^*_11(t) %*% b * gamma_{12}}
      arma::mat Eexpvj = (mean(w.each_col() % pb, 0)) % trans(haz.subvec(0, nj-1)); // this is a row vector
      arma::mat uEexpvj = Eexpvj % trans(u); //E lam_l(t) exp{} Z^*_11(t), this is row vector
      arma::mat u2Eexpvj = Eexpvj % trans(u) % trans(u); //E lam_l(t) exp{} Z^*_11(t) ^ 2
      arma::mat Eexpv = sum(Eexpvj, 1); // lambda0 x E[exp(v*gamma)]
      arma::mat uEexpv = sum(uEexpvj, 1); //sum(E lam_l(t) exp{} Z^*_11(t))
      arma::mat u2Eexpv = sum(u2Eexpvj, 1); //sum(E lam_l(t) exp{} Z^*_11(t) ^ 2)
      arma::mat hexpand = trans(repmat(haz.subvec(0, nj-1), K, 1)); // K reps of lambda0(tj)
      arma::mat outj = (mean(bztev.each_col() % pb, 0)) % hexpand;
      arma::mat uEbexpbj = (mean(bztev.each_col() % pb, 0)) % trans(u) % hexpand; //cross-prod
      arma::mat bzt2ev = bzt % bztev; // [b x t(z)]^2 . exp(v*gamma)
      arma::mat Ii_int_Kdiag = (mean(bzt2ev.each_col() % pb, 0)) % hexpand;
      
      arma::mat Eb = mean(b.each_col() % pb, 0);
      
      // loop of K longitudinal outcomes
      for (int k=0; k<K; k++) {
         // E[delta x v*(T_i)]
         Evstari(q+qq+k, i) = delta * arma::dot(z.col(nj*(k+1)-1), Eb);
         
         // score elements for K Zb's
         Si(q+qq+k, i) = arma::as_scalar(sum(outj.cols(nj*k, nj*(k+1)-1), 1));
         
         // cross-prod (diagonal) elements for K Zb's only
         Ii_int(q+qq+k, q+qq+k) = arma::as_scalar(sum(Ii_int_Kdiag.cols(nj*k, nj*(k+1)-1), 1));
         
         // cross-prod (off-diagonal) elements for K Zb's only
         for (int k2=k+1; k2<K; k2++) {
            arma::mat bztcross = bztev.cols(nj*k, nj*(k+1)-1) % bzt.cols(nj*k2, nj*(k2+1)-1);
            Ii_int(q+qq+k, k2+qq+q) = arma::as_scalar(sum((mean(bztcross.each_col() % pb, 0)) % trans(haz.subvec(0, nj-1)), 1));
            Ii_int(k2+qq+q, q+qq+k) = Ii_int(q+qq+k, k2+qq+q);
         }
         
         // cross-prod elements for q V_i's and K Zb's
         if (q > 0) {
            for (int j=0; j<q; j++) {
               Ii_int(j, q+qq+k) = Si(q+qq+k, i) * v(j);
               Ii_int(q+qq+k, j) = Ii_int(j, q+qq+k);
            }
            for(int jj = q; jj < (q+qq); jj++)
            {
               Ii_int(jj, q+qq+k) = arma::as_scalar(sum(uEbexpbj.cols(nj*k, nj*(k+1)-1), 1));//qq=1
               Ii_int(q+qq+k, jj) = Ii_int(jj, q+qq+k);
            }
            
         }
      } // end loop over outcomes k
      
      // Gamma_j vectors (for cross-product later)
      Gammaj.submat(q+qq, 0, q+qq+K-1, nj-1) += delta * trans(reshape(trans(outj), nj, K));
      
      // elements for q V_i's (only when q > 0)
      if (q > 0) {
         // E[delta x v]
         Evstari.submat(0, i, q-1, i) = delta * v;
         
         Evstari.submat(q, i, q+qq-1, i) = delta * u.row(nj-1); //delta_i * Z^*_{11}(u_i)
         // score elements
         Si.submat(0, i, q-1, i) = v * Eexpv;  // Esum(exp{}W_1)
         Si.submat(q, i, q+qq-1, i) = uEexpv; //Esum(exp{}Z^*_{11}(t))
         // cross-prod elements
         Ii_int.submat(0, 0, q-1, q-1) = (v * trans(v)) * arma::as_scalar(Eexpv);
         for(int jj = 0; jj < q; jj++)
         {
            Ii_int.submat(jj, q, jj, q+qq-1) = uEexpv * v(jj);
         }
         Ii_int.submat(q, 0, q+qq-1, q-1) = Ii_int.submat(0, q, q-1, q+qq-1);
         Ii_int.submat(q, q, q+qq-1, q+qq-1) = u2Eexpv;
         
         // Gamma_j elements
         Gammaj.submat(0, 0, q-1, nj-1) += delta * v * Eexpvj;
         //Gammaj.submat(q, 0, q+qq-1, nj-1) += uEbexpbj;//有问题
         Gammaj.submat(q, 0, q+qq-1, nj-1) += delta * uEexpvj;
      }
      
      S += (Evstari.col(i) - delta * Si.col(i)); // NB: actual score is sum(Evstari - Si)
      I += delta * Ii_int;
      
   } // end loop over subjects i
   
   // lambda0 x Gamma_j sum term (minus from information matrix)
   for (int t=0; t<nev; t++) {
      //Rcpp::Rcout << jcount(t) << std::endl;
      Gamma += Gammaj.col(t) * trans(Gammaj.col(t)) / jcount(t);
   }
   
   return List::create(
      Named("gDelta")  = solve(I, S),
      Named("I")  = I,
      Named("S")  = S
   );
   
}


//' @keywords internal
// [[Rcpp::export]]
arma::mat hazHat(const Rcpp::List& w_, const Rcpp::List& pb_, const arma::vec& nev) {
   
   // lambda0(t) for profile score function of beta
   
   arma::vec haz = arma::zeros<arma::vec>(nev.n_elem);
   
   // loop over subjects
   for (int i=0; i<w_.size(); i++) {
      
      // extract matrices from lists for subject i
      arma::mat w = Rcpp::as<arma::mat>(w_[i]);
      arma::vec pb = Rcpp::as<arma::vec>(pb_[i]);
      
      haz.subvec(0, w.n_cols-1) += arma::trans(mean(w.each_col() % pb, 0));
      
   }
   
   return(nev / haz);
   
}

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

//' @keywords internal
 // [[Rcpp::export]]
 arma::mat lambdaUpdate(const Rcpp::List& b_, const Rcpp::List& imat_,
                        const Rcpp::List& zt_, const Rcpp::List& pb_,
                        const Rcpp::List& v_, const Rcpp::List& u_,
                        const arma::mat& gam, const arma::vec& gam_vec, 
                        const int& q, const arma::vec& nev, const Rcpp::List& h_) {
   
   // Updates of lambda0 (E-step and M-step)
   
   arma::vec haz = arma::zeros<arma::vec>(nev.n_elem);
   
   // loop over subjects
   for (int i=0; i<b_.size(); i++) {
     
     // extract matrices from lists for subject i
     arma::mat b = Rcpp::as<arma::mat>(b_[i]);
     arma::mat I = Rcpp::as<arma::mat>(imat_[i]);
     arma::mat zt = Rcpp::as<arma::mat>(zt_[i]);
     arma::vec pb = Rcpp::as<arma::vec>(pb_[i]);
     arma::vec v = Rcpp::as<arma::vec>(v_[i]);
     arma::vec u = Rcpp::as<arma::vec>(u_[i]);
     Rcpp::DataFrame h = Rcpp::as<Rcpp::DataFrame>(h_[i]);
     
     // subjects who are censored before the first failure time
     // do not contribute towards \lambda estimation
     int tj_ind = h["tj.ind"];
     if (tj_ind == 0) continue;
     
     arma::mat expW_new = exp((b * gam) * trans(I * zt)) % repmat(trans(u), b.n_rows, 1); //nMC*tj_ind
     arma::mat EexpVstar = mean(expW_new.each_col() % pb, 0);
     if (q > 0) {
       EexpVstar *= arma::as_scalar(exp(v.t() * gam_vec.subvec(0, q-1)));
     }
     haz.subvec(0, EexpVstar.n_cols-1) += EexpVstar.t();
     
   }
   
   return(nev / haz);
   
 }
 
 
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

// [[Rcpp::export]]
List gammaUpdate_fixDelta(const Rcpp::List& b_, const Rcpp::List& z_, const Rcpp::List& w_,
                          const Rcpp::List& pb_, const arma::vec& haz,
                          const Rcpp::List& v_, const Rcpp::List& u_,
                          const Rcpp::List& h_,
                          const int& K, const int& q, const int& qq,
                          const int& nev, const arma::vec& jcount) {

  const int G = q + qq + K;
  const int nsubj = w_.size();

  arma::mat Si      = arma::zeros<arma::mat>(G, nsubj);
  arma::mat Evstari = arma::zeros<arma::mat>(G, nsubj);
  arma::vec S       = arma::zeros<arma::vec>(G);
  arma::mat I       = arma::zeros<arma::mat>(G, G);

  arma::mat Gammaj  = arma::zeros<arma::mat>(G, nev);
  arma::mat Gamma   = arma::zeros<arma::mat>(G, G);

  for (int i = 0; i < nsubj; i++) {

    Rcpp::checkUserInterrupt();

    arma::mat b   = Rcpp::as<arma::mat>(b_[i]);     // nMC x r
    arma::mat z   = Rcpp::as<arma::mat>(z_[i]);
    arma::mat w   = Rcpp::as<arma::mat>(w_[i]);     // nMC x nj
    arma::vec pb  = Rcpp::as<arma::vec>(pb_[i]);    // nMC
    arma::vec v   = Rcpp::as<arma::vec>(v_[i]);     // q
    arma::mat u   = Rcpp::as<arma::mat>(u_[i]);     // SHOULD be nj x qq, but may come in as 1 x nj
    Rcpp::DataFrame h = Rcpp::as<Rcpp::DataFrame>(h_[i]);

    int tj_ind = h["tj.ind"];
    if (tj_ind == 0) continue;

    int nj    = w.n_cols;
    int delta = h["delta"];

    // ---- critical shape fix for u ----
    // If R passes Zt.1 (1 x nj), convert to (nj x 1).
    if (u.n_rows == 1 && u.n_cols == nj) {
      u = u.t();
    }
    // If qq==1 but u is a vector-like mat, enforce nj x 1
    if (qq == 1 && u.n_cols != 1 && u.n_rows == nj) {
      // e.g., u accidentally nj x something -> take first col
      u = u.col(0);
    }
    // sanity: require u has nj rows
    if ((int)u.n_rows != nj) {
      stop("u_ has incompatible shape: expected nj rows.");
    }
    if ((int)u.n_cols != qq) {
      stop("u_ has incompatible shape: expected qq columns.");
    }
    // ----------------------------------

    arma::mat Ii_int = arma::zeros<arma::mat>(G, G);

    arma::mat bzt   = b * z;
    arma::mat bztev = bzt % repmat(w, 1, K);

    arma::rowvec Eexpvj = mean(w.each_col() % pb, 0) % trans(haz.subvec(0, nj-1)); // 1 x nj
    double Eexpv = arma::as_scalar(sum(Eexpvj));

    arma::rowvec uEexpvj_row(qq, arma::fill::zeros);
    arma::mat    u2Eexpv_mat(qq, qq, arma::fill::zeros);
    for (int l = 0; l < nj; l++) {
      arma::rowvec ul = u.row(l);      // 1 x qq
      double el = Eexpvj(l);
      uEexpvj_row += el * ul;
      u2Eexpv_mat += el * (ul.t() * ul);
    }

    arma::mat hexpand = trans(repmat(haz.subvec(0, nj-1), K, 1)); // 1 x (nj*K)
    arma::mat outj = (mean(bztev.each_col() % pb, 0)) % hexpand;

    arma::mat bzt2ev = bzt % bztev;
    arma::mat Ii_int_Kdiag = (mean(bzt2ev.each_col() % pb, 0)) % hexpand;

    arma::rowvec Eb = mean(b.each_col() % pb, 0);

    // ----- K blocks -----
    for (int k = 0; k < K; k++) {
      Evstari(q + qq + k, i) = delta * arma::dot(z.col(nj*(k+1)-1), Eb);
      Si(q + qq + k, i) = arma::as_scalar(sum(outj.cols(nj*k, nj*(k+1)-1), 1));
      Ii_int(q + qq + k, q + qq + k) =
        arma::as_scalar(sum(Ii_int_Kdiag.cols(nj*k, nj*(k+1)-1), 1));

      for (int k2 = k + 1; k2 < K; k2++) {
        arma::mat bztcross =
          bztev.cols(nj*k, nj*(k+1)-1) % bzt.cols(nj*k2, nj*(k2+1)-1);
        Ii_int(q + qq + k, q + qq + k2) =
          arma::as_scalar(sum((mean(bztcross.each_col() % pb, 0)) % trans(haz.subvec(0, nj-1)), 1));
        Ii_int(q + qq + k2, q + qq + k) = Ii_int(q + qq + k, q + qq + k2);
      }

      if (q > 0) {
        for (int j = 0; j < q; j++) {
          Ii_int(j, q + qq + k) = Si(q + qq + k, i) * v(j);
          Ii_int(q + qq + k, j) = Ii_int(j, q + qq + k);
        }
      }
      for (int jj = 0; jj < qq; jj++) {
        double tmp = 0.0;
        for (int l = 0; l < nj; l++) {
          tmp += haz(l) * arma::as_scalar(mean(bztev.col(nj*k + l) % pb)) * u(l, jj);
        }
        Ii_int(q + jj, q + qq + k) = tmp;
        Ii_int(q + qq + k, q + jj) = tmp;
      }
    }

    // Gammaj for K block (NO delta)
    Gammaj.submat(q + qq, 0, q + qq + K - 1, nj - 1) += trans(reshape(trans(outj), nj, K));

    // ----- v and u blocks -----
    if (q > 0) {
      Evstari.submat(0, i, q - 1, i) = delta * v;
      Evstari.submat(q, i, q + qq - 1, i) = delta * trans(u.row(nj - 1));

      Si.submat(0, i, q - 1, i) = v * Eexpv;
      Si.submat(q, i, q + qq - 1, i) = trans(uEexpvj_row);

      Ii_int.submat(0, 0, q - 1, q - 1) = (v * trans(v)) * Eexpv;

      for (int jj = 0; jj < q; jj++) {
        Ii_int.submat(jj, q, jj, q + qq - 1) = trans(uEexpvj_row) * v(jj);
      }
      Ii_int.submat(q, 0, q + qq - 1, q - 1) = Ii_int.submat(0, q, q - 1, q + qq - 1);
      Ii_int.submat(q, q, q + qq - 1, q + qq - 1) = u2Eexpv_mat;

      // Gammaj for v and u blocks (NO delta)
      for (int l = 0; l < nj; l++) {
        Gammaj.submat(0, l, q - 1, l) += v * Eexpvj(l);
        Gammaj.submat(q, l, q + qq - 1, l) += trans(u.row(l)) * Eexpvj(l);
      }
    }

    // ---- critical delta placement fix ----
    S += (Evstari.col(i) - Si.col(i));
    I += Ii_int;
  }

  for (int t = 0; t < nev; t++) {
    Gamma += Gammaj.col(t) * trans(Gammaj.col(t)) / jcount(t);
  }

  arma::vec gDelta = solve(I - Gamma, S);

  return List::create(
    Named("gDelta") = gDelta,
    Named("I")      = I,
    Named("Gamma")  = Gamma,
    Named("S")      = S
  );
}