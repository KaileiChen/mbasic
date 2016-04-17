#include "mcmc.h"
#include <math.h>

//this file is under testing

SEXP madbayes_zeroInflated(SEXP _clusterLabels, SEXP _Theta, SEXP _Mu, SEXP _D,
	      SEXP _Gamma, SEXP _Y, SEXP _lambdaw, SEXP _lambda,
	      SEXP _maxitr, SEXP _tol, SEXP _inflation, SEXP _p0, SEXP _p1) {
	
	// The following values are 1updated in MCMC iterations
	IntegerVector clusterLabels(_clusterLabels); // length I
	IntegerMatrix Theta(_Theta); // K by I
	NumericMatrix Mu(_Mu); // N by S
	IntegerVector D(_D); // Length N, valued in {0, 1, ..., K-1}
	double lambda = as<double>(_lambda);

	// The following values are piror parameters and are fixed
	double lambdaw = as<double>(_lambdaw);

	double tol = as<double>(_tol);
	int maxitr = as<int>(_maxitr);

 //      Rcout << "This line is ok!" << std::endl;

	// extract the dimensions
	int I = Theta.ncol();
	int S = Mu.ncol();
	int K = Theta.nrow();
	int N = D.size();

//    Rcout << "I = " << I << std::endl;

	// Distance from each unit to each cluster
	NumericMatrix Dist(I, I+1);
//	Rcout << "This line is ok3!" << std::endl;
	IntegerVector clusterSizes(I + 1);
//	Rcout << "This line is ok4!" << std::endl;

	IntegerMatrix inflation(_inflation); // N by I

/*	for (int n=0; n < N; n++){
		printf("%d\t",inflation(n,0));
	}
	printf("\n");

	for (int n=0; n < N; n++){
		printf("%d\t",inflation(n,I-1));
	}
	printf("\n");

			Rcout << "This line is ok!" << std::endl;
*/
	NumericVector p0(_p0); // Length N, valued in [0,1]
    NumericVector p1(_p1); // Length N, valued in [0,1]

/*	for (int n=0; n < N; n++){
		printf("%3.3f\t", p0[n]);
	}
	printf("\n");

	for (int n=0; n < N; n++){
		printf("%3.3f\t",p1[n]);
	}
	printf("\n");
*/
//		Rcout << "This line is ok!" << std::endl;
	

	// The following is the external information.
	NumericMatrix Gamma(_Gamma); // N*S by I
	NumericMatrix Y(_Y); // N by I

//       	Rcout << "This line is ok1!" << std::endl;

	// The following will be computed
	NumericMatrix W(K * S, I + 1);
	NumericMatrix P(I, S);
	NumericMatrix Sigma(N, S);

//	Rcout << "This line is ok2!" << std::endl;

	// iterators
	int i, j, k = 0, s = 0, n, i1, j1, itr;//, likid;
	int firstLabel, lastLabel;
	double loss = 0, oldloss;

	double _LOW = 1e-10;

//	Rcout << "This line 2 is ok!" << std::endl;

	for(i = 0; i < I + 1; i ++) {
		clusterSizes[i] = 0;
		/*	if (i % 100 == 0) {
		  Rcout << "i:" << i << std::endl;
		  }*/
	}

	//Rcout << "This line 2 is ok!" << std::endl;

	// Compute J
	int J = 0;
	int nextClusterLabel = 0;	
	for(i = 0; i < I; i ++) {
		if(J < clusterLabels[i]) {
			J = clusterLabels[i];
		}
		clusterSizes[clusterLabels[i]] ++;
	}
	J ++;
	nextClusterLabel = J;
	printf("start for lambda %3.3f\n", lambda / lambdaw);


	IntegerVector zero(N); // Length N, number of zeros in each data set
	IntegerVector one(N); // Length N, number of ones in each data set
	IntegerVector two(N); // Length N, number of ones in each data set
	
	for (n=0; n < N; n++){
		zero[n] = 0;
		one[n] = 0;
		for (i=0; i < I; i++){
			if (inflation(n, i) == 1){
				zero[n]++;
			} else if (inflation(n, i) == 2){
				one[n]++;
			} else if (inflation(n, i) == 4){
				two[n]++;
			} 
		}
	}

	for(itr = 0; itr < maxitr; itr ++) {
		oldloss = loss;
		// Update W
		NumericVector counts(S);
		NumericVector w_sol(S);
		for(j = 0; j < J; j ++) {
			if(clusterSizes[j] > 0) {
				for(k = 0; k < K; k ++) {
					for(s = 0; s < S; s ++) {
						counts[s] = 0;
					}
					for(i = 0; i < I; i ++) {
						if(clusterLabels[i] == j) {
							for(s = 0; s < S; s ++) {
								if(Theta(k, i) == s) {
									counts[s] ++;
								}
							}
						}
					}
					for(s = 0; s < S; s ++) {
						w_sol[s] = (counts[s] + _LOW) / ((double) clusterSizes[j] + S * _LOW);
					}
					for(s = 0; s < S; s ++) {
						W(s * K + k, j) = w_sol(s);
					}
				}
			}
		}

//		loss = ComputeLoss_zeroInflated(D, Theta, Y, Mu, Gamma, W, clusterLabels, clusterSizes, lambda, lambdaw, inflation, p0, p1);
//		printf("update W, Loss function = %3.3f, number of clusters = %d\n", loss, J);

		// update Theta
		for(i = 0; i < I; i ++) {
			for(k = 0; k < K; k ++) {
				double tmp[S];
				// Omit the Ones and Zeros in non enriched state
				s = 0;
				tmp[s] = 0.0;
				for(n = 0; n < N; n ++) {
					if(D[n] == k) {
						if (inflation(n , i) != 1 && inflation(n , i) != 2){
							tmp[s] += (Y(n, i) - Mu(n, s) * Gamma(s * N + n, i)) * (Y(n, i) - Mu(n, s) * Gamma(s * N + n, i));
						}
					}
				}

				// use all loci in the other states
				for(s = 1; s < S; s ++) {
					// initialize
					tmp[s] = 0.0;
					for(n = 0; n < N; n ++) {
						if(D[n] == k) {
							tmp[s] += (Y(n, i) - Mu(n, s) * Gamma(s * N + n, i)) * (Y(n, i) - Mu(n, s) * Gamma(s * N + n, i));
						}
					}
				}

				// Assign new values
				Theta(k, i) = 0;
				for(s = 1; s < S; s ++) {
					if(tmp[s] < tmp[Theta(k, i)])
						Theta(k, i) = s;
				}
			}
		}

//		loss = ComputeLoss_zeroInflated(D, Theta, Y, Mu, Gamma, W, clusterLabels, clusterSizes, lambda, lambdaw, inflation, p0, p1);
//		printf("update states, Loss function = %3.3f, number of clusters = %d\n", loss, J);


		
//		double loss1 = ComputeLoss_zeroInflated(D, Theta, Y, Mu, Gamma, W, clusterLabels, clusterSizes, lambda, lambdaw, inflation, p0, p1);
//			printf("update mu, Loss function = %3.3f\n", loss1);

		
	       	// update p0, p1
			NumericVector pnew0(N);
			NumericVector pnew1(N);
		for (n = 0; n < N; n ++) {
			double poisson_mean = Mu(n, 0) * Mu(n, 0) + Sigma(n, 0) * Sigma(n, 0) - 0.25;
			poisson_mean = fmax(2.0, poisson_mean);
			poisson_mean = fmin(8.0, poisson_mean);
			int NonState = 0;
			for (i = 0; i < I; i ++) {
				if (Theta(D[n], i) == 0) {
					NonState++;
				}
			}
			double pStar = (double)(two[n] + 1.0) / (double)(NonState+ 10.0);
			pStar = (pStar + pStar) * exp(poisson_mean) / (poisson_mean * poisson_mean);
			pStar = fmin(pStar, 0.99);
			
			pnew0[n] = (double)(zero[n] + pnew0[n]) / (double)(NonState + 1.0) - pStar * exp(-poisson_mean);
			pnew1[n] = (double)(one[n] + pnew1[n]) / (double)(NonState + 1.0) - pStar * poisson_mean * exp(-poisson_mean);
			double temp = pStar + pnew0[n] + pnew1[n];
			pnew0[n] /= temp;
			pnew1[n] /= temp;
			pnew0[n] = fmax(pnew0[n], 0.05);
			pnew1[n] = fmax(pnew1[n], 0.05);
		}
			
		p0 = 0.9*p0 + 0.1 * pnew0;
		p1 = 0.9* p0 + 0.1 * pnew1;

				
		//update clusters
		for(i = 0; i < I; i++) {
			int oldClusterLabel = clusterLabels[i];
			clusterSizes[clusterLabels[i]] --;
			double mintmp = lambda;// cost of starting a new cluster
			int newClusterLabel = nextClusterLabel;
			for(j = 0; j < J; j ++) {
				if(clusterSizes[j] > 0) {
					Dist(i, j) = 0;
					for(k = 0; k < K; k ++) {
						for(s = 0; s < S; s ++) {
							if(Theta(k, i) == s) {
								Dist(i, j) += (1 - W(s * K + k, j)) * (1 - W(s * K + k, j));
							} else {
								Dist(i, j) += W(s * K + k, j) * W(s * K + k, j);
							}
						}
					}
			
					// assign the cluster label as well as the outlier label
					if(Dist(i, j) * lambdaw < mintmp) {
						mintmp = Dist(i, j) * lambdaw;
						newClusterLabel = j;
					}
				}
			}

			clusterLabels[i] = newClusterLabel;
			clusterSizes[newClusterLabel] ++;
			if(mintmp >= lambda) {
				// a new cluster is formed
				// if(J != newClusterLabel) {printf("Error: new state is not J = %d.", J);exit(1);}
				for(s = 0; s < S; s ++) {
					for(k = 0; k < K; k ++) {
						if(Theta(k, i) != s) {
							W(s * K + k, newClusterLabel) = _LOW;
						} else {
							W(s * K + k, newClusterLabel) = 1 - (S - 1) * _LOW;
						}
					}
				}
				for(i1 = 0; i1 < I; i1 ++) {
					Dist(i1, newClusterLabel) = 0;
					for(k = 0; k < K; k ++) {
						if(Theta(k, i) != Theta(k, i1)) {
							Dist(i1, newClusterLabel) += 2;
						}
					}
				}
				J ++;
			}
			if(clusterSizes[oldClusterLabel] == 0) {
				// an old cluster should be removed
				if(false) {
					W(_, oldClusterLabel) = W(_, J - 1);
					Dist(_, oldClusterLabel) = Dist(_, J - 1);
					for(k = 0; k < K * S; k ++) {
						W(k, J - 1) = 0;
						Dist(k, J - 1) = 0;
					}
					for(i1 = 0; i1 < I; i1 ++) {
						if(clusterLabels[i1] == J - 1)
							clusterLabels[i1] = oldClusterLabel;
					}
					clusterSizes[oldClusterLabel] = clusterSizes[J - 1];
					clusterSizes[J - 1] = 0;
				}
				J --;
				nextClusterLabel = oldClusterLabel;
			} else {
				// find the next cluster label
				while(clusterSizes[nextClusterLabel] > 0) {
					nextClusterLabel ++;
					if(nextClusterLabel > I) {
						nextClusterLabel -= I + 1;
					}
				}
			}
			//too many clusters
			if(J >= sqrt(I) * 3) {
				break;
			}
		}

//		loss = ComputeLoss_zeroInflated(D, Theta, Y, Mu, Gamma, W, clusterLabels, clusterSizes, lambda, lambdaw, inflation, p0, p1);
//		printf("update cluster, Loss function = %3.3f, number of clusters = %d\n", loss, J);
		// test if each cluster can be dismissed
		double dismissLoss = 0, bestLoss = 0;
		int newLabels[I];
		bool canDismiss = false;
		for(j = 0; j < I; j ++) {
			if(clusterSizes[j] > 0) {
				dismissLoss = 0;
				canDismiss = true;
				for(i = 0; i < I; i ++) {
					if(clusterLabels[i] == j) {
						// assign to the next closest cluster
						newLabels[i] = 0;
						bestLoss = 2 * K;
						for(j1 = 0; j1 < I + 1; j1 ++) {
							if(clusterSizes[j1] > 0 && j1 != j && bestLoss > Dist(i, j1)) {
								newLabels[i] = j1;
								bestLoss = Dist(i, j1);
							}
						}
						dismissLoss += bestLoss - Dist(i, j);
						if(dismissLoss * lambdaw > lambda) {
							canDismiss = false;
						}
					}
					if(! canDismiss) {
						break;
					}
				}
				if(canDismiss) {
//				 printf("Dismiss cluster %d, %3.3f > %3.3f.\n", j, lambda / lambdaw, dismissLoss);
					for(i = 0; i < I; i ++) {
						if(clusterLabels[i] == j) {
							clusterLabels[i] = newLabels[i];
							clusterSizes[newLabels[i]] ++;
						}
					}
					if(false) {
					if(j < J - 1) {
						// an old cluster should be removed
						W(_, j) = W(_, J - 1);
						Dist(_, j) = Dist(_, J - 1);
						for(i1 = 0; i1 < I; i1 ++) {
							if(clusterLabels[i1] == J - 1)
								clusterLabels[i1] = j;
						}
						clusterSizes[j] = clusterSizes[J - 1];
						j --;// need to recheck this new cluster since it is relabelled
					}
					}
					clusterSizes[j] = 0;
					if(false) {
						for(k = 0; k < K * S; k ++) {
							W(k, J - 1) = 0;
							Dist(k, J - 1) = 0;
						}
						clusterSizes[J - 1] = 0;
					}
					J --;
					nextClusterLabel = j;
				}
			}
		}

//		loss = ComputeLoss_zeroInflated(D, Theta, Y, Mu, Gamma, W, clusterLabels, clusterSizes, lambda, lambdaw, inflation, p0, p1);
//		printf("dismissed clusters, Loss function = %3.3f, number of clusters = %d\n", loss, J);

		for(n = 0; n < N; n ++) {
			double denom = 0, numer = 0;
			int n0 = I * p0[n];
			int n1 = I * p1[n];
			int count0 = 0;
			int count1 = 0;
			//Initialize Mu as the pooled mean to avoid 0/0
			for(i = 0; i < I; i ++) {
				if (inflation(n, i) == 1 && count0 < n0) {
					count0++;
				} else if (inflation(n, i) == 2 && count1 < n1) {
					count1++;
				} else {
					for(s = 0; s < S; s ++) {
						denom += Gamma(n + s * N, i);
						numer += Y(n, i);
					}
				}	
			}
			for(s = 0; s < S; s ++) {
				Mu(n, s) = numer / denom;
			}

			s = 0;
			count0 = 0;
			count1 = 0;
			denom = 0;
			numer = 0;
			for(i = 0; i < I; i ++) {
				if (inflation(n, i) == 1 && count0 < n0) {
					count0++;
				} else if (inflation(n, i) == 2 && count1 < n1) {
					count1++;
				} else {
					if(Theta(D[n], i) == s) {
						numer += Y(n, i);
						denom += Gamma(n + s * N, i);
					}
				}
			}
			if(denom > 0) {
				Mu(n, s) = numer / denom;
			}
			
		
			for(s = 1; s < S; s ++) {
				denom = 0;
				numer = 0;
				for(i = 0; i < I; i ++) {
					if(Theta(D[n], i) == s) {
						numer += Y(n, i);
						denom += Gamma(n + s * N, i);
					}
				}
				if(denom > 0) {
					Mu(n, s) = numer / denom;
				}	
			}

			// Initialize Sigma as the pooled variance to avoid 0/0
			denom = 0;
			numer = 0;
			count0 = 0;
			count1 = 0;
			for(i = 0; i < I; i ++) {
				for(s = 0; s < S; s ++) {
					if (inflation(n, i) == 1 && count0 < n0) {
						count0++;
					} else if (inflation(n, i) == 2 && count1 < n1) {
						count1++;
					} else {
						if(Theta(D[n], i) == s) {
							numer += (Y(n, i) - Gamma(n + s * N, i) * Mu(n, s)) * (Y(n, i) - Gamma(n + s * N, i) * Mu(n, s));
							denom ++;
						}
					}
				}
			}
			for(s = 0; s < S; s ++) {
				Sigma(n, s) = numer / denom;
			}

			s = 0;
			denom = 0;
			numer = 0;
			count0 = 0;
			count1 = 0;
			for(i = 0; i < I; i ++) {
				if (inflation(n, i) == 1 && count0 < n0) {
					count0++;
				} else if (inflation(n, i) == 2 && count1 < n1) {
					count1++;
				} else {
					if(Theta(D[n], i) == s) {
						numer += (Y(n, i) - Gamma(n + s * N, i) * Mu(n, s)) * (Y(n, i) - Gamma(n + s * N, i) * Mu(n, s));
						denom ++;
					}
				}
			}
			if(denom > 0) {
				Sigma(n, s) = numer / denom;
			}

			for(s = 1; s < S; s ++) {
				denom = 0;
				numer = 0;
				for(i = 0; i < I; i ++) {
					if(Theta(D[n], i) == s) {
						numer += (Y(n, i) - Gamma(n + s * N, i) * Mu(n, s)) * (Y(n, i) - Gamma(n + s * N, i) * Mu(n, s));
						denom ++;
					}
				}
				if(denom > 0) {
					Sigma(n, s) = numer / denom;
				} 
			}
		}
		

		// calculate the loss function

//		loss = ComputeLoss_zeroInflated(D, Theta, Y, Mu, Gamma, W, clusterLabels, clusterSizes, lambda, lambdaw, inflation, p0, p1);
//		printf("update distribution, Loss function = %3.3f, number of clusters = %d\n", loss, J);

		if(abs(oldloss - loss) < tol * oldloss && itr > 1) {
			break;
		}

		if(J >= sqrt(I) * 3) {
			printf("Warning: too many clusters. Consider using larger lambda.\n");
			break;
		}

		lastLabel = I;
		while(clusterSizes[lastLabel] == 0) {
			lastLabel --;
		}

		if(lastLabel > J * 2) {
			//printf("relabel clusters\n");
			firstLabel = 0, lastLabel = I;
			while(firstLabel < lastLabel) {
				while(clusterSizes[firstLabel] > 0 && firstLabel < I) {
					firstLabel ++;
				}
				while(clusterSizes[lastLabel] == 0 && lastLabel > 0) {
					lastLabel --;
				}
				if(firstLabel < lastLabel) {
					for(i = 0; i < I; i ++) {
						if(clusterLabels[i] == lastLabel) {
							clusterLabels[i] = firstLabel;
						}
					}
					W(_, firstLabel) = W(_, lastLabel);
					Dist(_, firstLabel) = Dist(_, lastLabel);
					clusterSizes[firstLabel] = clusterSizes[lastLabel];
					clusterSizes[lastLabel] = 0;
					nextClusterLabel = lastLabel;
				}
			}
		}
	}

	//printf("relabel clusters\n");
	// Reorder the clusterLabels and the columns of W
	firstLabel = 0, lastLabel = I;
	while(firstLabel < lastLabel) {
		while(clusterSizes[firstLabel] > 0 && firstLabel < I) {
			firstLabel ++;
		}
		while(clusterSizes[lastLabel] == 0 && lastLabel > 0) {
			lastLabel --;
		}
		if(firstLabel < lastLabel) {
			for(i = 0; i < I; i ++) {
				if(clusterLabels[i] == lastLabel) {
					clusterLabels[i] = firstLabel;
				}
			}
			W(_, firstLabel) = W(_, lastLabel);
			for(k = 0; k < K * S; k ++) {
				W(k, lastLabel) = 0;
			}
			clusterSizes[firstLabel] = clusterSizes[lastLabel];
			clusterSizes[lastLabel] = 0;
		}
	}

	Rcpp::List ret = Rcpp::List::create(
					    Rcpp::Named("Theta") = Theta,
					    Rcpp::Named("clusterLabels") = clusterLabels,
					    Rcpp::Named("W") = W,
					    Rcpp::Named("Mu") = Mu,
					    Rcpp::Named("Sigma") = Sigma,
					    Rcpp::Named("loss") = loss,
					    Rcpp::Named("Iter") = itr,
					    Rcpp::Named("p0") = p0,
					    Rcpp::Named("p1") = p1
					    );
	Rcout << "C++ finished again!" << std::endl;
	return( ret );
	
}

double ComputeLoss_zeroInflated(IntegerVector D, IntegerMatrix Theta, NumericMatrix Y, NumericMatrix Mu, NumericMatrix Gamma, NumericMatrix W, IntegerVector clusterLabels, IntegerVector clusterSizes,
		 double lambda, double lambdaw, IntegerMatrix inflation, NumericVector p0, NumericVector p1) {
	int I = Theta.ncol();
	int K = Theta.nrow();
	int N = D.size();
	int S = Mu.ncol();

	int i, n, k, s, J;
	double loss = 0, tmp;

	for(n = 0; n < N; n ++) {	
		int n0 = I * p0[n];
		int n1 = I * p1[n];
		int count0 = 0;
		int count1 = 0;
		for(i = 0; i < I; i ++) {
			k = D[n];
			s = Theta(k, i);
			if (s == 0 && inflation(n, i) == 1 && count0 < n0) {
				count0++;
			} else if (s == 0 && inflation(n, i) == 2 && count1 < n1) {
				count1++;
			} else {
				loss += (Y(n, i) - Mu(n, s) * Gamma(s * N + n, i)) *
					(Y(n, i) - Mu(n, s) * Gamma(s * N + n, i));
			}
		}
	}
	//	printf("loss after data %3.3f\n", loss);
		
	for(i = 0; i < I; i ++) {
		tmp = 0;
		for(k = 0; k < K; k ++) {
			for(s = 0; s < S; s ++) {
				if(Theta(k, i) == s) {
					tmp += (1 - W(s * K + k, clusterLabels[i])) * (1 - W(s * K + k, clusterLabels[i])); 
				} else {
					tmp += W(s * K + k, clusterLabels[i]) * W(s * K + k, clusterLabels[i]);
				}
			}
		}
		//		printf("i = %d, loss = %3.3f, b = %d\n", i, tmp, b[i]);
		loss += lambdaw * tmp;
	}
	//	printf("loss after w %3.3f\n", loss);
	J = 0;
	for(i = 0; i < I + 1; i ++) {
		if(clusterSizes[i] > 0) {
			J ++;
		}
	}
	loss += lambda * (J - 1);
	//	printf("loss for all terms  %3.3f\n", loss);
	return(loss);
}

double ComputeLoss_theta_zeroInflated(IntegerVector D, IntegerMatrix Theta, NumericMatrix Y, NumericMatrix Mu, NumericMatrix Gamma, IntegerMatrix inflation, NumericVector p0, NumericVector p1) {
	int I = Theta.ncol();
	int N = D.size();

	int i, n, k, s;
	double loss = 0;

	for(n = 0; n < N; n ++) {
		int n0 = I * p0[n];
		int n1 = I * p1[n];
		int count0 = 0;
		int count1 = 0;
		for(i = 0; i < I; i ++) {
			k = D[n];
			s = Theta(k, i);
			if (s == 0 && inflation(n, i) == 1 && count0 < n0) {
				count0++;
			} else if (s == 0 && inflation(n, i) == 2 && count1 < n1) {
				count1++;
			} else {
				loss += (Y(n, i) - Mu(n, s) * Gamma(s * N + n, i)) *
					(Y(n, i) - Mu(n, s) * Gamma(s * N + n, i));
			}
		}
	}
	
	return(loss);
}


SEXP madbayes_theta_zeroInflated(SEXP _Theta, SEXP _Mu, SEXP _D, SEXP _Gamma, SEXP _Y, SEXP _maxitr, SEXP _tol, SEXP _inflation, SEXP _p0, SEXP _p1) {
	
	// The following values are 1updated in MCMC iterations
	IntegerMatrix Theta(_Theta); // K by I
	NumericMatrix Mu(_Mu); // N by S
	IntegerVector D(_D); // Length N, valued in {0, 1, ..., K-1}

	// The following is the external information.
	NumericMatrix Gamma(_Gamma); // N*S by I
	NumericMatrix Y(_Y); // N by I

	int maxitr = as<int>(_maxitr);
	double tol = as<double>(_tol);

	// Ones and Zeros
	IntegerMatrix inflation(_inflation); // N by I

	NumericVector p0(_p0); // Length N, valued in [0,1]
    NumericVector p1(_p1); // Length N, valued in [0,1]

	// extract the dimensions
	int I = Theta.ncol();
	int S = Mu.ncol();
	int K = Theta.nrow();
	int N = D.size();


	// The following will be computed
	NumericMatrix Sigma(N, S);

	// iterators
	int i, k = 0, s = 0, n, itr;//, likid;
	double loss = 0, oldloss;

	IntegerVector zero(N); // Length N, number of zeros in each data set
	IntegerVector one(N); // Length N, number of ones in each data set
	IntegerVector two(N); // Length N, number of ones in each data set
	
	for (n=0; n < N; n++){
		zero[n] = 0;
		one[n] = 0;
		for (i=0; i < I; i++){
			if (inflation(n, i) == 1){
				zero[n]++;
			} else if (inflation(n, i) == 2){
				one[n]++;
			} else if (inflation(n, i) == 4){
				two[n]++;
			} 
		}
	}

//	Rcout << "n: " << n << std::endl;

	for(itr = 0; itr < maxitr; itr ++) {
		oldloss = loss;
//		Rcout << "itr: " << itr << std::endl;
		// update Theta
		for(i = 0; i < I; i ++) {
			for(k = 0; k < K; k ++) {
				double tmp[S];
				
				// Omit the Ones and Zeros in non enriched state by randomization to account for inflation
				s = 0;
				tmp[s] = 0.0;
				for(n = 0; n < N; n ++) {
					if(D[n] == k) {
						if (inflation(n , i) != 1 && inflation(n , i) != 2){
							tmp[s] += (Y(n, i) - Mu(n, s) * Gamma(s * N + n, i)) * (Y(n, i) - Mu(n, s) * Gamma(s * N + n, i));
						}
					}
				}

				// use all loci in the other states
				for(s = 1; s < S; s ++) {
					// initialize
					tmp[s] = 0.0;
					for(n = 0; n < N; n ++) {
						if(D[n] == k) {
							tmp[s] += (Y(n, i) - Mu(n, s) * Gamma(s * N + n, i)) * (Y(n, i) - Mu(n, s) * Gamma(s * N + n, i));
						}
					}
				}

				// Assign new values
				Theta(k, i) = 0;
				for(s = 1; s < S; s ++) {
					if(tmp[s] < tmp[Theta(k, i)])
						Theta(k, i) = s;
				}
			}
		}
		//Rcout << "theta updated" << std::endl;
//		loss = ComputeLoss_theta_zeroInflated(D, Theta, Y, Mu, Gamma, inflation, p0, p1);
//			printf("update states, Loss function = %3.3f\n", loss);
		
		// update mu
		for(n = 0; n < N; n ++) {
			double denom = 0, numer = 0;
			int n0 = I * p0[n];
			int n1 = I * p1[n];
			int count0 = 0;
			int count1 = 0;
			//Initialize Mu as the pooled mean to avoid 0/0
			for(i = 0; i < I; i ++) {
				if (inflation(n, i) == 1 && count0 < n0) {
					count0++;
				} else if (inflation(n, i) == 2 && count1 < n1) {
					count1++;
				} else {
					for(s = 0; s < S; s ++) {
						denom += Gamma(n + s * N, i);
						numer += Y(n, i);
					}
				}	
			}
			for(s = 0; s < S; s ++) {
				Mu(n, s) = numer / denom;
			}

			s = 0;
			count0 = 0;
			count1 = 0;
			denom = 0;
			numer = 0;
			for(i = 0; i < I; i ++) {
				if (inflation(n, i) == 1 && count0 < n0) {
					count0++;
				} else if (inflation(n, i) == 2 && count1 < n1) {
					count1++;
				} else {
					if(Theta(D[n], i) == s) {
						numer += Y(n, i);
						denom += Gamma(n + s * N, i);
					}
				}
			}
			if(denom > 0) {
				Mu(n, s) = numer / denom;
			}
			
		
			for(s = 1; s < S; s ++) {
				denom = 0;
				numer = 0;
				for(i = 0; i < I; i ++) {
					if(Theta(D[n], i) == s) {
						numer += Y(n, i);
						denom += Gamma(n + s * N, i);
					}
				}
				if(denom > 0) {
					Mu(n, s) = numer / denom;
				}	
			}

			// Initialize Sigma as the pooled variance to avoid 0/0
			denom = 0;
			numer = 0;
			count0 = 0;
			count1 = 0;
			for(i = 0; i < I; i ++) {
				for(s = 0; s < S; s ++) {
					if (inflation(n, i) == 1 && count0 < n0) {
						count0++;
					} else if (inflation(n, i) == 2 && count1 < n1) {
						count1++;
					} else {
						if(Theta(D[n], i) == s) {
							numer += (Y(n, i) - Gamma(n + s * N, i) * Mu(n, s)) * (Y(n, i) - Gamma(n + s * N, i) * Mu(n, s));
							denom ++;
						}
					}
				}
			}
			for(s = 0; s < S; s ++) {
				Sigma(n, s) = numer / denom;
			}


			s = 0;
			denom = 0;
			numer = 0;
			count0 = 0;
			count1 = 0;
			for(i = 0; i < I; i ++) {
				if (inflation(n, i) == 1 && count0 < n0) {
					count0++;
				} else if (inflation(n, i) == 2 && count1 < n1) {
					count1++;
				} else {
					if(Theta(D[n], i) == s) {
						numer += (Y(n, i) - Gamma(n + s * N, i) * Mu(n, s)) * (Y(n, i) - Gamma(n + s * N, i) * Mu(n, s));
						denom ++;
					}
				}
			}
			if(denom > 0) {
				Sigma(n, s) = numer / denom;
			}

			for(s = 1; s < S; s ++) {
				denom = 0;
				numer = 0;
				for(i = 0; i < I; i ++) {
					if(Theta(D[n], i) == s) {
						numer += (Y(n, i) - Gamma(n + s * N, i) * Mu(n, s)) * (Y(n, i) - Gamma(n + s * N, i) * Mu(n, s));
						denom ++;
					}
				}
				if(denom > 0) {
					Sigma(n, s) = numer / denom;
				}
			}
		}

//		double loss1 = ComputeLoss_theta_zeroInflated(D, Theta, Y, Mu, Gamma, inflation, p0, p1);
//			printf("update mu, Loss function = %3.3f\n", loss1);
		// update p0, p1
			NumericVector pnew0(N);
			NumericVector pnew1(N);
		for (n = 0; n < N; n ++) {
			double poisson_mean = Mu(n, 0) * Mu(n, 0) + Sigma(n, 0) * Sigma(n, 0) - 0.25;
			poisson_mean = fmax(2.0, poisson_mean);
			poisson_mean = fmin(10.0, poisson_mean);
			int NonState = 0;
			for (i = 0; i < I; i ++) {
				if (Theta(D[n], i) == 0) {
					NonState++;
				}
			}
			double pStar = (double)(two[n]) / (double)(NonState);
			pStar = (pStar + pStar) * exp(poisson_mean) / (poisson_mean * poisson_mean);
			pStar = fmin(pStar, 0.99);
			
			pnew0[n] = (double)(zero[n] + pnew0[n]) / (double)(NonState + 1.0) - pStar * exp(-poisson_mean);
			pnew1[n] = (double)(one[n] + pnew1[n]) / (double)(NonState + 1.0) - pStar * exp(-poisson_mean);
			double temp = pStar + pnew0[n] + pnew1[n];
			pnew0[n] /= temp;
			pnew1[n] /= temp;
			pnew0[n] = fmax(pnew0[n], 0.05);
			pnew1[n] = fmax(pnew1[n], 0.05);
		}
			
		// calculate the loss function
		
		loss = ComputeLoss_theta_zeroInflated(D, Theta, Y, Mu, Gamma, inflation, pnew0, pnew1);
//		printf("initialize distribution, Loss function = %3.3f\n", loss);
		
//		if (loss1 < loss){
//		  loss = loss1;
//		} else {
//		  for (n=0; n < N; n++){
//		    p0[n] = pnew0[n];
//		    p1[n] = pnew1[n];
//		  }
//		}

		p0 = 0.9*p0 + 0.1 * pnew0;
		p1 = 0.9* p0 + 0.1 * pnew1;

		if(itr > 1 && abs(oldloss - loss) < tol * oldloss) {
			break;
		}
	}

		printf("initialize distribution, Loss function = %3.3f\n", loss);
	
	Rcpp::List ret = Rcpp::List::create(
					    Rcpp::Named("Theta") = Theta,
					    Rcpp::Named("Mu") = Mu,
					    Rcpp::Named("Sigma") = Sigma,
					    Rcpp::Named("loss") = loss,
					    Rcpp::Named("Iter") = itr,
					    Rcpp::Named("p0") = p0,
					    Rcpp::Named("p1") = p1
					    );
	Rcout << "C++ finished!" << std::endl;
	return( ret );

}



