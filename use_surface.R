cppFunction(
'NumericMatrix use_surfaceC(NumericMatrix xy_wolv, NumericMatrix xy, int n_grid, NumericVector MoveP)
{	
  int pixels = xy.nrow();
  int N_wolv = xy_wolv.nrow();
  double sd_x = MoveP[0];
  double sd_y = MoveP[1];
  
  double * use  = new double[pixels];
  double * gtot = new double[n_grid];
  double tot;
	double term1, term2, term3;
  int gid;
          
  NumericMatrix out(n_grid,2);
  
  for(int wolv = 0; wolv < N_wolv; wolv++) //Builds probability of use surface for each individual.
	{
		tot = 0.0;  
    for(int g = 0; g < n_grid; g++)
      gtot[g] = 0.0;
		for(int px = 0; px < pixels; px++)
		{
			//Bivariate normal formula, broken up.
			term1 = pow((xy(px,0) - xy_wolv(wolv,0))/ sd_x, 2) + pow((xy(px,1) - xy_wolv(wolv,1))/ sd_y, 2); 
			term2 = 2* PI * sd_x * sd_y;
			term3 = log(term2);

			use[px] = exp((term1 + 2 * term3)/(-2));
			if(MoveP[3]>0)
			{
				if(term1 > pow(MoveP[2],2)) //truncates home ranges
					if(term1 < pow(MoveP[2] + MoveP[3],2)){
           use[px] = (1 - MoveP[4]) / term2;    // 1 /(2pi) = dbinorm(0,0)
        } else { use[px] = 0; }
			}
			use[px] = xy(px,2) * use[px]; // rescale by snow values
			tot += use[px]; //keep track of total probability over entire surface.
		}
		for(int px = 0; px < pixels; px++){
     gid = (int) xy(px,3);
     if( gid > 0 ){                          // Exclude non-grid areas 
       use[px] = use[px] / tot;                    //Rescale to 1 for landscape
       gtot[gid] = 1 - (1 - gtot[gid]) * (1 - use[px]);  // Accumulate cell total for ind
       } 
    }
    for(int g=1; g < n_grid; g++){
      if( gtot[g] > 0.01 ) 
          out(g,0)++;
      out(g,1) = 1 - (1- out(g,1))*(1 - gtot[g]);
    }
         
	}

  delete [] gtot;
  delete [] use;

    return out;
  }')
